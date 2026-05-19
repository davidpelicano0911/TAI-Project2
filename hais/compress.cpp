// HAIS — Hybrid Astronomical Image compressor  (version 5)
//
// Pipeline por bloco 500×500:
//   Modos: 0=raw  1=avg  2=LS(W,N,NW)  3=global_mean  4=LS+bias(W,N,NW,1)
//   Escolhe modo com menor byte_cost (H(hi8)+H(lo8))
//   Comprime com FSE/tANS (2 streams: hi8, lo8)
//
// Formato header: MAGIC(4) + W(4) + H(4) + bs(4) + gmean(2)
//
// Usage: ./compress <input> <output.hais> [width height]

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>
#include <cstring>

static constexpr uint8_t  MAGIC[4]   = {'H','A','I','S'};
static constexpr uint32_t SCALE_BITS = 12;
static constexpr uint32_t SCALE      = 1u << SCALE_BITS;
static constexpr int      BLOCK_SIZE = 500;

// ---------------------------------------------------------------------------
// FSE / tANS
// ---------------------------------------------------------------------------

struct FseDecodeEntry {
    uint8_t  sym;
    uint8_t  nb_bits;
    uint16_t base;
};

struct FseTable {
    uint32_t freq[256] = {};
    FseDecodeEntry dec[SCALE];
    uint16_t enc_flat[SCALE];
    uint32_t enc_offset[256];
    int      nb_base[256];

    void build() {
        uint32_t pos = 0;
        const uint32_t step = (SCALE >> 1) + (SCALE >> 3) + 3;
        uint8_t spread[SCALE];
        for (int s = 0; s < 256; s++)
            for (uint32_t n = 0; n < freq[s]; n++) {
                spread[pos] = (uint8_t)s;
                pos = (pos + step) & (SCALE - 1);
            }
        uint32_t cumul = 0;
        for (int s = 0; s < 256; s++) {
            enc_offset[s] = cumul;
            cumul += freq[s];
            nb_base[s] = (freq[s] > 0) ? std::max(0, (int)(__builtin_clz(freq[s]) + SCALE_BITS - 32)) : 0;
        }
        uint32_t next[256];
        for (int s = 0; s < 256; s++) next[s] = freq[s];
        for (uint32_t state = 0; state < SCALE; state++) {
            uint8_t sym = spread[state];
            uint32_t x  = next[sym]++;
            uint8_t nb  = (uint8_t)(SCALE_BITS - (31u - __builtin_clz(x)));
            uint32_t base = (x << nb) - SCALE;
            dec[state]  = {sym, nb, (uint16_t)base};
            enc_flat[enc_offset[sym] + (x - freq[sym])] = (uint16_t)state;
        }
    }
};

static void fit_freqs(const uint64_t* cnt, uint32_t* freq) {
    uint64_t total = 0;
    for (int i = 0; i < 256; i++) total += cnt[i];
    uint32_t used = 0;
    for (int i = 0; i < 256; i++) {
        if (!cnt[i]) { freq[i] = 0; continue; }
        freq[i] = (uint32_t)std::max((uint64_t)1, cnt[i] * (uint64_t)SCALE / total);
        used += freq[i];
    }
    int peak = 0;
    for (int i = 1; i < 256; i++) if (cnt[i] > cnt[peak]) peak = i;
    if (used < SCALE) freq[peak] += SCALE - used;
    else              freq[peak] -= used - SCALE;
}

struct BitWriter {
    struct Op { uint16_t bits; uint8_t nbits; };
    std::vector<Op> ops;
    std::vector<uint8_t> bytes;

    void put(uint32_t bits, int n) { ops.push_back({(uint16_t)bits, (uint8_t)n}); }

    void flush_reverse() {
        uint64_t bitbuf = 0; int bitcnt = 0;
        for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
            bitbuf |= (uint64_t)it->bits << bitcnt;
            bitcnt += it->nbits;
            while (bitcnt >= 8) { bytes.push_back((uint8_t)(bitbuf & 0xFF)); bitbuf >>= 8; bitcnt -= 8; }
        }
        if (bitcnt > 0) bytes.push_back((uint8_t)(bitbuf & 0xFF));
    }
};

struct FseEncoder {
    uint16_t state = 0;
    BitWriter bw;

    void push(uint8_t sym, const FseTable& t) {
        uint32_t xs = (uint32_t)state + SCALE;
        int nb = t.nb_base[sym];
        if ((xs >> nb) >= 2 * t.freq[sym]) ++nb;
        bw.put(xs & ((1u << nb) - 1u), nb);
        state = t.enc_flat[t.enc_offset[sym] + (xs >> nb) - t.freq[sym]];
    }

    void finish(std::vector<uint8_t>& out) {
        bw.flush_reverse();
        out.push_back((uint8_t)(state & 0xFF));
        out.push_back((uint8_t)(state >> 8));
        out.insert(out.end(), bw.bytes.begin(), bw.bytes.end());
    }
};

// ---------------------------------------------------------------------------
// Buffer / file helpers
// ---------------------------------------------------------------------------

static void bput8  (std::vector<uint8_t>& b, uint8_t  v) { b.push_back(v); }
static void bput32 (std::vector<uint8_t>& b, uint32_t v) {
    for (int i = 0; i < 4; i++) { b.push_back(v & 0xFF); v >>= 8; }
}
static void bput64 (std::vector<uint8_t>& b, uint64_t v) {
    for (int i = 0; i < 8; i++) { b.push_back(v & 0xFF); v >>= 8; }
}
static void write_u8   (FILE* f, uint8_t  v) { fwrite(&v, 1, 1, f); }
static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4]; for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 4, f);
}

// ---------------------------------------------------------------------------
// Encode a byte stream: [flag][table][size8][state+bits]
// ---------------------------------------------------------------------------

static void encode_stream(std::vector<uint8_t>& buf, const std::vector<uint8_t>& bytes) {
    uint64_t cnt[256] = {};
    for (uint8_t b : bytes) cnt[b]++;
    int nnz = 0;
    for (int i = 0; i < 256; i++) if (cnt[i]) nnz++;

    FseTable tab;
    if (nnz == 256) {
        for (int i = 0; i < 256; i++) tab.freq[i] = SCALE / 256;
        tab.build();
        bput8(buf, 0xFF);
    } else {
        fit_freqs(cnt, tab.freq);
        tab.build();
        if (nnz <= 254) {
            bput8(buf, (uint8_t)nnz);
            for (int i = 0; i < 256; i++) { if (!tab.freq[i]) continue; bput8(buf, (uint8_t)i); bput32(buf, tab.freq[i]); }
        } else {
            bput8(buf, 0);
            for (int i = 0; i < 256; i++) bput32(buf, tab.freq[i]);
        }
    }

    FseEncoder enc;
    enc.bw.ops.reserve(bytes.size());
    std::vector<uint8_t> stream;
    for (int i = (int)bytes.size() - 1; i >= 0; i--) enc.push(bytes[i], tab);
    enc.finish(stream);
    bput64(buf, (uint64_t)stream.size());
    buf.insert(buf.end(), stream.begin(), stream.end());
}

// ---------------------------------------------------------------------------
// Predictors
// ---------------------------------------------------------------------------

static inline uint16_t zigzag(int16_t v) {
    return (v >= 0) ? (uint16_t)(v * 2) : (uint16_t)((-v) * 2 - 1);
}

static inline uint16_t avg_pred(const std::vector<uint16_t>& blk, int x, int y, int bw) {
    if (y == 0 && x == 0) return 0;
    if (y == 0)            return blk[x - 1];
    if (x == 0)            return blk[(y - 1) * bw + x];
    int A = blk[y * bw + (x - 1)];
    int B = blk[(y - 1) * bw + x];
    int C = blk[(y - 1) * bw + (x - 1)];
    return (uint16_t)((A + B + C + 1) / 3);
}

// ---------------------------------------------------------------------------
// Generic N×N Gauss-Jordan solver
// ---------------------------------------------------------------------------

template<int N>
static void ls_solve_n(double XtX[N][N], double Xty[N], float w[N]) {
    double M[N][N+1];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) M[i][j] = XtX[i][j];
        M[i][N] = Xty[i];
    }
    for (int col = 0; col < N; col++) {
        int piv = col;
        for (int r = col+1; r < N; r++)
            if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
        if (piv != col)
            for (int k = 0; k <= N; k++) std::swap(M[col][k], M[piv][k]);
        double d = M[col][col];
        if (std::abs(d) < 1e-8) { w[col] = 0.0f; continue; }
        for (int r = 0; r < N; r++) {
            if (r == col) continue;
            double f = M[r][col] / d;
            for (int k = col; k <= N; k++) M[r][k] -= f * M[col][k];
        }
    }
    for (int i = 0; i < N; i++)
        w[i] = (std::abs(M[i][i]) > 1e-8) ? (float)(M[i][N] / M[i][i]) : 0.0f;
}

// ---------------------------------------------------------------------------
// Mode 2 — LS(W, N, NW): 3 float32 weights stored per block
// ---------------------------------------------------------------------------

static void ls_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                       float w[3], std::vector<uint16_t>& syms) {
    double XtX[3][3] = {}, Xty[3] = {};
    for (int y = 1; y < bh; y++)
        for (int x = 1; x < bw; x++) {
            double f[3] = { (double)blk[y*bw+(x-1)], (double)blk[(y-1)*bw+x], (double)blk[(y-1)*bw+(x-1)] };
            double t = blk[y*bw+x];
            for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) XtX[i][j] += f[i]*f[j]; Xty[i] += f[i]*t; }
        }
    ls_solve_n<3>(XtX, Xty, w);
    int npix = bw * bh;
    syms.resize(npix);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int idx = y*bw+x;
            uint16_t pred;
            if (y == 0 && x == 0) pred = 0;
            else if (y == 0)      pred = blk[x-1];
            else if (x == 0)      pred = blk[(y-1)*bw+x];
            else {
                float p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x] + w[2]*blk[(y-1)*bw+(x-1)];
                pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
            }
            syms[idx] = zigzag((int16_t)(blk[idx]-pred));
        }
}

// ---------------------------------------------------------------------------
// Mode 4 — LS+bias(W, N, NW, 1): 4 float32 weights stored per block
// ---------------------------------------------------------------------------

static void ls4_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                        float w[4], std::vector<uint16_t>& syms) {
    double XtX[4][4] = {}, Xty[4] = {};
    for (int y = 1; y < bh; y++)
        for (int x = 1; x < bw; x++) {
            double f[4] = { (double)blk[y*bw+(x-1)], (double)blk[(y-1)*bw+x], (double)blk[(y-1)*bw+(x-1)], 1.0 };
            double t = blk[y*bw+x];
            for (int i = 0; i < 4; i++) { for (int j = 0; j < 4; j++) XtX[i][j] += f[i]*f[j]; Xty[i] += f[i]*t; }
        }
    ls_solve_n<4>(XtX, Xty, w);
    int npix = bw * bh;
    syms.resize(npix);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int idx = y*bw+x;
            uint16_t pred;
            if (y == 0 && x == 0) pred = (uint16_t)std::max(0.0f,std::min(65535.0f,w[3]+0.5f));
            else if (y == 0)      pred = (uint16_t)std::max(0.0f,std::min(65535.0f,w[0]*blk[x-1]+w[3]+0.5f));
            else if (x == 0)      pred = (uint16_t)std::max(0.0f,std::min(65535.0f,w[1]*blk[(y-1)*bw+x]+w[3]+0.5f));
            else {
                float p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x] + w[2]*blk[(y-1)*bw+(x-1)] + w[3];
                pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
            }
            syms[idx] = zigzag((int16_t)(blk[idx]-pred));
        }
}

// ---------------------------------------------------------------------------
// Symbol generation per mode (modes 0-1)
// ---------------------------------------------------------------------------

static void make_syms(const std::vector<uint16_t>& blk, int mode, int bw, int bh,
                      std::vector<uint16_t>& syms, uint16_t gmean) {
    int npix = bw * bh;
    syms.resize(npix);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int idx = y * bw + x;
            uint16_t pix = blk[idx];
            if (mode == 0) { syms[idx] = pix; continue; }
            if (mode == 3) { syms[idx] = zigzag((int16_t)(pix - (int)gmean)); continue; }
            syms[idx] = zigzag((int16_t)(pix - avg_pred(blk, x, y, bw)));
        }
}

// ---------------------------------------------------------------------------
// Byte-level cost estimate: H(hi8) + H(lo8)
// ---------------------------------------------------------------------------

static double byte_cost(const std::vector<uint16_t>& syms) {
    uint64_t hfreq[256] = {}, lfreq[256] = {};
    for (uint16_t s : syms) { hfreq[s >> 8]++; lfreq[s & 0xFF]++; }
    double H = 0.0, N = (double)syms.size();
    for (int i = 0; i < 256; i++) {
        if (hfreq[i]) { double p = hfreq[i] / N; H -= p * std::log2(p); }
        if (lfreq[i]) { double p = lfreq[i] / N; H -= p * std::log2(p); }
    }
    return H;
}

// ---------------------------------------------------------------------------
// Per-block compression result
// ---------------------------------------------------------------------------

struct BlockResult {
    uint8_t              mode;
    std::vector<uint8_t> payload;
};

static BlockResult compress_block(const std::vector<uint16_t>& image,
                                   int W, int H, int bs, int bx, int by, uint16_t gmean) {
    int bw   = std::min(bs, W - bx * bs);
    int bh   = std::min(bs, H - by * bs);
    int npix = bw * bh;

    std::vector<uint16_t> blk(npix);
    for (int y = 0; y < bh; y++) {
        int gy = by * bs + y;
        for (int x = 0; x < bw; x++)
            blk[y * bw + x] = image[gy * W + bx * bs + x];
    }

    std::vector<uint16_t> s0, s1, s2, s3, s4;
    float ls_w[3], ls4_w[4];

    make_syms(blk, 0, bw, bh, s0, gmean);   // raw
    make_syms(blk, 1, bw, bh, s1, gmean);   // avg
    ls_compute (blk, bw, bh, ls_w,  s2);    // LS(W,N,NW)
    make_syms(blk, 3, bw, bh, s3, gmean);   // global_mean
    ls4_compute(blk, bw, bh, ls4_w, s4);    // LS+bias

    const double ls_overhead  = (12.0 * 8.0) / npix;
    const double ls4_overhead = (16.0 * 8.0) / npix;
    const std::vector<uint16_t>* sp[5] = {&s0, &s1, &s2, &s3, &s4};
    double costs[5] = { byte_cost(s0), byte_cost(s1),
                        byte_cost(s2) + ls_overhead,
                        byte_cost(s3),
                        byte_cost(s4) + ls4_overhead };
    int best = 0;
    for (int m = 1; m < 5; m++)
        if (costs[m] < costs[best]) best = m;

    const auto& sbest = *sp[best];
    std::vector<uint8_t> hi(npix), lo(npix);
    for (int i = 0; i < npix; i++) { hi[i] = sbest[i] >> 8; lo[i] = sbest[i] & 0xFF; }

    BlockResult r;
    r.mode = (uint8_t)best;
    if (best == 2) {
        uint32_t bits;
        for (int i = 0; i < 3; i++) { memcpy(&bits, &ls_w[i], 4); bput32(r.payload, bits); }
    } else if (best == 4) {
        uint32_t bits;
        for (int i = 0; i < 4; i++) { memcpy(&bits, &ls4_w[i], 4); bput32(r.payload, bits); }
    }
    encode_stream(r.payload, hi);
    encode_stream(r.payload, lo);
    return r;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 3 || argc == 4 || argc > 5) {
        fprintf(stderr, "Usage: %s <input> <output.hais> [width height]\n", argv[0]);
        return 1;
    }

    int W = 1500, H = 1500;
    if (argc == 5) {
        W = std::atoi(argv[3]); H = std::atoi(argv[4]);
        if (W <= 0 || H <= 0) { fprintf(stderr, "Invalid dimensions\n"); return 1; }
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    fseek(fin, 0, SEEK_END); long fsize = ftell(fin); rewind(fin);
    if (fsize != (long)W * H * 2) {
        fprintf(stderr, "File size mismatch: got %ld, expected %ld\n", fsize, (long)W*H*2);
        fclose(fin); return 1;
    }
    std::vector<uint16_t> image(W * H);
    for (int i = 0; i < W * H; i++) {
        uint8_t b[2];
        if (fread(b, 1, 2, fin) != 2) { fprintf(stderr, "Short read\n"); return 1; }
        image[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    uint64_t sum = 0;
    for (auto px : image) sum += px;
    uint16_t gmean = (uint16_t)(sum / (W * H));

    int bs       = BLOCK_SIZE;
    int blocks_x = (W + bs - 1) / bs;
    int blocks_y = (H + bs - 1) / bs;
    int total    = blocks_x * blocks_y;

    std::vector<BlockResult> results(total);
    std::atomic<int> next{0};
    int nthreads = std::max(1, (int)std::thread::hardware_concurrency());

    auto worker = [&]() {
        int idx;
        while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < total)
            results[idx] = compress_block(image, W, H, bs, idx % blocks_x, idx / blocks_x, gmean);
    };
    std::vector<std::thread> pool(nthreads);
    for (auto& t : pool) t = std::thread(worker);
    for (auto& t : pool) t.join();

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u32le(fout, (uint32_t)W);
    write_u32le(fout, (uint32_t)H);
    write_u32le(fout, (uint32_t)bs);
    write_u16le(fout, gmean);

    int mode_count[5] = {};
    for (int idx = 0; idx < total; idx++) {
        const auto& r = results[idx];
        mode_count[r.mode]++;
        write_u8(fout, r.mode);
        fwrite(r.payload.data(), 1, r.payload.size(), fout);
    }
    fclose(fout);

    long in_bytes = (long)W * H * 2;
    FILE* ft = fopen(argv[2], "rb"); fseek(ft, 0, SEEK_END); long out_bytes = ftell(ft); fclose(ft);
    fprintf(stderr, "Dimensões: %dx%d  blocos: %dx%d (%d total)\n", W, H, blocks_x, blocks_y, total);
    fprintf(stderr, "Modos: raw=%d  avg=%d  ls=%d  mean=%d  ls4=%d\n",
            mode_count[0], mode_count[1], mode_count[2], mode_count[3], mode_count[4]);
    fprintf(stderr, "Comprimido: %ld -> %ld bytes  (%.4f bits/byte)\n",
            in_bytes, out_bytes, out_bytes * 8.0 / in_bytes);
    return 0;
}
