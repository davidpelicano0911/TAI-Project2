// HAIS — Hybrid Astronomical Image compressor  (version 4)
//
// Pipeline por bloco 150×150:
//   1. Gerar símbolos com 3 modos: raw / MED / avg
//   2. Escolher modo com menor byte_cost (H(hi8)+H(lo8))
//   3. Comprimir com FSE/tANS (split hi8+lo8)
//
// Usage: ./compress <input> <output.hais> [width height]
//        Default width=1500 height=1500 (backward compatible with benchmark)

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
static constexpr int      BLOCK_SIZE = 150;

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
    uint16_t enc_flat[SCALE];   // new state for each (sym, slot)
    uint32_t enc_offset[256];   // cumulative freq offsets into enc_flat
    int      nb_base[256];      // __builtin_clz(freq[s]) + SCALE_BITS - 32

    void build() {
        uint32_t pos = 0;
        const uint32_t step = (SCALE >> 1) + (SCALE >> 3) + 3;
        uint8_t spread[SCALE];
        for (int s = 0; s < 256; s++) {
            for (uint32_t n = 0; n < freq[s]; n++) {
                spread[pos] = (uint8_t)s;
                pos = (pos + step) & (SCALE - 1);
            }
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
    struct Op {
        uint16_t bits;
        uint8_t  nbits;
    };

    std::vector<Op> ops;
    std::vector<uint8_t> bytes;

    void put(uint32_t bits, int n) {
        ops.push_back({(uint16_t)bits, (uint8_t)n});
    }

    void flush_reverse() {
        uint64_t bitbuf = 0;
        int bitcnt = 0;
        for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
            bitbuf |= (uint64_t)it->bits << bitcnt;
            bitcnt += it->nbits;
            while (bitcnt >= 8) {
                bytes.push_back((uint8_t)(bitbuf & 0xFF));
                bitbuf >>= 8;
                bitcnt -= 8;
            }
        }
        if (bitcnt > 0) bytes.push_back((uint8_t)(bitbuf & 0xFF));
    }
};

struct FseEncoder {
    uint16_t state = 0;
    BitWriter bw;

    void push(uint8_t sym, const FseTable& t) {
        uint32_t xs = (uint32_t)state + SCALE;   // xs in [SCALE, 2*SCALE)
        int nb = t.nb_base[sym];
        if ((xs >> nb) >= 2 * t.freq[sym]) ++nb; // O(1) edge-case adjustment
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
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4]; for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 4, f);
}

// ---------------------------------------------------------------------------
// Encode a byte stream into buf: [table][stream_size][final_state+bits]
//
// Table format:
//   0x00          dense: 256 × uint32
//   0x01..0xFE    sparse: nnz × (uint8 sym + uint32 freq)
//   0xFF          implicit uniform (freq[i]=SCALE/256 for all i) - no table bytes
//                 used when nnz==256: saves 1024 bytes with no compression loss
//                 for nearly-uniform distributions (CCD readout noise)
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
        // 0xFF is reserved for implicit-uniform; 0x00 for dense.
        // Sparse flag range is 0x01..0xFE (nnz 1..254).
        // nnz == 255 must use dense to avoid the 0xFF collision.
        if (nnz <= 254) {                 // sparse
            bput8(buf, (uint8_t)nnz);
            for (int i = 0; i < 256; i++) {
                if (!tab.freq[i]) continue;
                bput8(buf, (uint8_t)i);
                bput32(buf, tab.freq[i]);
            }
        } else {                          // dense (nnz == 255)
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

// Mode 1 — avg: dampens noise in flat regions where extrapolation amplifies fluctuations.
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
// Mode 2 — LS (Least Squares per block)
// Finds optimal weights w[W, N, NW] by solving the normal equations on the
// interior pixels of the block. Weights are stored in the bitstream (12 B
// overhead) so the decoder can reproduce the exact same prediction.
// ---------------------------------------------------------------------------

static void ls_solve(double XtX[3][3], double Xty[3], float w[3]) {
    double M[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) M[i][j] = XtX[i][j];
        M[i][3] = Xty[i];
    }
    for (int col = 0; col < 3; col++) {
        int piv = col;
        for (int r = col + 1; r < 3; r++)
            if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
        if (piv != col)
            for (int k = 0; k <= 3; k++) std::swap(M[col][k], M[piv][k]);
        double d = M[col][col];
        if (std::abs(d) < 1e-8) { w[col] = 0.0f; continue; }
        for (int r = 0; r < 3; r++) {
            if (r == col) continue;
            double f = M[r][col] / d;
            for (int k = col; k <= 3; k++) M[r][k] -= f * M[col][k];
        }
    }
    for (int i = 0; i < 3; i++)
        w[i] = (std::abs(M[i][i]) > 1e-8) ? (float)(M[i][3] / M[i][i]) : 0.0f;
}

static void ls_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                       float w[3], std::vector<uint16_t>& syms) {
    double XtX[3][3] = {}, Xty[3] = {};
    for (int y = 1; y < bh; y++) {
        for (int x = 1; x < bw; x++) {
            double f[3] = { (double)blk[y*bw+(x-1)],
                            (double)blk[(y-1)*bw+x],
                            (double)blk[(y-1)*bw+(x-1)] };
            double t = blk[y*bw+x];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) XtX[i][j] += f[i]*f[j];
                Xty[i] += f[i]*t;
            }
        }
    }
    ls_solve(XtX, Xty, w);

    int npix = bw * bh;
    syms.resize(npix);
    for (int y = 0; y < bh; y++) {
        for (int x = 0; x < bw; x++) {
            int idx = y*bw+x;
            uint16_t pix = blk[idx];
            uint16_t pred;
            if (y == 0 && x == 0) pred = 0;
            else if (y == 0)      pred = blk[x-1];
            else if (x == 0)      pred = blk[(y-1)*bw+x];
            else {
                float p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x] + w[2]*blk[(y-1)*bw+(x-1)];
                pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
            }
            syms[idx] = zigzag((int16_t)(pix - pred));
        }
    }
}

// ---------------------------------------------------------------------------
// Symbol generation per mode (modes 0-1)
// ---------------------------------------------------------------------------

static void make_syms(const std::vector<uint16_t>& blk, int mode, int bw, int bh,
                      std::vector<uint16_t>& syms) {
    int npix = bw * bh;
    syms.resize(npix);
    for (int y = 0; y < bh; y++) {
        for (int x = 0; x < bw; x++) {
            int idx = y * bw + x;
            uint16_t pix = blk[idx];
            if (mode == 0) { syms[idx] = pix; continue; }
            syms[idx] = zigzag((int16_t)(pix - avg_pred(blk, x, y, bw)));
        }
    }
}

// ---------------------------------------------------------------------------
// Byte-level cost estimate: H(hi8) + H(lo8)
// Measures what the encoder actually pays (two independent byte streams).
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

// W, H = full image dimensions; bs = block_size; bx, by = block column/row index
static BlockResult compress_block(const std::vector<uint16_t>& image,
                                   int W, int H, int bs, int bx, int by) {
    int bw   = std::min(bs, W - bx * bs);
    int bh   = std::min(bs, H - by * bs);
    int npix = bw * bh;

    std::vector<uint16_t> blk(npix);
    for (int y = 0; y < bh; y++) {
        int gy = by * bs + y;
        for (int x = 0; x < bw; x++)
            blk[y * bw + x] = image[gy * W + bx * bs + x];
    }

    std::vector<uint16_t> s0, s1, s2;
    float ls_w[3];
    make_syms(blk, 0, bw, bh, s0);
    make_syms(blk, 1, bw, bh, s1);
    ls_compute(blk, bw, bh, ls_w, s2);

    double ls_overhead = (12.0 * 8.0) / npix;
    const std::vector<uint16_t>* sp[3] = {&s0, &s1, &s2};
    double costs[3] = { byte_cost(s0), byte_cost(s1),
                        byte_cost(s2) + ls_overhead };
    int best = 0;
    for (int m = 1; m < 3; m++)
        if (costs[m] < costs[best]) best = m;

    const auto& sbest = *sp[best];
    std::vector<uint8_t> hi(npix), lo(npix);
    for (int i = 0; i < npix; i++) {
        hi[i] = (uint8_t)(sbest[i] >> 8);
        lo[i] = (uint8_t)(sbest[i] & 0xFF);
    }

    BlockResult r;
    r.mode = (uint8_t)best;
    if (best == 2) {
        uint32_t bits;
        for (int i = 0; i < 3; i++) {
            memcpy(&bits, &ls_w[i], 4);
            bput32(r.payload, bits);
        }
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
        W = std::atoi(argv[3]);
        H = std::atoi(argv[4]);
        if (W <= 0 || H <= 0) {
            fprintf(stderr, "Invalid dimensions %d x %d\n", W, H);
            return 1;
        }
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    rewind(fin);
    if (fsize != (long)W * H * 2) {
        fprintf(stderr, "Tamanho do ficheiro (%ld bytes) não corresponde a %dx%d pixels (%ld bytes esperados)\n",
                fsize, W, H, (long)W * H * 2);
        fclose(fin); return 1;
    }
    std::vector<uint16_t> image(W * H);
    for (int i = 0; i < W * H; i++) {
        uint8_t b[2];
        if (fread(b, 1, 2, fin) != 2) { fprintf(stderr, "Short read at pixel %d\n", i); return 1; }
        image[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

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
            results[idx] = compress_block(image, W, H, bs, idx % blocks_x, idx / blocks_x);
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

    int mode_count[3] = {};
    for (int idx = 0; idx < total; idx++) {
        const auto& r = results[idx];
        mode_count[r.mode]++;
        write_u8(fout, r.mode);
        fwrite(r.payload.data(), 1, r.payload.size(), fout);
    }
    fclose(fout);

    long in_bytes = (long)W * H * 2;
    FILE* ft = fopen(argv[2], "rb");
    fseek(ft, 0, SEEK_END);
    long out_bytes = ftell(ft);
    fclose(ft);

    fprintf(stderr, "Dimensões: %dx%d  blocos: %dx%d (%d total)\n",
            W, H, blocks_x, blocks_y, total);
    fprintf(stderr, "Modos: raw=%d  avg=%d  ls=%d\n",
            mode_count[0], mode_count[1], mode_count[2]);
    fprintf(stderr, "Comprimido: %ld -> %ld bytes  (%.4f bits/byte)\n",
            in_bytes, out_bytes, out_bytes * 8.0 / in_bytes);
    return 0;
}
