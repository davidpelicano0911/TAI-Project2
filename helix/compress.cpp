//
// HELIX v7 — per-block predictor + 3-level JPEG-LS context + per-context FSE
//
// v7 improvements over v6:
//   - 3-level gradient quantization (T1 = p50 of |gradient|, adaptive for 16-bit)
//     → 27 contexts (14 valid after sign canonicalization)
//     → All contexts well-populated for 16-bit astronomical images (was 86% in 1 bin)
//   - MAE-based predictor selection (reverted from ctx_cost which was worse)
//
// File format (little-endian):
//   [4]  magic "HELX"
//   [2]  version = 7
//   [4]  width
//   [4]  height
//   [2]  flags
//   [2]  gmean
//   [4]  T1  (gradient threshold for 3-level quantization, uint32)
//   [4]  block_size
//   [4]  nblocks_x
//   [4]  nblocks_y
//   [nblocks]     block mode bytes (0=MED,1=avg,2=mean,3=LS6)
//   [nblocks*24]  float LS6 weights (6 per block, zeros for non-LS6)
//   [4]  N_CTX = 27
//   for each context q (0..26):
//     [4] n (bit31 = lo_split)
//     if n>0: [FSE hi][FSE top] then [FSE lo] or [FSE lo_zero][FSE lo_nonzero]
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

// ---------------------------------------------------------------------------
// FSE / tANS  (inlined)
// ---------------------------------------------------------------------------

static constexpr uint32_t FSE_SCALE_BITS = 16;
static constexpr uint32_t FSE_SCALE      = 1u << FSE_SCALE_BITS;

struct FseDecodeEntry { uint8_t sym; uint8_t nb_bits; uint16_t base; };

struct FseTable {
    uint32_t       freq[256]  = {};
    FseDecodeEntry dec[FSE_SCALE];
    uint16_t       enc_flat[FSE_SCALE];
    uint32_t       enc_offset[256];
    int            nb_base[256];

    void build() {
        const uint32_t step = (FSE_SCALE >> 1) + (FSE_SCALE >> 3) + 3;
        uint8_t spread[FSE_SCALE];
        uint32_t pos = 0;
        for (int s = 0; s < 256; s++)
            for (uint32_t n = 0; n < freq[s]; n++) {
                spread[pos] = (uint8_t)s;
                pos = (pos + step) & (FSE_SCALE - 1);
            }
        uint32_t cumul = 0;
        for (int s = 0; s < 256; s++) {
            enc_offset[s] = cumul;
            cumul += freq[s];
            nb_base[s] = (freq[s] > 0)
                ? std::max(0, (int)(__builtin_clz(freq[s]) + FSE_SCALE_BITS - 32))
                : 0;
        }
        uint32_t next[256];
        for (int s = 0; s < 256; s++) next[s] = freq[s];
        for (uint32_t state = 0; state < FSE_SCALE; state++) {
            uint8_t  sym  = spread[state];
            uint32_t x    = next[sym]++;
            uint8_t  nb   = (uint8_t)(FSE_SCALE_BITS - (31u - __builtin_clz(x)));
            uint32_t base = (x << nb) - FSE_SCALE;
            dec[state]    = {sym, nb, (uint16_t)base};
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
        freq[i] = (uint32_t)std::max((uint64_t)1, cnt[i] * (uint64_t)FSE_SCALE / total);
        used += freq[i];
    }
    int peak = 0;
    for (int i = 1; i < 256; i++) if (cnt[i] > cnt[peak]) peak = i;
    if (used < FSE_SCALE) freq[peak] += FSE_SCALE - used;
    else                  freq[peak] -= used - FSE_SCALE;
}

static inline void bput8 (std::vector<uint8_t>& b, uint8_t  v) { b.push_back(v); }
static inline void bput32(std::vector<uint8_t>& b, uint32_t v) {
    for (int i = 0; i < 4; i++) { b.push_back(v & 0xFF); v >>= 8; }
}
static inline void bput64(std::vector<uint8_t>& b, uint64_t v) {
    for (int i = 0; i < 8; i++) { b.push_back(v & 0xFF); v >>= 8; }
}

struct FseBitWriter {
    struct Op { uint16_t bits; uint8_t nbits; };
    std::vector<Op> ops;
    std::vector<uint8_t> bytes;
    void put(uint32_t bits, int n) { ops.push_back({(uint16_t)bits, (uint8_t)n}); }
    void flush_reverse() {
        uint64_t bitbuf = 0; int bitcnt = 0;
        for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
            bitbuf |= (uint64_t)it->bits << bitcnt;
            bitcnt += it->nbits;
            while (bitcnt >= 8) {
                bytes.push_back((uint8_t)(bitbuf & 0xFF));
                bitbuf >>= 8; bitcnt -= 8;
            }
        }
        if (bitcnt > 0) bytes.push_back((uint8_t)(bitbuf & 0xFF));
    }
};

struct FseEncoder {
    uint16_t state = 0;
    FseBitWriter bw;
    void push(uint8_t sym, const FseTable& t) {
        uint32_t xs = (uint32_t)state + FSE_SCALE;
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

static void encode_stream(std::vector<uint8_t>& buf, const std::vector<uint8_t>& bytes) {
    if (bytes.empty()) return;
    uint64_t cnt[256] = {};
    for (uint8_t b : bytes) cnt[b]++;
    int nnz = 0;
    for (int i = 0; i < 256; i++) if (cnt[i]) nnz++;

    if (nnz == 1) {
        int sym = 0; while (!cnt[sym]) sym++;
        bput8(buf, 0x01); bput8(buf, (uint8_t)sym);
        return;
    }
    FseTable tab;
    fit_freqs(cnt, tab.freq);
    tab.build();
    if (nnz <= 254) {
        bput8(buf, (uint8_t)nnz);
        for (int i = 0; i < 256; i++) {
            if (!tab.freq[i]) continue;
            bput8(buf, (uint8_t)i); bput32(buf, tab.freq[i]);
        }
    } else {
        bput8(buf, 0);
        for (int i = 0; i < 256; i++) bput32(buf, tab.freq[i]);
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
// Constants
// ---------------------------------------------------------------------------

static constexpr uint8_t  MAGIC[4]   = {'H','E','L','X'};
static constexpr uint16_t VERSION    = 7;
static constexpr int      N_CTX      = 27;   // 3^3 (14 valid after sign canonicalization)
static constexpr int      BLOCK_SIZE = 500;

// ---------------------------------------------------------------------------
// LS6 solver: Gauss-Jordan for 6×6
// ---------------------------------------------------------------------------

static void ls6_solve(double XtX[6][6], double Xty[6], float w[6]) {
    double M[6][7];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) M[i][j] = XtX[i][j];
        M[i][6] = Xty[i];
    }
    for (int col = 0; col < 6; col++) {
        int piv = col;
        for (int r = col+1; r < 6; r++)
            if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
        if (piv != col)
            for (int k = 0; k <= 6; k++) std::swap(M[col][k], M[piv][k]);
        double d = M[col][col];
        if (std::abs(d) < 1e-8) { w[col] = 0.0f; continue; }
        for (int r = 0; r < 6; r++) {
            if (r == col) continue;
            double f = M[r][col] / d;
            for (int k = col; k <= 6; k++) M[r][k] -= f * M[col][k];
        }
    }
    for (int i = 0; i < 6; i++)
        w[i] = (std::abs(M[i][i]) > 1e-8) ? (float)(M[i][6] / M[i][i]) : 0.0f;
}

// ---------------------------------------------------------------------------
// Predictors
// ---------------------------------------------------------------------------

static inline int32_t med_predict(int32_t a, int32_t b, int32_t c) {
    int32_t lo = std::min(a, b);
    int32_t hi = std::max(a, b);
    if      (c >= hi) return lo;
    else if (c <= lo) return hi;
    else              return a + b - c;
}

static inline int32_t avg_predict(int32_t a, int32_t b, int32_t c,
                                   int x, int y) {
    if (y == 0 && x == 0) return 0;
    if (y == 0)            return a;
    if (x == 0)            return b;
    return (a + b + c + 1) / 3;
}

// ---------------------------------------------------------------------------
// 3-level gradient quantization (adaptive T1 = p50 of |gradient|)
// ---------------------------------------------------------------------------

static inline int quantize_grad3(int32_t g, int32_t T1) {
    if (g < -T1) return -1;
    if (g <=  T1) return  0;
    return  1;
}

static inline int ctx_index3(int q1, int q2, int q3) {
    return (q1 + 1) * 9 + (q2 + 1) * 3 + (q3 + 1);
}

// Compute T1 = median of |gradient| values across the whole image
static int32_t compute_threshold(const std::vector<uint16_t>& raw,
                                  int width, int height) {
    std::vector<int32_t> mags;
    mags.reserve((size_t)(width - 2) * (height - 1) * 3);
    for (int y = 1; y < height; y++) {
        for (int x = 1; x < width - 1; x++) {
            int32_t d  = raw[(y-1)*width + (x+1)];
            int32_t b2 = raw[(y-1)*width + x];
            int32_t c2 = raw[(y-1)*width + (x-1)];
            int32_t a  = raw[y*width + (x-1)];
            mags.push_back(std::abs(d  - b2));
            mags.push_back(std::abs(b2 - c2));
            mags.push_back(std::abs(c2 - a));
        }
    }
    std::sort(mags.begin(), mags.end());
    return std::max((int32_t)1, mags[mags.size() / 2]);
}

static inline uint32_t map_signed(int32_t v) {
    return (v >= 0) ? (uint32_t)(2 * v) : (uint32_t)(-2 * v - 1);
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 4, f);
}

// ---------------------------------------------------------------------------
// Compress
// ---------------------------------------------------------------------------

static int compress(const char* in_path, const char* out_path) {
    FILE* fin = fopen(in_path, "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", in_path); return 1; }
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    rewind(fin);

    const int width  = 1500;
    const int height = 1500;
    const int npix   = width * height;

    if (fsize != (long)npix * 2) {
        fprintf(stderr, "Unexpected file size %ld (expected %d)\n", fsize, npix * 2);
        fclose(fin); return 1;
    }

    std::vector<uint16_t> raw(npix);
    for (int i = 0; i < npix; i++) {
        uint8_t b[2]; fread(b, 1, 2, fin);
        raw[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    // Compute global mean and adaptive gradient threshold
    uint64_t sum = 0;
    for (auto px : raw) sum += px;
    uint16_t gmean = (uint16_t)(sum / npix);

    int32_t T1 = compute_threshold(raw, width, height);
    fprintf(stderr, "T1 = %d (gradient threshold)\n", (int)T1);

    // Work on int32 values
    std::vector<int32_t> img(npix);
    for (int i = 0; i < npix; i++) img[i] = (int32_t)raw[i];
    uint16_t flags = 0x0004;

    // Per-block: fit LS6, select best mode from {MED(0), avg(1), mean(2), LS6(3)}
    const int BS   = BLOCK_SIZE;
    const int NBX  = (width  + BS - 1) / BS;  // 3
    const int NBY  = (height + BS - 1) / BS;  // 3
    const int NBLK = NBX * NBY;               // 9

    std::vector<std::array<float, 6>> bw(NBLK);
    std::vector<uint8_t> blk_mode(NBLK, 0);

    for (int by = 0; by < NBY; by++) {
        for (int bx = 0; bx < NBX; bx++) {
            int bidx = by * NBX + bx;
            int y0 = by * BS, y1 = std::min(y0 + BS, height);
            int x0 = bx * BS, x1 = std::min(x0 + BS, width);

            // Fit LS6: uses W, N, NW, WW, NN, bias
            double XtX[6][6] = {}, Xty[6] = {};
            long count = 0;
            for (int y = y0+2; y < y1; y++) {
                for (int x = x0+2; x < x1; x++) {
                    double f[6] = {
                        (double)img[y*width+(x-1)],
                        (double)img[(y-1)*width+x],
                        (double)img[(y-1)*width+(x-1)],
                        (double)img[y*width+(x-2)],
                        (double)img[(y-2)*width+x],
                        1.0
                    };
                    double t = img[y*width+x];
                    for (int i = 0; i < 6; i++) {
                        for (int j = 0; j < 6; j++) XtX[i][j] += f[i]*f[j];
                        Xty[i] += f[i]*t;
                    }
                    ++count;
                }
            }
            if (count >= 6) {
                ls6_solve(XtX, Xty, bw[bidx].data());
            } else {
                bw[bidx] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            }

            // Evaluate modes using MAE (correlates better with entropy for helix's context model)
            double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
            int nr01 = 0, nr3 = 0;

            for (int y = y0+1; y < y1; y++) {
                for (int x = x0+1; x < x1; x++) {
                    int32_t pix = img[y*width+x];
                    int32_t a   = img[y*width+(x-1)];
                    int32_t b2  = img[(y-1)*width+x];
                    int32_t c   = img[(y-1)*width+(x-1)];

                    sum0 += std::abs(pix - med_predict(a, b2, c));
                    sum1 += std::abs(pix - avg_predict(a, b2, c, x, y));
                    sum2 += std::abs(pix - (int32_t)gmean);
                    nr01++;

                    if (y >= y0+2 && x >= x0+2) {
                        const auto& w = bw[bidx];
                        float p = w[0]*img[y*width+(x-1)] + w[1]*img[(y-1)*width+x]
                                + w[2]*img[(y-1)*width+(x-1)] + w[3]*img[y*width+(x-2)]
                                + w[4]*img[(y-2)*width+x] + w[5];
                        sum3 += std::abs(pix - std::max(0, std::min(65535, (int32_t)(p+0.5f))));
                        nr3++;
                    }
                }
            }

            double cost0 = nr01 > 0 ? sum0 / nr01 : 1e18;
            double cost1 = nr01 > 0 ? sum1 / nr01 : 1e18;
            double cost2 = nr01 > 0 ? sum2 / nr01 : 1e18;
            double cost3 = (nr3 > 0 && count >= 6) ? sum3 / nr3 : 1e18;

            double costs[4] = {cost0, cost1, cost2, cost3};
            int best = 0;
            for (int m = 1; m < 4; m++)
                if (costs[m] < costs[best]) best = m;
            blk_mode[bidx] = (uint8_t)best;
        }
    }

    // Log mode selection
    int mode_cnt[4] = {};
    for (int b = 0; b < NBLK; b++) mode_cnt[blk_mode[b]]++;
    fprintf(stderr, "Modes: MED=%d avg=%d mean=%d LS6=%d\n",
            mode_cnt[0], mode_cnt[1], mode_cnt[2], mode_cnt[3]);

    // Main pixel loop: predict → JPEG-LS 3-level context → collect FSE symbols
    std::vector<std::vector<uint8_t>> ctx_lo(N_CTX), ctx_hi(N_CTX), ctx_top(N_CTX);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int32_t a = (x > 0)                  ? img[y * width + (x - 1)]        : 0;
            int32_t b = (y > 0)                  ? img[(y - 1) * width + x]        : 0;
            int32_t c = (x > 0 && y > 0)         ? img[(y - 1) * width + (x - 1)] : 0;
            int32_t d = (y > 0 && x + 1 < width) ? img[(y - 1) * width + (x + 1)] : b;

            // 3-level JPEG-LS context with adaptive T1
            int32_t g1 = d - b, g2 = b - c, g3 = c - a;
            int q1 = quantize_grad3(g1, T1);
            int q2 = quantize_grad3(g2, T1);
            int q3 = quantize_grad3(g3, T1);
            int sign = 1;
            if (q1 < 0 || (q1 == 0 && q2 < 0) || (q1 == 0 && q2 == 0 && q3 < 0)) {
                sign = -1;
                q1 = -q1; q2 = -q2; q3 = -q3;
            }
            int Q = ctx_index3(q1, q2, q3);

            // Block-selected prediction
            int bidx_cur = (y / BS) * NBX + (x / BS);
            int mode = blk_mode[bidx_cur];
            int32_t pred;
            if (mode == 0) {
                if (y == 0 && x == 0) pred = 0;
                else if (y == 0)      pred = a;
                else if (x == 0)      pred = b;
                else                  pred = med_predict(a, b, c);
            } else if (mode == 1) {
                pred = avg_predict(a, b, c, x, y);
            } else if (mode == 2) {
                pred = (int32_t)gmean;
            } else {
                // LS6 with border fallback to MED
                if (y < 2 || x < 2) {
                    if (y == 0 && x == 0) pred = 0;
                    else if (y == 0)      pred = a;
                    else if (x == 0)      pred = b;
                    else                  pred = med_predict(a, b, c);
                } else {
                    const auto& w = bw[bidx_cur];
                    float p = w[0]*img[y*width+(x-1)]
                            + w[1]*img[(y-1)*width+x]
                            + w[2]*img[(y-1)*width+(x-1)]
                            + w[3]*img[y*width+(x-2)]
                            + w[4]*img[(y-2)*width+x]
                            + w[5];
                    pred = std::max(0, std::min(65535, (int32_t)(p + 0.5f)));
                }
            }

            int32_t e_raw = img[y * width + x] - pred;
            int32_t e_c   = sign * e_raw;
            uint32_t u    = map_signed(e_c);

            ctx_lo[Q].push_back((uint8_t)(u & 0xFF));
            ctx_hi[Q].push_back((uint8_t)((u >> 8) & 0xFF));
            ctx_top[Q].push_back((uint8_t)(u >> 16));
        }
    }

    // FSE encode per-context with lo-split optimization
    std::vector<uint8_t> data_buf;
    data_buf.reserve(npix);

    for (int q = 0; q < N_CTX; q++) {
        uint32_t n = (uint32_t)ctx_lo[q].size();
        if (n == 0) {
            for (int i = 0; i < 4; i++) data_buf.push_back(0);
            continue;
        }

        std::vector<uint8_t> lo_zero, lo_nonzero;
        lo_zero.reserve(n); lo_nonzero.reserve(n);
        for (uint32_t i = 0; i < n; i++) {
            if (ctx_hi[q][i] == 0 && ctx_top[q][i] == 0)
                lo_zero.push_back(ctx_lo[q][i]);
            else
                lo_nonzero.push_back(ctx_lo[q][i]);
        }

        std::vector<uint8_t> plain_hi, plain_top, plain_lo;
        encode_stream(plain_hi,  ctx_hi[q]);
        encode_stream(plain_top, ctx_top[q]);
        encode_stream(plain_lo,  ctx_lo[q]);
        size_t plain_size = plain_hi.size() + plain_top.size() + plain_lo.size();

        std::vector<uint8_t> split_lo_zero, split_lo_nonzero;
        encode_stream(split_lo_zero,    lo_zero);
        encode_stream(split_lo_nonzero, lo_nonzero);
        size_t split_size = plain_hi.size() + plain_top.size()
                          + split_lo_zero.size() + split_lo_nonzero.size();

        bool use_split = (split_size < plain_size);
        uint32_t n_flag = use_split ? (n | 0x80000000u) : n;
        for (int i = 0; i < 4; i++) { data_buf.push_back(n_flag & 0xFF); n_flag >>= 8; }

        data_buf.insert(data_buf.end(), plain_hi.begin(),  plain_hi.end());
        data_buf.insert(data_buf.end(), plain_top.begin(), plain_top.end());
        if (use_split) {
            data_buf.insert(data_buf.end(), split_lo_zero.begin(),    split_lo_zero.end());
            data_buf.insert(data_buf.end(), split_lo_nonzero.begin(), split_lo_nonzero.end());
        } else {
            data_buf.insert(data_buf.end(), plain_lo.begin(), plain_lo.end());
        }
    }

    // Write output
    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u16le(fout, VERSION);
    write_u32le(fout, (uint32_t)width);
    write_u32le(fout, (uint32_t)height);
    write_u16le(fout, flags);
    write_u16le(fout, gmean);
    write_u32le(fout, (uint32_t)T1);
    write_u32le(fout, (uint32_t)BS);
    write_u32le(fout, (uint32_t)NBX);
    write_u32le(fout, (uint32_t)NBY);
    fwrite(blk_mode.data(), 1, NBLK, fout);
    for (int b = 0; b < NBLK; b++) {
        for (int k = 0; k < 6; k++) {
            uint32_t bits;
            memcpy(&bits, &bw[b][k], 4);
            write_u32le(fout, bits);
        }
    }
    write_u32le(fout, (uint32_t)N_CTX);
    fwrite(data_buf.data(), 1, data_buf.size(), fout);

    long out_size = (long)(4+2+4+4+2+2+4 + 4+4+4 + NBLK + NBLK*6*4 + 4 + data_buf.size());
    fprintf(stderr, "Compressed: %ld -> %ld bytes (%.4f bpp)\n",
            fsize, out_size, (double)out_size * 8.0 / npix);
    fclose(fout);
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 3) { fprintf(stderr, "Usage: %s <input> <output.hlx>\n", argv[0]); return 1; }
    return compress(argv[1], argv[2]);
}
