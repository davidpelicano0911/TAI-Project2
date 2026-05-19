//
// HELIX — context-adaptive lossless compressor for raw astronomical images
//
// Pipeline:
//   1. Per-row debiasing      (subtract row median, store medians in header)
//   2. MED predictor          (LOCO-I / JPEG-LS)
//   3. JPEG-LS context model  (3 gradients × 5 levels with sign symmetry → 63 ctx)
//   4. Adaptive bias correction per context (LOCO-I §4.3)
//   5. Adaptive Golomb-Rice per context (k derived from running |residual| mean)
//
// File format (little-endian):
//   [4]   magic "HELX"
//   [2]   version = 1
//   [4]   width
//   [4]   height
//   [2]   flags  (bit 0 = row_debias_used)
//   [W*2] per-row medians (uint16 LE)            -- if flag 0 set
//   [8]   bitstream byte count
//   [N]   Golomb-Rice bitstream (MSB-first)
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

static constexpr int      WIDTH_DEFAULT  = 1500;
static constexpr int      HEIGHT_DEFAULT = 1500;
static constexpr uint8_t  MAGIC[4]       = {'H','E','L','X'};
static constexpr uint16_t VERSION        = 1;

// JPEG-LS-style quantization: 5 levels per gradient, thresholds T1=3, T2=21
// → q ∈ {-2,-1,0,+1,+2}, 5³=125 raw contexts (sparse, 63 canonical after sign)
static constexpr int N_CTX       = 125;
static constexpr int RICE_LIMIT  = 32;  // unary quotient cap before escape
static constexpr int RESET_N     = 64;  // halve context counters when N hits this

// ---------------------------------------------------------------------------
// Bit writer — MSB-first, batched in uint64 buffer
// ---------------------------------------------------------------------------

struct BitWriter {
    std::vector<uint8_t> out;
    uint64_t buf   = 0;
    int      nbits = 0;

    void write(uint64_t value, int n) {
        // assumes n <= 56 and only low n bits of `value` are meaningful
        buf   = (buf << n) | (value & ((n == 64) ? ~0ULL : ((1ULL << n) - 1)));
        nbits += n;
        while (nbits >= 8) {
            nbits -= 8;
            out.push_back((uint8_t)((buf >> nbits) & 0xFF));
            buf &= (nbits == 0 ? 0ULL : ((1ULL << nbits) - 1));
        }
    }

    void write_unary(uint32_t q) {
        // q ones followed by a single zero terminator
        while (q >= 56) {
            write((1ULL << 56) - 1, 56);
            q -= 56;
        }
        if (q > 0) write((1ULL << q) - 1, q);
        write(0, 1);
    }

    void flush() {
        if (nbits > 0) {
            buf <<= (8 - nbits);
            out.push_back((uint8_t)(buf & 0xFF));
            buf = 0; nbits = 0;
        }
    }
};

// ---------------------------------------------------------------------------
// MED predictor (LOCO-I / JPEG-LS)
// Works on int32 because debiased values can be negative.
// ---------------------------------------------------------------------------

static inline int32_t med_predict(int32_t a, int32_t b, int32_t c) {
    int32_t lo = std::min(a, b);
    int32_t hi = std::max(a, b);
    if      (c >= hi) return lo;
    else if (c <= lo) return hi;
    else              return a + b - c;
}

// ---------------------------------------------------------------------------
// Gradient quantization: signed → {-2,-1,0,1,2}
// ---------------------------------------------------------------------------

static inline int quantize_grad(int32_t g) {
    if (g <= -21) return -2;
    if (g <= -3)  return -1;
    if (g <   3)  return  0;
    if (g <  21)  return  1;
    return  2;
}

// Pack canonical (q1,q2,q3) → context index in [0, 124]
static inline int ctx_index(int q1, int q2, int q3) {
    return (q1 + 2) * 25 + (q2 + 2) * 5 + (q3 + 2);
}

// ---------------------------------------------------------------------------
// Per-context adaptive state (Golomb-Rice + bias correction)
// ---------------------------------------------------------------------------

struct CtxState {
    int32_t A;  // sum of |corrected residual|
    int32_t B;  // accumulated bias (running signed residual)
    int32_t N;  // sample count
    int32_t C;  // current integer bias correction
};

static void init_ctx(std::array<CtxState, N_CTX>& s) {
    for (auto& c : s) { c.A = 4; c.B = 0; c.N = 1; c.C = 0; }
}

static inline int rice_k(int32_t A, int32_t N) {
    int k = 0;
    while ((N << k) < A) k++;
    return k;
}

// LOCO-I bias-update (§4.3): adjusts C so the mean residual stays near zero.
static inline void update_bias(CtxState& c, int32_t e_corrected) {
    c.B += e_corrected;
    c.A += (e_corrected < 0) ? -e_corrected : e_corrected;
    c.N += 1;
    if (c.N == RESET_N) { c.A >>= 1; c.B >>= 1; c.N >>= 1; }

    if (c.B <= -c.N) {
        if (c.C > -128) c.C--;
        c.B += c.N;
        if (c.B <= -c.N) c.B = -c.N + 1;
    } else if (c.B > 0) {
        if (c.C <  127) c.C++;
        c.B -= c.N;
        if (c.B > 0) c.B = 0;
    }
}

// ---------------------------------------------------------------------------
// Map signed residual → unsigned (zigzag)
// ---------------------------------------------------------------------------

static inline uint32_t map_signed(int32_t v) {
    return (v >= 0) ? (uint32_t)(2 * v) : (uint32_t)(-2 * v - 1);
}

// ---------------------------------------------------------------------------
// Adaptive Golomb-Rice with escape (LIMIT)
// ---------------------------------------------------------------------------

static inline void rice_encode(BitWriter& bw, uint32_t u, int k) {
    uint32_t q = u >> k;
    if (q < RICE_LIMIT) {
        bw.write_unary(q);
        if (k > 0) bw.write((uint64_t)(u & ((1u << k) - 1)), k);
    } else {
        // Escape: LIMIT ones + zero terminator + 32-bit raw u
        bw.write_unary(RICE_LIMIT);
        bw.write((uint64_t)u, 32);
    }
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
static void write_u64le(FILE* f, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 8, f);
}

// ---------------------------------------------------------------------------
// Row median (raw, exact)
// ---------------------------------------------------------------------------

static uint16_t row_median(const uint16_t* row, int width) {
    std::vector<uint16_t> tmp(row, row + width);
    std::nth_element(tmp.begin(), tmp.begin() + width / 2, tmp.end());
    return tmp[width / 2];
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

    int width  = WIDTH_DEFAULT;
    int height = HEIGHT_DEFAULT;
    int npix   = width * height;

    if (fsize != (long)npix * 2) {
        fprintf(stderr, "Unexpected file size %ld (expected %d)\n", fsize, npix * 2);
        fclose(fin);
        return 1;
    }

    // Read raw big-endian uint16
    std::vector<uint16_t> raw(npix);
    for (int i = 0; i < npix; i++) {
        uint8_t b[2]; fread(b, 1, 2, fin);
        raw[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    // --- Stage 1: per-row debiasing ---
    // Store medians, build int32 debiased image.
    std::vector<uint16_t> med(height);
    std::vector<int32_t>  img(npix);
    for (int y = 0; y < height; y++) {
        med[y] = row_median(raw.data() + y * width, width);
        for (int x = 0; x < width; x++) {
            img[y * width + x] = (int32_t)raw[y * width + x] - (int32_t)med[y];
        }
    }
    // Always enabled — medians are 3 KB header overhead.
    uint16_t flags = 0x0001;

    // --- Stages 2-5: predict + classify + bias-correct + rice-encode ---
    std::array<CtxState, N_CTX> ctx;
    init_ctx(ctx);

    BitWriter bw;
    bw.out.reserve(npix);  // ~1 byte/pixel rough upper bound

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int32_t a = (x > 0)                  ? img[y * width + (x - 1)]       : 0;
            int32_t b = (y > 0)                  ? img[(y - 1) * width + x]       : 0;
            int32_t c = (x > 0 && y > 0)         ? img[(y - 1) * width + (x - 1)] : 0;
            int32_t d = (y > 0 && x + 1 < width) ? img[(y - 1) * width + (x + 1)] : b;

            // Border fallbacks for predictor: replicate available neighbour.
            int32_t pred;
            if (y == 0 && x == 0)      pred = 0;
            else if (y == 0)           pred = a;
            else if (x == 0)           pred = b;
            else                       pred = med_predict(a, b, c);

            // Gradients for context
            int32_t g1 = d - b;
            int32_t g2 = b - c;
            int32_t g3 = c - a;

            int q1 = quantize_grad(g1);
            int q2 = quantize_grad(g2);
            int q3 = quantize_grad(g3);

            // Sign symmetry: canonicalise to leading non-negative gradient
            int sign = 1;
            if (q1 < 0 || (q1 == 0 && q2 < 0) || (q1 == 0 && q2 == 0 && q3 < 0)) {
                sign = -1;
                q1 = -q1; q2 = -q2; q3 = -q3;
            }
            int Q = ctx_index(q1, q2, q3);

            int32_t actual = img[y * width + x];

            // Bias-corrected, sign-flipped residual
            int32_t e_raw = actual - pred;
            int32_t e_s   = sign * e_raw;          // flip sign for canonical ctx
            int32_t e_c   = e_s - ctx[Q].C;        // bias correction

            // Map to unsigned and encode
            uint32_t u = map_signed(e_c);
            int k = rice_k(ctx[Q].A, ctx[Q].N);
            rice_encode(bw, u, k);

            // Update context state with the corrected residual
            update_bias(ctx[Q], e_c);
        }
    }
    bw.flush();

    // --- Write output ---
    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u16le(fout, VERSION);
    write_u32le(fout, (uint32_t)width);
    write_u32le(fout, (uint32_t)height);
    write_u16le(fout, flags);

    // Row medians
    for (int y = 0; y < height; y++) write_u16le(fout, med[y]);

    write_u64le(fout, (uint64_t)bw.out.size());
    fwrite(bw.out.data(), 1, bw.out.size(), fout);

    fclose(fout);

    long out_size = (long)(4 + 2 + 4 + 4 + 2 + height * 2 + 8 + bw.out.size());
    fprintf(stderr, "Compressed: %ld -> %ld bytes (%.4f bpp)\n",
            fsize, out_size, (double)out_size * 8.0 / npix);
    return 0;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input> <output.hlx>\n", argv[0]);
        return 1;
    }
    return compress(argv[1], argv[2]);
}
