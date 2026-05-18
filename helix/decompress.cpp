//
// HELIX decompressor — exact mirror of compress.cpp
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

static constexpr uint8_t  MAGIC[4]   = {'H','E','L','X'};
static constexpr uint16_t VERSION    = 1;
static constexpr int      N_CTX      = 125;
static constexpr int      RICE_LIMIT = 32;
static constexpr int      RESET_N    = 64;

// ---------------------------------------------------------------------------
// Bit reader — MSB-first, batched
// ---------------------------------------------------------------------------

struct BitReader {
    const uint8_t* ptr;
    const uint8_t* end;
    uint64_t buf   = 0;
    int      nbits = 0;

    BitReader(const uint8_t* p, const uint8_t* e) : ptr(p), end(e) {}

    inline void refill() {
        while (nbits <= 56 && ptr < end) {
            buf = (buf << 8) | (uint64_t)(*ptr++);
            nbits += 8;
        }
    }

    inline uint64_t read(int n) {
        if (nbits < n) refill();
        uint64_t v = (n == 64) ? buf
                               : ((buf >> (nbits - n)) & ((1ULL << n) - 1));
        nbits -= n;
        buf &= (nbits == 0 ? 0ULL : ((1ULL << nbits) - 1));
        return v;
    }

    inline uint32_t read_unary() {
        uint32_t total = 0;
        // Loop only when an entire refill is all-ones (extremely rare).
        while (true) {
            if (nbits == 0) refill();
            if (nbits == 0) return total;
            // Align the next bit to bit 63 of a 64-bit word, then count leading 1s
            // as the leading-zero count of the bitwise NOT.
            uint64_t aligned = (nbits == 64) ? buf : (buf << (64 - nbits));
            uint32_t ones = (uint32_t)__builtin_clzll(~aligned);
            if (ones >= (uint32_t)nbits) {
                // All available bits are 1s — consume them, refill and continue.
                total += (uint32_t)nbits;
                nbits = 0; buf = 0;
                continue;
            }
            // Consume `ones` 1-bits + 1 zero terminator.
            int consumed = (int)ones + 1;
            nbits -= consumed;
            buf &= (nbits == 0 ? 0ULL : ((1ULL << nbits) - 1));
            return total + ones;
        }
    }
};

// ---------------------------------------------------------------------------
// Predictor + context utilities (identical to encoder)
// ---------------------------------------------------------------------------

static inline int32_t med_predict(int32_t a, int32_t b, int32_t c) {
    int32_t lo = std::min(a, b);
    int32_t hi = std::max(a, b);
    if      (c >= hi) return lo;
    else if (c <= lo) return hi;
    else              return a + b - c;
}

static inline int quantize_grad(int32_t g) {
    if (g <= -21) return -2;
    if (g <= -3)  return -1;
    if (g <   3)  return  0;
    if (g <  21)  return  1;
    return  2;
}

static inline int ctx_index(int q1, int q2, int q3) {
    return (q1 + 2) * 25 + (q2 + 2) * 5 + (q3 + 2);
}

struct CtxState {
    int32_t A;
    int32_t B;
    int32_t N;
    int32_t C;
};

static void init_ctx(std::array<CtxState, N_CTX>& s) {
    for (auto& c : s) { c.A = 4; c.B = 0; c.N = 1; c.C = 0; }
}

static inline int rice_k(int32_t A, int32_t N) {
    int k = 0;
    while ((N << k) < A) k++;
    return k;
}

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

// inverse of map_signed
static inline int32_t unmap_signed(uint32_t u) {
    return (u & 1) ? -(int32_t)((u + 1) >> 1) : (int32_t)(u >> 1);
}

static inline uint32_t rice_decode(BitReader& br, int k) {
    uint32_t q = br.read_unary();
    if (q < (uint32_t)RICE_LIMIT) {
        uint32_t r = (k > 0) ? (uint32_t)br.read(k) : 0;
        return (q << k) | r;
    } else {
        // Escape: 32 raw bits follow
        return (uint32_t)br.read(32);
    }
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static uint16_t read_u16le(FILE* f) {
    uint8_t b[2]; fread(b, 1, 2, f);
    return (uint16_t)(b[0] | ((uint16_t)b[1] << 8));
}
static uint32_t read_u32le(FILE* f) {
    uint8_t b[4]; fread(b, 1, 4, f);
    uint32_t v = 0;
    for (int i = 3; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}
static uint64_t read_u64le(FILE* f) {
    uint8_t b[8]; fread(b, 1, 8, f);
    uint64_t v = 0;
    for (int i = 7; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.hlx> <output>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open: %s\n", argv[1]); return 1; }

    uint8_t magic[4]; fread(magic, 1, 4, fin);
    if (memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n"); fclose(fin); return 1;
    }
    uint16_t version = read_u16le(fin);
    if (version != VERSION) {
        fprintf(stderr, "Unknown version %u\n", version); fclose(fin); return 1;
    }
    uint32_t width  = read_u32le(fin);
    uint32_t height = read_u32le(fin);
    uint16_t flags  = read_u16le(fin);
    int npix = (int)(width * height);

    std::vector<uint16_t> med(height);
    if (flags & 0x0001) {
        for (uint32_t y = 0; y < height; y++) med[y] = read_u16le(fin);
    }

    uint64_t nbytes = read_u64le(fin);
    std::vector<uint8_t> bits(nbytes);
    fread(bits.data(), 1, nbytes, fin);
    fclose(fin);

    BitReader br(bits.data(), bits.data() + bits.size());

    // --- Decode ---
    std::array<CtxState, N_CTX> ctx;
    init_ctx(ctx);

    std::vector<int32_t> img(npix);

    for (uint32_t y = 0; y < height; y++) {
        for (uint32_t x = 0; x < width; x++) {
            int32_t a = (x > 0)                  ? img[y * width + (x - 1)]       : 0;
            int32_t b = (y > 0)                  ? img[(y - 1) * width + x]       : 0;
            int32_t c = (x > 0 && y > 0)         ? img[(y - 1) * width + (x - 1)] : 0;
            int32_t d = (y > 0 && x + 1 < width) ? img[(y - 1) * width + (x + 1)] : b;

            int32_t pred;
            if (y == 0 && x == 0)      pred = 0;
            else if (y == 0)           pred = a;
            else if (x == 0)           pred = b;
            else                       pred = med_predict(a, b, c);

            int32_t g1 = d - b;
            int32_t g2 = b - c;
            int32_t g3 = c - a;
            int q1 = quantize_grad(g1);
            int q2 = quantize_grad(g2);
            int q3 = quantize_grad(g3);

            int sign = 1;
            if (q1 < 0 || (q1 == 0 && q2 < 0) || (q1 == 0 && q2 == 0 && q3 < 0)) {
                sign = -1;
                q1 = -q1; q2 = -q2; q3 = -q3;
            }
            int Q = ctx_index(q1, q2, q3);

            int k = rice_k(ctx[Q].A, ctx[Q].N);
            uint32_t u  = rice_decode(br, k);
            int32_t  e_c = unmap_signed(u);

            // Reconstruct signed residual: e_s = e_c + C, then e_raw = sign*e_s
            int32_t e_s   = e_c + ctx[Q].C;
            int32_t e_raw = sign * e_s;
            int32_t value = e_raw + pred;
            img[y * width + x] = value;

            update_bias(ctx[Q], e_c);
        }
    }

    // --- Re-add row medians and write big-endian output ---
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", argv[2]); return 1; }

    for (uint32_t y = 0; y < height; y++) {
        int32_t bias = (flags & 0x0001) ? (int32_t)med[y] : 0;
        for (uint32_t x = 0; x < width; x++) {
            int32_t v = img[y * width + x] + bias;
            uint16_t px = (uint16_t)v;
            uint8_t b2[2] = {(uint8_t)(px >> 8), (uint8_t)(px & 0xFF)};
            fwrite(b2, 1, 2, fout);
        }
    }
    fclose(fout);
    return 0;
}
