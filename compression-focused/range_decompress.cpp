#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "../hais2/model/Predictors.hpp"

// Range encoder / decoder — LZMA-style carry propagation.
// Lifted from prism/compress.cpp & prism/decompress.cpp with minor cleanup.

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <algorithm>

namespace apex {

// ---------------------------------------------------------------------------
// Forward declaration — the AdaptModel lives in Fenwick.hpp
// ---------------------------------------------------------------------------

struct AdaptModel;

// ---------------------------------------------------------------------------
// Range encoder (write to byte vector, LZMA carry chain)
// ---------------------------------------------------------------------------

struct RangeEncoder {
    std::vector<uint8_t> out;
    uint64_t low    = 0;
    uint32_t range  = 0xFFFFFFFFu;
    uint8_t  cache  = 0;
    uint32_t pending = 0;

    void shift() {
        bool carry = (low >> 32) != 0;
        uint8_t top = (uint8_t)((uint32_t)low >> 24);
        if (top < 0xFF || carry) {
            out.push_back(cache + (carry ? 1u : 0u));
            uint8_t fill = carry ? 0x00u : 0xFFu;
            for (uint32_t i = 0; i < pending; i++) out.push_back(fill);
            pending = 0;
            cache = top;
        } else {
            pending++;
        }
        low = (uint64_t)((uint32_t)low << 8);
    }

    // Encode one symbol using the given adaptive model.
    // Defined out-of-line after AdaptModel is complete (below, or in Fenwick.hpp).
    inline void encode(AdaptModel& m, uint8_t sym);

    void finish() {
        for (int i = 0; i < 5; i++) shift();
    }
};

// ---------------------------------------------------------------------------
// Range decoder (read from byte buffer, LZMA carry chain)
// ---------------------------------------------------------------------------

struct RangeDecoder {
    const uint8_t* ptr;
    uint32_t range = 0xFFFFFFFFu;
    uint32_t code  = 0;

    void init(const uint8_t* data) {
        ptr = data;
        for (int i = 0; i < 5; i++) code = (code << 8) | (*ptr++);
    }

    // Decode one symbol.  Defined out-of-line after AdaptModel is complete.
    inline uint8_t decode(AdaptModel& m);
};

} // namespace apex
// Adaptive 256-symbol frequency model with Fenwick tree for cumulative lookup.
// Lifted from prism/compress.cpp with minor cleanup.
// Also defines the out-of-line RangeEncoder::encode / RangeDecoder::decode.

#include <cstdint>
#include <algorithm>

namespace apex {

static constexpr uint32_t MAX_TOTAL = 1u << 16;  // halve when total exceeds this
static constexpr int      N_SYMS    = 256;

struct AdaptModel {
    uint32_t freq[N_SYMS];
    uint32_t tree[N_SYMS + 1];
    uint32_t total;

    void init() {
        for (int i = 0; i < N_SYMS; i++) freq[i] = 1;
        total = N_SYMS;
        rebuild_tree();
    }

    void add_tree(int sym, uint32_t delta) {
        for (int i = sym + 1; i <= N_SYMS; i += i & -i) tree[i] += delta;
    }

    void rebuild_tree() {
        for (int i = 0; i <= N_SYMS; i++) tree[i] = 0;
        for (int i = 0; i < N_SYMS; i++) add_tree(i, freq[i]);
    }

    uint32_t prefix_less(int sym) const {
        uint32_t sum = 0;
        for (int i = sym; i > 0; i -= i & -i) sum += tree[i];
        return sum;
    }

    // Find symbol whose cumulative frequency bracket contains `slot`.
    // Also returns cumul (sum of frequencies before the found symbol).
    uint8_t find(uint32_t slot, uint32_t& cumul) const {
        uint32_t idx = 0;
        uint32_t sum = 0;
        for (uint32_t bit = N_SYMS >> 1; bit != 0; bit >>= 1) {
            uint32_t next = idx + bit;
            if (next <= (uint32_t)N_SYMS && sum + tree[next] <= slot) {
                idx = next;
                sum += tree[next];
            }
        }
        cumul = sum;
        return (uint8_t)idx;
    }

    // Encoder-side find (no cumul output needed)
    uint8_t find(uint32_t slot) const {
        uint32_t idx = 0;
        uint32_t sum = 0;
        for (uint32_t bit = N_SYMS >> 1; bit != 0; bit >>= 1) {
            uint32_t next = idx + bit;
            if (next <= (uint32_t)N_SYMS && sum + tree[next] <= slot) {
                idx = next;
                sum += tree[next];
            }
        }
        return (uint8_t)idx;
    }

    void update(uint8_t sym) {
        freq[sym]++;
        total++;
        if (total >= MAX_TOTAL) {
            total = 0;
            for (int i = 0; i < N_SYMS; i++) {
                freq[i] = (freq[i] + 1) >> 1;  // halve, min 1
                total += freq[i];
            }
            rebuild_tree();
            return;
        }
        add_tree(sym, 1);
    }
};

// ---------------------------------------------------------------------------
// Out-of-line definitions for RangeEncoder / RangeDecoder methods
// (they need the full AdaptModel definition)
// ---------------------------------------------------------------------------

inline void RangeEncoder::encode(AdaptModel& m, uint8_t sym) {
    uint32_t r = range / m.total;
    uint32_t cumul = m.prefix_less(sym);
    low += (uint64_t)cumul * r;
    range = (sym < 255) ? m.freq[sym] * r
                        : range - cumul * r;
    while (range < (1u << 24)) { range <<= 8; shift(); }
    m.update(sym);
}

inline uint8_t RangeDecoder::decode(AdaptModel& m) {
    uint32_t r    = range / m.total;
    uint32_t slot = std::min(code / r, m.total - 1);
    uint32_t cumul = 0;
    uint8_t  sym  = m.find(slot, cumul);
    code  -= cumul * r;
    range  = (sym < 255) ? m.freq[sym] * r
                         : range - cumul * r;
    while (range < (1u << 24)) { code = (code << 8) | (*ptr++); range <<= 8; }
    m.update(sym);
    return sym;
}

} // namespace apex

using namespace apex;




static constexpr uint8_t MAGIC[4] = {'A','P','R','1'};

static bool read_exact(FILE* f, void* data, size_t n) {
    return fread(data, 1, n, f) == n;
}
static uint8_t read_u8(FILE* f) {
    int c = fgetc(f);
    return c < 0 ? 0 : (uint8_t)c;
}
static uint16_t read_u16le(FILE* f) {
    uint8_t b[2];
    if (!read_exact(f, b, 2)) return 0;
    return (uint16_t)(b[0] | ((uint16_t)b[1] << 8));
}
static uint32_t read_u32le(FILE* f) {
    uint8_t b[4];
    if (!read_exact(f, b, 4)) return 0;
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}
static float read_float(FILE* f) {
    uint32_t bits = read_u32le(f);
    float v;
    memcpy(&v, &bits, 4);
    return v;
}

static int lo_ctx(uint8_t hi, int mode) {
    if (mode == 1) {
        if (hi == 0) return 0;
        if (hi <= 2) return 1;
        if (hi <= 8) return 2;
        return 3;
    }
    if (mode == 2) return hi;
    return 0;
}

static uint16_t decode_sym(RangeDecoder& dec,
                           std::vector<AdaptModel>& hi_models,
                           std::vector<AdaptModel>& lo_models,
                           int ctx_mode, uint8_t& prev_hi) {
    bool hi_ctx = (ctx_mode >= 3);
    int base_ctx = hi_ctx ? ctx_mode - 3 : ctx_mode;
    int hi_idx = hi_ctx ? lo_ctx(prev_hi, 1) : 0;
    uint8_t hi = dec.decode(hi_models[hi_idx]);
    uint8_t lo = dec.decode(lo_models[lo_ctx(hi, base_ctx)]);
    prev_hi = hi;
    return (uint16_t)(((uint16_t)hi << 8) | lo);
}

static uint16_t clamp_float(float p) {
    if (p < 0.0f) return 0;
    if (p > 65535.0f) return 65535;
    return (uint16_t)(int)(p + 0.5f);
}

static uint16_t predict_pixel(const std::vector<uint16_t>& blk,
                              int bw, int x, int y, int mode,
                              uint16_t gmean,
                              const float* w4,
                              const float* w5,
                              const float* w6,
                              const float* w24) {
    switch (mode) {
    case 1:
        return avg_pred(blk, x, y, bw);
    case 3:
        return gmean;
    case 7:
        return gap_pred(blk, x, y, bw);
    case 8:
        return med_pred(blk, x, y, bw);
    case 4:
        if (y == 0 && x == 0) return clamp_float(w4[3]);
        if (y == 0) return clamp_float(w4[0] * blk[x - 1] + w4[3]);
        if (x == 0) return clamp_float(w4[1] * blk[(y - 1) * bw + x] + w4[3]);
        return clamp_float(w4[0] * blk[y * bw + x - 1] +
                           w4[1] * blk[(y - 1) * bw + x] +
                           w4[2] * blk[(y - 1) * bw + x - 1] +
                           w4[3]);
    case 6:
        if (y == 0 && x == 0) return 0;
        if (y == 0) return blk[x - 1];
        if (x == 0) return blk[(y - 1) * bw + x];
        if (y == 1) {
            float ne = (x < bw - 1) ? (float)blk[(y - 1) * bw + x + 1]
                                    : (float)blk[(y - 1) * bw + x];
            return clamp_float(w5[0] * blk[y * bw + x - 1] +
                               w5[1] * blk[(y - 1) * bw + x] +
                               w5[2] * blk[(y - 1) * bw + x - 1] +
                               w5[3] * ne);
        } else {
            float ne = (x < bw - 1) ? (float)blk[(y - 1) * bw + x + 1]
                                    : (float)blk[(y - 1) * bw + x];
            return clamp_float(w5[0] * blk[y * bw + x - 1] +
                               w5[1] * blk[(y - 1) * bw + x] +
                               w5[2] * blk[(y - 1) * bw + x - 1] +
                               w5[3] * ne +
                               w5[4] * blk[(y - 2) * bw + x]);
        }
    case 9:
        if (y < 2 || x < 2) return gap_pred(blk, x, y, bw);
        return clamp_float(w6[0] * blk[y * bw + x - 1] +
                           w6[1] * blk[(y - 1) * bw + x] +
                           w6[2] * blk[(y - 1) * bw + x - 1] +
                           w6[3] * blk[y * bw + x - 2] +
                           w6[4] * blk[(y - 2) * bw + x] +
                           w6[5]);
    case 10:
        if (y < 3 || x < 4 || x >= bw - 3) return gap_pred(blk, x, y, bw);
        return clamp_float(w6[0] * blk[y * bw + x - 1] +
                           w6[1] * blk[y * bw + x - 2] +
                           w6[2] * blk[y * bw + x - 3] +
                           w6[3] * blk[y * bw + x - 4] +
                           w6[4] * blk[(y - 1) * bw + x] +
                           w6[5] * blk[(y - 1) * bw + x - 1] +
                           w6[6] * blk[(y - 1) * bw + x + 1] +
                           w6[7] * blk[(y - 1) * bw + x - 2] +
                           w6[8] * blk[(y - 1) * bw + x + 2] +
                           w6[9] * blk[(y - 1) * bw + x - 3] +
                           w6[10] * blk[(y - 1) * bw + x + 3] +
                           w6[11] * blk[(y - 2) * bw + x] +
                           w6[12] * blk[(y - 2) * bw + x - 1] +
                           w6[13] * blk[(y - 2) * bw + x + 1] +
                           w6[14] * blk[(y - 3) * bw + x] +
                           w6[15]);
    case 11:
        if (y < 3 || x < 5 || x >= bw - 4) return gap_pred(blk, x, y, bw);
        return clamp_float(w24[0]  * blk[y * bw + x - 1] +
                           w24[1]  * blk[y * bw + x - 2] +
                           w24[2]  * blk[y * bw + x - 3] +
                           w24[3]  * blk[y * bw + x - 4] +
                           w24[4]  * blk[y * bw + x - 5] +
                           w24[5]  * blk[(y - 1) * bw + x] +
                           w24[6]  * blk[(y - 1) * bw + x - 1] +
                           w24[7]  * blk[(y - 1) * bw + x + 1] +
                           w24[8]  * blk[(y - 1) * bw + x - 2] +
                           w24[9]  * blk[(y - 1) * bw + x + 2] +
                           w24[10] * blk[(y - 1) * bw + x - 3] +
                           w24[11] * blk[(y - 1) * bw + x + 3] +
                           w24[12] * blk[(y - 1) * bw + x - 4] +
                           w24[13] * blk[(y - 1) * bw + x + 4] +
                           w24[14] * blk[(y - 2) * bw + x] +
                           w24[15] * blk[(y - 2) * bw + x - 1] +
                           w24[16] * blk[(y - 2) * bw + x + 1] +
                           w24[17] * blk[(y - 2) * bw + x - 2] +
                           w24[18] * blk[(y - 2) * bw + x + 2] +
                           w24[19] * blk[(y - 3) * bw + x] +
                           w24[20] * blk[(y - 3) * bw + x - 1] +
                           w24[21] * blk[(y - 3) * bw + x + 1] +
                           w24[22] * blk[(y - 3) * bw + x - 2] +
                           w24[23]);
    default:
        return 0;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.nvr> <output_raw>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    uint8_t magic[4];
    if (!read_exact(fin, magic, 4) || memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad APR magic\n");
        fclose(fin);
        return 1;
    }
    int W = (int)read_u32le(fin);
    int H = (int)read_u32le(fin);
    int bs = (int)read_u32le(fin);
    uint16_t gmean = read_u16le(fin);
    int blocks_x = (W + bs - 1) / bs;
    int blocks_y = (H + bs - 1) / bs;

    std::vector<uint16_t> image((size_t)W * H);
    for (int by = 0; by < blocks_y; by++) {
        for (int bx = 0; bx < blocks_x; bx++) {
            int bw = std::min(bs, W - bx * bs);
            int bh = std::min(bs, H - by * bs);
            int npix = bw * bh;

            int mode = read_u8(fin);
            int ctx_mode = read_u8(fin);
            uint32_t payload_len = read_u32le(fin);
            float w4[4] = {}, w5[5] = {}, w6[16] = {}, w24[24] = {};
            if (mode == 4) for (float& v : w4) v = read_float(fin);
            else if (mode == 6) for (float& v : w5) v = read_float(fin);
            else if (mode == 9) for (int i = 0; i < 6; i++) w6[i] = read_float(fin);
            else if (mode == 10) for (float& v : w6) v = read_float(fin);
            else if (mode == 11) for (float& v : w24) v = read_float(fin);

            std::vector<uint8_t> payload(payload_len + 16, 0);
            if (!read_exact(fin, payload.data(), payload_len)) {
                fprintf(stderr, "Truncated block payload\n");
                fclose(fin);
                return 1;
            }

            bool hi_ctx = (ctx_mode >= 3);
            int base_ctx = hi_ctx ? ctx_mode - 3 : ctx_mode;
            std::vector<AdaptModel> hi_models(hi_ctx ? 4 : 1);
            for (auto& m : hi_models) m.init();
            std::vector<AdaptModel> lo_models(base_ctx == 0 ? 1 : (base_ctx == 1 ? 4 : 256));
            for (auto& m : lo_models) m.init();
            RangeDecoder dec;
            dec.init(payload.data());

            std::vector<uint16_t> blk(npix);
            uint8_t prev_hi = 0;
            for (int y = 0; y < bh; y++) {
                for (int x = 0; x < bw; x++) {
                    int idx = y * bw + x;
                    uint16_t sym = decode_sym(dec, hi_models, lo_models, ctx_mode, prev_hi);
                    uint16_t pixel;
                    if (mode == 0) {
                        pixel = sym;
                    } else {
                        uint16_t pred = predict_pixel(blk, bw, x, y, mode, gmean, w4, w5, w6, w24);
                        int16_t rc = zagzig(sym);
                        pixel = (uint16_t)((int)pred + rc);
                    }
                    blk[idx] = pixel;
                }
            }

            for (int y = 0; y < bh; y++)
                for (int x = 0; x < bw; x++)
                    image[(by * bs + y) * W + bx * bs + x] = blk[y * bw + x];
        }
    }
    fclose(fin);

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }
    for (uint16_t px : image) {
        uint8_t b[2] = {(uint8_t)(px >> 8), (uint8_t)(px & 255)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
