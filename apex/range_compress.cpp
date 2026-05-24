#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "../hais2/model/CompressorModel.hpp"

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

static void write_u8(FILE* f, uint8_t v) { fwrite(&v, 1, 1, f); }
static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 255), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) { b[i] = (uint8_t)(v & 255); v >>= 8; }
    fwrite(b, 1, 4, f);
}
static void put_float(std::vector<uint8_t>& out, float v) {
    uint32_t bits;
    memcpy(&bits, &v, 4);
    for (int i = 0; i < 4; i++) out.push_back((uint8_t)((bits >> (8 * i)) & 255));
}

static bool infer_dims(long fsize, int& W, int& H) {
    if (fsize <= 0 || (fsize & 1)) return false;
    long npix = fsize / 2;
    long sq = (long)std::sqrt((double)npix);
    if (sq * sq == npix) { W = H = (int)sq; return true; }
    if (npix % 1500 == 0) { W = 1500; H = (int)(npix / 1500); return true; }
    return false;
}

static int default_block_size(int W, int H) {
    for (int bs = 512; bs >= 32; --bs)
        if (W % bs == 0 && H % bs == 0 && (W / bs) * (H / bs) >= 4)
            return bs;
    return 500;
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

static std::vector<uint8_t> encode_syms(const std::vector<uint16_t>& syms,
                                        int ctx_mode) {
    std::vector<AdaptModel> hi_models(1);
    std::vector<AdaptModel> lo_models(ctx_mode == 0 ? 1 : (ctx_mode == 1 ? 4 : 256));
    hi_models[0].init();
    for (auto& m : lo_models) m.init();

    RangeEncoder enc;
    enc.out.reserve(syms.size() * 2 / 3);
    for (uint16_t s : syms) {
        uint8_t hi = (uint8_t)(s >> 8);
        uint8_t lo = (uint8_t)(s & 255);
        enc.encode(hi_models[0], hi);
        enc.encode(lo_models[lo_ctx(hi, ctx_mode)], lo);
    }
    enc.finish();
    return enc.out;
}

struct BlockResult {
    uint8_t mode = 0;
    uint8_t ctx_mode = 0;
    std::vector<uint8_t> weights;
    std::vector<uint8_t> payload;
};

static void add_weights(std::vector<uint8_t>& out, const float* w, int n) {
    out.clear();
    for (int i = 0; i < n; i++) put_float(out, w[i]);
}

static uint16_t clamp_float(float p) {
    if (p < 0.0f) return 0;
    if (p > 65535.0f) return 65535;
    return (uint16_t)(int)(p + 0.5f);
}

static void ls16_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                         float w[16], std::vector<uint16_t>& syms) {
    double XtX[16][16] = {}, Xty[16] = {};
    for (int y = 3; y < bh; y++) {
        for (int x = 4; x < bw - 3; x++) {
            double f[16] = {
                (double)blk[y * bw + x - 1],
                (double)blk[y * bw + x - 2],
                (double)blk[y * bw + x - 3],
                (double)blk[y * bw + x - 4],
                (double)blk[(y - 1) * bw + x],
                (double)blk[(y - 1) * bw + x - 1],
                (double)blk[(y - 1) * bw + x + 1],
                (double)blk[(y - 1) * bw + x - 2],
                (double)blk[(y - 1) * bw + x + 2],
                (double)blk[(y - 1) * bw + x - 3],
                (double)blk[(y - 1) * bw + x + 3],
                (double)blk[(y - 2) * bw + x],
                (double)blk[(y - 2) * bw + x - 1],
                (double)blk[(y - 2) * bw + x + 1],
                (double)blk[(y - 3) * bw + x],
                1.0
            };
            double t = blk[y * bw + x];
            for (int i = 0; i < 16; i++) {
                Xty[i] += f[i] * t;
                for (int j = 0; j < 16; j++) XtX[i][j] += f[i] * f[j];
            }
        }
    }
    ls_solve<16>(XtX, Xty, w);

    syms.resize((size_t)bw * bh);
    for (int y = 0; y < bh; y++) {
        for (int x = 0; x < bw; x++) {
            int idx = y * bw + x;
            uint16_t pred;
            if (y < 3 || x < 4 || x >= bw - 3) {
                pred = gap_pred(blk, x, y, bw);
            } else {
                float p = w[0] * blk[y * bw + x - 1]
                        + w[1] * blk[y * bw + x - 2]
                        + w[2] * blk[y * bw + x - 3]
                        + w[3] * blk[y * bw + x - 4]
                        + w[4] * blk[(y - 1) * bw + x]
                        + w[5] * blk[(y - 1) * bw + x - 1]
                        + w[6] * blk[(y - 1) * bw + x + 1]
                        + w[7] * blk[(y - 1) * bw + x - 2]
                        + w[8] * blk[(y - 1) * bw + x + 2]
                        + w[9] * blk[(y - 1) * bw + x - 3]
                        + w[10] * blk[(y - 1) * bw + x + 3]
                        + w[11] * blk[(y - 2) * bw + x]
                        + w[12] * blk[(y - 2) * bw + x - 1]
                        + w[13] * blk[(y - 2) * bw + x + 1]
                        + w[14] * blk[(y - 3) * bw + x]
                        + w[15];
                pred = clamp_float(p);
            }
            syms[idx] = zigzag((int16_t)(blk[idx] - pred));
        }
    }
}

static BlockResult compress_block(const std::vector<uint16_t>& image,
                                  int W, int H, int bs,
                                  int bx, int by, uint16_t gmean) {
    int bw = std::min(bs, W - bx * bs);
    int bh = std::min(bs, H - by * bs);
    int npix = bw * bh;

    std::vector<uint16_t> blk(npix);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++)
            blk[y * bw + x] = image[(by * bs + y) * W + bx * bs + x];

    std::vector<uint16_t> s0, s1, s3, s4, s6, s9, s10;
    float w4[4], w5[5], w6[6], w16[16];
    make_syms(blk, 0, bw, bh, s0, gmean);
    make_syms(blk, 1, bw, bh, s1, gmean);
    make_syms(blk, 3, bw, bh, s3, gmean);
    ls4_compute(blk, bw, bh, w4, s4);
    ls5_compute(blk, bw, bh, w5, s6);
    ls6_compute(blk, bw, bh, w6, s9);
    ls16_compute(blk, bw, bh, w16, s10);

    struct Candidate {
        uint8_t mode;
        const std::vector<uint16_t>* syms;
        const float* weights;
        int nweights;
    };
    const Candidate candidates[] = {
        {0, &s0, nullptr, 0},
        {1, &s1, nullptr, 0},
        {3, &s3, nullptr, 0},
        {4, &s4, w4, 4},
        {6, &s6, w5, 5},
        {9, &s9, w6, 6},
        {10, &s10, w16, 16},
    };

    BlockResult best;
    size_t best_size = (size_t)-1;
    for (const auto& c : candidates) {
        for (int ctx = 0; ctx <= 2; ctx++) {
            std::vector<uint8_t> payload = encode_syms(*c.syms, ctx);
            size_t total = payload.size() + (size_t)c.nweights * 4;
            if (total < best_size) {
                best_size = total;
                best.mode = c.mode;
                best.ctx_mode = (uint8_t)ctx;
                if (c.nweights) add_weights(best.weights, c.weights, c.nweights);
                else best.weights.clear();
                best.payload = std::move(payload);
            }
        }
    }
    return best;
}

int main(int argc, char* argv[]) {
    // Suporta duas formas:
    //   range_compress n_rows n_cols input output [block_size]  (argc 5 ou 6)
    //   range_compress input output [block_size]                (argc 3 ou 4, backward compat)
    if (argc < 3 || argc > 6 || argc == 2) {
        fprintf(stderr, "Usage: %s [n_rows n_cols] <input_raw> <output.nvr> [block_size]\n", argv[0]);
        return 1;
    }

    bool has_dims = (argc == 5 || argc == 6) && std::atoi(argv[1]) > 0;
    int  n_rows_arg = has_dims ? std::atoi(argv[1]) : 0;
    int  n_cols_arg = has_dims ? std::atoi(argv[2]) : 0;
    const char* src_arg = has_dims ? argv[3] : argv[1];
    const char* dst_arg = has_dims ? argv[4] : argv[2];
    int  bs_arg = (has_dims && argc == 6) ? std::atoi(argv[5]) :
                  (!has_dims && argc == 4) ? std::atoi(argv[3]) : 0;

    FILE* fin = fopen(src_arg, "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", src_arg); return 1; }
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    rewind(fin);

    int W = 0, H = 0;
    if (has_dims) {
        H = n_rows_arg;
        W = n_cols_arg;
        if (W <= 0 || H <= 0 || fsize != (long)W * H * 2) {
            fprintf(stderr, "Dimensions %dx%d don't match file size %ld\n", W, H, fsize);
            fclose(fin);
            return 1;
        }
    } else if (!infer_dims(fsize, W, H)) {
        fprintf(stderr, "Cannot infer image dimensions — use: %s n_rows n_cols ...\n", argv[0]);
        fclose(fin);
        return 1;
    }

    std::vector<uint16_t> image((size_t)W * H);
    for (int i = 0; i < W * H; i++) {
        int a = fgetc(fin), b = fgetc(fin);
        if (a < 0 || b < 0) { fclose(fin); return 1; }
        image[i] = (uint16_t)((a << 8) | b);
    }
    fclose(fin);

    uint64_t sum = 0;
    for (uint16_t px : image) sum += px;
    uint16_t gmean = (uint16_t)(sum / image.size());

    int bs = (bs_arg > 0) ? bs_arg : default_block_size(W, H);
    if (bs <= 0) bs = default_block_size(W, H);
    int blocks_x = (W + bs - 1) / bs;
    int blocks_y = (H + bs - 1) / bs;

    FILE* fout = fopen(dst_arg, "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", dst_arg); return 1; }
    fwrite(MAGIC, 1, 4, fout);
    write_u32le(fout, (uint32_t)W);
    write_u32le(fout, (uint32_t)H);
    write_u32le(fout, (uint32_t)bs);
    write_u16le(fout, gmean);

    for (int by = 0; by < blocks_y; by++) {
        for (int bx = 0; bx < blocks_x; bx++) {
            BlockResult r = compress_block(image, W, H, bs, bx, by, gmean);
            write_u8(fout, r.mode);
            write_u8(fout, r.ctx_mode);
            write_u32le(fout, (uint32_t)r.payload.size());
            fwrite(r.weights.data(), 1, r.weights.size(), fout);
            fwrite(r.payload.data(), 1, r.payload.size(), fout);
        }
    }
    fclose(fout);
    return 0;
}
