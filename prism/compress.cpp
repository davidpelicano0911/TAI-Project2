// PRISM v3 — adaptive range coder edition
//
// Same GAP predictor and 8 lo-byte context classes as v2, but the static
// two-pass rANS is replaced by a single-pass adaptive range coder that
// updates its probability model online after each symbol.  No pre-scan, no
// frequency tables stored in the file — smaller header overhead and better
// adaptation to local image statistics.
//
// Adaptive model: each of the 9 contexts (1 hi + 8 lo) maintains 256 symbol
// counts plus a Fenwick tree for cumulative frequencies.  After encoding each
// symbol the count is incremented.  When the total exceeds MAX_TOTAL the whole
// table is halved (floor) with a minimum of 1, keeping the model from going stale.
//
// Format (little-endian):
//   [4]  magic  "PRM3"
//   [4]  width
//   [4]  height
//   [4]  n_lo_ctx  (= N_LO_CTX = 8)
//   [N]  interleaved range-coded stream:
//          for each pixel in raster order:
//            encode hi byte  in hi-byte model
//            encode lo byte  in lo-byte model[ctx]
//
// Usage: ./compress <input_raw> <output.prism>

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

static constexpr int     WIDTH     = 1500;
static constexpr int     HEIGHT    = 1500;
static constexpr int     N_LO_CTX  = 8;
static constexpr int     HEADER_SIZE = 16;

static constexpr uint8_t MAGIC[4]  = {'P','R','M','3'};

// ---------------------------------------------------------------------------
// GAP predictor (W, N, NW weighted blend)
// ---------------------------------------------------------------------------

static inline uint16_t gap_predict(int W, int N, int NW) {
    int dh = std::abs(W - NW);
    int dv = std::abs(N - NW);
    if (dh > dv * 2) return (uint16_t)N;
    if (dv > dh * 2) return (uint16_t)W;
    int wW = dv + 1, wN = dh + 1;
    int pred = (wW * W + wN * N) / (wW + wN);
    int lo = std::min(W, N), hi = std::max(W, N);
    if (pred < lo) pred = lo;
    if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

// ---------------------------------------------------------------------------
// Modular zigzag
// ---------------------------------------------------------------------------

static inline uint16_t zigzag_enc(uint16_t u) {
    return (u <= 32767u) ? (uint16_t)(u * 2u) : (uint16_t)((65536u - u) * 2u - 1u);
}

// ---------------------------------------------------------------------------
// Adaptive frequency model (256 symbols, Fenwick cumulative table)
// ---------------------------------------------------------------------------

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

    uint8_t find(uint32_t slot) const {
        uint32_t idx = 0;
        uint32_t sum = 0;
        for (uint32_t bit = N_SYMS >> 1; bit != 0; bit >>= 1) {
            uint32_t next = idx + bit;
            if (next <= N_SYMS && sum + tree[next] <= slot) {
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
// Range encoder (LZMA-style carry propagation, write to byte vector)
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

    void encode(AdaptModel& m, uint8_t sym) {
        uint32_t r = range / m.total;
        uint32_t cumul = m.prefix_less(sym);
        low += (uint64_t)cumul * r;
        range = (sym < 255) ? m.freq[sym] * r
                            : range - cumul * r;
        while (range < (1u << 24)) { range <<= 8; shift(); }
        m.update(sym);
    }

    void finish() {
        for (int i = 0; i < 5; i++) shift();
    }
};

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4]; for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; } fwrite(b,1,4,f);
}

static bool read_exact(FILE* f, void* data, size_t n) {
    return fread(data, 1, n, f) == n;
}

// ---------------------------------------------------------------------------
// Main compress
// ---------------------------------------------------------------------------

static int compress(const char* in_path, const char* out_path) {
    FILE* fin = fopen(in_path, "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", in_path); return 1; }
    fseek(fin, 0, SEEK_END); long fsize = ftell(fin); rewind(fin);
    int npix = WIDTH * HEIGHT;
    if (fsize != (long)npix * 2) {
        fprintf(stderr, "Unexpected file size %ld\n", fsize); fclose(fin); return 1;
    }
    std::vector<uint8_t> input((size_t)fsize);
    if (!read_exact(fin, input.data(), input.size())) {
        fprintf(stderr, "Input read failed\n");
        fclose(fin);
        return 1;
    }
    fclose(fin);

    std::vector<uint16_t> image(npix);
    for (int i = 0; i < npix; i++) {
        image[i] = (uint16_t)((input[(size_t)i * 2] << 8) | input[(size_t)i * 2 + 1]);
    }

    // Adaptive models: 1 hi + N_LO_CTX lo
    AdaptModel m_hi;
    m_hi.init();
    AdaptModel m_lo[N_LO_CTX];
    for (int c = 0; c < N_LO_CTX; c++) m_lo[c].init();

    RangeEncoder enc;
    std::vector<int16_t> prev_res(WIDTH, 0), curr_res(WIDTH, 0);

    for (int gy = 0; gy < HEIGHT; gy++) {
        for (int gx = 0; gx < WIDTH; gx++) {
            int W  = (gx > 0)           ? image[gy*WIDTH+gx-1]     : 0;
            int N  = (gy > 0)           ? image[(gy-1)*WIDTH+gx]   : W;
            int NW = (gy > 0 && gx > 0) ? image[(gy-1)*WIDTH+gx-1] : W;

            uint16_t pred = (gy == 0 && gx == 0) ? 0u : gap_predict(W, N, NW);
            uint16_t u    = (uint16_t)(image[gy*WIDTH+gx] - pred);
            uint16_t zz   = zigzag_enc(u);
            uint8_t  hi   = (uint8_t)(zz >> 8);
            uint8_t  lo   = (uint8_t)(zz & 0xFF);

            int16_t res = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
            curr_res[gx] = res;

            // Context for lo byte
            int mW = (gx > 0) ? std::abs((int)curr_res[gx-1]) : 0;
            int mN = (gy > 0) ? std::abs((int)prev_res[gx]) : 0;
            int mg = mW + mN;
            int ctx;
            if      (mg ==  0) ctx = 0;
            else if (mg <=  2) ctx = 1;
            else if (mg <=  8) ctx = 2;
            else if (mg <= 32) ctx = 3;
            else if (mg <= 64) ctx = 4;
            else if (mg <= 128) ctx = 5;
            else if (mg <= 512) ctx = 6;
            else                ctx = 7;

            enc.encode(m_hi,     hi);
            enc.encode(m_lo[ctx], lo);
        }
        prev_res.swap(curr_res);
    }
    enc.finish();

    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u32le(fout, (uint32_t)WIDTH);
    write_u32le(fout, (uint32_t)HEIGHT);
    write_u32le(fout, (uint32_t)N_LO_CTX);
    fwrite(enc.out.data(), 1, enc.out.size(), fout);
    fclose(fout);

    long out_size = HEADER_SIZE + (long)enc.out.size();
    fprintf(stderr, "Compressed: %ld -> %ld bytes (%.4f b/B)\n",
            fsize, out_size, (double)out_size * 8.0 / fsize);
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 3) { fprintf(stderr, "Usage: %s <input> <output.prism>\n", argv[0]); return 1; }
    return compress(argv[1], argv[2]);
}
