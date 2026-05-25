// PRISM Block compressor - rANS variant.
// Same predictor and adaptive model as the range-coder version,
// but uses tabled rANS instead of a range encoder.
//
// rANS state: x in [L, 2L), L = 1<<23.
// Symbol s with freq f, cumul c, total M:
//   encode: x' = (x / f) * M + c + (x % f)
//   decode: s = find(x % M);  x' = f * (x / M) + (x % M) - c
//
// Usage: ./compress_rans <input_raw> <output.prmb> [width height]

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>

static constexpr int     N_LO_CTX     = 8;
static constexpr int     DEFAULT_BLOCK_ROWS = 150;
static constexpr uint8_t MAGIC[4]     = {'P','R','M','1'};

static constexpr int     N_SYMS       = 256;
// rANS parameters
static constexpr uint32_t RANS_L      = 1u << 23;   // lower bound of state
static constexpr uint32_t RANS_M_BITS = 16;          // freq table precision
static constexpr uint32_t RANS_M      = 1u << RANS_M_BITS; // 65536

// ---------------------------------------------------------------------------
// Predictor
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

static inline uint16_t zigzag_enc(uint16_t u) {
    return (u <= 32767u) ? (uint16_t)(u * 2u) : (uint16_t)((65536u - u) * 2u - 1u);
}

// ---------------------------------------------------------------------------
// Adaptive frequency model (Fenwick tree, rescales at 2^16)
// ---------------------------------------------------------------------------

static constexpr uint32_t MAX_TOTAL = 1u << 16;

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

    void update(uint8_t sym) {
        freq[sym]++;
        total++;
        if (total >= MAX_TOTAL) {
            total = 0;
            for (int i = 0; i < N_SYMS; i++) {
                freq[i] = (freq[i] + 1) >> 1;
                total += freq[i];
            }
            rebuild_tree();
            return;
        }
        add_tree(sym, 1);
    }
};

// ---------------------------------------------------------------------------
// rANS encoder (streaming, no pre-scan)
//
// We use a "streaming rANS" approach that works with adaptive models:
// each symbol is encoded with the model state AT THAT SYMBOL.
// Because rANS encodes in reverse we buffer (sym, freq, cumul, total)
// tuples and process them backwards at flush time.
// ---------------------------------------------------------------------------

struct RansEntry {
    uint32_t freq;
    uint32_t cumul;
    uint32_t total;
};

struct RansEncoder {
    std::vector<RansEntry> entries;  // forward order
    std::vector<uint8_t>  out;       // final byte stream (after flush)

    void reserve(size_t n) { entries.reserve(n); }

    // Record the encoding parameters for one symbol (model already updated by caller).
    // freq/cumul/total are the values BEFORE the update (snapshot taken by caller).
    void encode(uint32_t freq, uint32_t cumul, uint32_t total) {
        entries.push_back({freq, cumul, total});
    }

    // Emit the byte stream by processing entries in reverse (rANS property).
    void flush() {
        uint32_t x = RANS_L;
        // We write 32-bit words little-endian; decoder reads them in order.
        // Use a temporary word buffer (reversed at end).
        std::vector<uint32_t> words;
        words.reserve(entries.size() / 2 + 4);

        for (int i = (int)entries.size() - 1; i >= 0; i--) {
            const RansEntry& e = entries[i];
            uint32_t f = e.freq, c = e.cumul, M = e.total;
            // Normalise: push bytes until x is in [L*f/M, 2L*f/M)
            uint32_t x_max = ((RANS_L / M) * f) << 1;  // 2*L*f/M rounded down
            while (x >= x_max) {
                words.push_back(x & 0xFFFFu);
                x >>= 16;
            }
            x = (x / f) * M + c + (x % f);
        }
        // Flush final state (2 x 16-bit words)
        words.push_back(x & 0xFFFFu);
        words.push_back(x >> 16);

        // Reverse so decoder reads forward
        std::reverse(words.begin(), words.end());

        out.reserve(words.size() * 2);
        for (uint32_t w : words) {
            out.push_back((uint8_t)(w & 0xFF));
            out.push_back((uint8_t)(w >> 8));
        }
        entries.clear();
    }
};

// ---------------------------------------------------------------------------
// Block / file I/O helpers
// ---------------------------------------------------------------------------

struct Block {
    uint32_t row0 = 0;
    uint32_t rows = 0;
    std::vector<uint8_t> payload;
};

static bool read_exact(FILE* f, void* data, size_t n) {
    return fread(data, 1, n, f) == n;
}

static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) { b[i] = (uint8_t)(v & 0xFFu); v >>= 8; }
    fwrite(b, 1, 4, f);
}

static void write_u64le(FILE* f, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) { b[i] = (uint8_t)(v & 0xFFu); v >>= 8; }
    fwrite(b, 1, 8, f);
}

// ---------------------------------------------------------------------------
// Context lookup table
// ---------------------------------------------------------------------------

static int LO_CTX_TAB[1024];
static void build_lo_ctx_tab() {
    for (int mg = 0; mg < 1024; mg++) {
        if      (mg == 0)   LO_CTX_TAB[mg] = 0;
        else if (mg <= 2)   LO_CTX_TAB[mg] = 1;
        else if (mg <= 8)   LO_CTX_TAB[mg] = 2;
        else if (mg <= 32)  LO_CTX_TAB[mg] = 3;
        else if (mg <= 64)  LO_CTX_TAB[mg] = 4;
        else if (mg <= 128) LO_CTX_TAB[mg] = 5;
        else if (mg <= 512) LO_CTX_TAB[mg] = 6;
        else                LO_CTX_TAB[mg] = 7;
    }
}

static inline int lo_context(int mW, int mN) {
    int mg = mW + mN;
    return LO_CTX_TAB[mg < 1024 ? mg : 1023];
}

// ---------------------------------------------------------------------------
// Per-pixel encode helper
// ---------------------------------------------------------------------------

static inline void encode_pixel(RansEncoder& enc,
                                 AdaptModel& m_hi, AdaptModel* m_lo,
                                 uint16_t pix, uint16_t pred,
                                 int16_t* curr_res, const int16_t* prev_res,
                                 int gx, bool has_above) {
    uint16_t u   = (uint16_t)(pix - pred);
    uint16_t zz  = zigzag_enc(u);
    uint8_t  hi  = (uint8_t)(zz >> 8);
    uint8_t  lo  = (uint8_t)(zz & 0xFF);
    int16_t  res = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
    curr_res[gx] = res;

    int mW = (gx > 0) ? std::abs((int)curr_res[gx - 1]) : 0;
    int mN = has_above ? std::abs((int)prev_res[gx]) : 0;
    int ctx = lo_context(mW, mN);

    // Snapshot freq/cumul/total before update, then update model.
    uint32_t hi_f = m_hi.freq[hi];
    uint32_t hi_c = m_hi.prefix_less(hi);
    uint32_t hi_t = m_hi.total;
    m_hi.update(hi);
    enc.encode(hi_f, hi_c, hi_t);

    uint32_t lo_f = m_lo[ctx].freq[lo];
    uint32_t lo_c = m_lo[ctx].prefix_less(lo);
    uint32_t lo_t = m_lo[ctx].total;
    m_lo[ctx].update(lo);
    enc.encode(lo_f, lo_c, lo_t);
}

// ---------------------------------------------------------------------------
// Block encoder
// ---------------------------------------------------------------------------

static Block encode_block(const std::vector<uint16_t>& image,
                          int width, int row0, int rows) {
    Block block;
    block.row0 = (uint32_t)row0;
    block.rows = (uint32_t)rows;

    AdaptModel m_hi;
    AdaptModel m_lo[N_LO_CTX];
    m_hi.init();
    for (int c = 0; c < N_LO_CTX; c++) m_lo[c].init();

    RansEncoder enc;
    enc.reserve((size_t)rows * width * 2);

    std::vector<int16_t> prev_res(width, 0), curr_res(width, 0);

    for (int by = 0; by < rows; by++) {
        int gy = row0 + by;
        const uint16_t* row  = image.data() + gy * width;
        const uint16_t* prow = (by > 0) ? (image.data() + (gy - 1) * width) : nullptr;
        bool has_above = (by > 0);

        {
            uint16_t pred = has_above ? prow[0] : 0u;
            encode_pixel(enc, m_hi, m_lo, row[0], pred,
                         curr_res.data(), prev_res.data(), 0, has_above);
        }
        if (has_above) {
            for (int gx = 1; gx < width; gx++) {
                encode_pixel(enc, m_hi, m_lo, row[gx],
                             gap_predict(row[gx-1], prow[gx], prow[gx-1]),
                             curr_res.data(), prev_res.data(), gx, true);
            }
        } else {
            for (int gx = 1; gx < width; gx++) {
                encode_pixel(enc, m_hi, m_lo, row[gx],
                             row[gx - 1],
                             curr_res.data(), prev_res.data(), gx, false);
            }
        }
        prev_res.swap(curr_res);
    }

    enc.flush();
    block.payload = std::move(enc.out);
    return block;
}

// ---------------------------------------------------------------------------
// Dimension inference
// ---------------------------------------------------------------------------

static bool infer_dims(long fsize, int& W, int& H) {
    if (fsize <= 0 || fsize % 2 != 0) return false;
    long npix = fsize / 2;
    long sq = (long)std::sqrt((double)npix);
    if (sq * sq == npix) { W = H = (int)sq; return true; }
    for (int w : {4096, 3000, 2048, 2000, 1920, 1500, 1024, 512, 256}) {
        if (npix % w == 0) { W = w; H = (int)(npix / w); return true; }
    }
    W = (int)npix; H = 1;
    return true;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    build_lo_ctx_tab();
    if (argc != 3 && argc != 5) {
        fprintf(stderr, "Usage: %s <input_raw> <output.prmb> [width height]\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", argv[1]); return 1; }
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    rewind(fin);

    int W, H;
    if (argc == 5) {
        W = std::atoi(argv[3]);
        H = std::atoi(argv[4]);
        if (W <= 0 || H <= 0) { fprintf(stderr, "Invalid dimensions\n"); fclose(fin); return 1; }
        if (fsize != (long)W * H * 2) {
            fprintf(stderr, "File size mismatch: got %ld, expected %d\n", fsize, W * H * 2);
            fclose(fin); return 1;
        }
    } else {
        if (!infer_dims(fsize, W, H)) {
            fprintf(stderr, "Cannot infer dims from size %ld\n", fsize);
            fclose(fin); return 1;
        }
    }

    const int npix = W * H;
    std::vector<uint8_t> input((size_t)fsize);
    if (!read_exact(fin, input.data(), fsize)) {
        fprintf(stderr, "Read failed\n"); fclose(fin); return 1;
    }
    fclose(fin);

    std::vector<uint16_t> image(npix);
    for (int i = 0; i < npix; i++)
        image[i] = (uint16_t)((input[(size_t)i*2] << 8) | input[(size_t)i*2+1]);

    const int block_rows = DEFAULT_BLOCK_ROWS;
    const int nblocks    = (H + block_rows - 1) / block_rows;
    std::vector<Block> blocks(nblocks);

    unsigned nth = std::thread::hardware_concurrency();
    if (nth == 0) nth = 1;
    if (nth > (unsigned)nblocks) nth = (unsigned)nblocks;
    std::atomic<int> next(0);
    std::vector<std::thread> workers;
    workers.reserve(nth);
    for (unsigned t = 0; t < nth; t++) {
        workers.emplace_back([&]() {
            while (true) {
                int b = next.fetch_add(1, std::memory_order_relaxed);
                if (b >= nblocks) break;
                int r0   = b * block_rows;
                int rows = std::min(block_rows, H - r0);
                blocks[b] = encode_block(image, W, r0, rows);
            }
        });
    }
    for (auto& th : workers) th.join();

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", argv[2]); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u32le(fout, (uint32_t)W);
    write_u32le(fout, (uint32_t)H);
    write_u32le(fout, (uint32_t)block_rows);
    write_u32le(fout, (uint32_t)nblocks);
    for (const Block& b : blocks) {
        write_u32le(fout, b.row0);
        write_u32le(fout, b.rows);
        write_u64le(fout, (uint64_t)b.payload.size());
        fwrite(b.payload.data(), 1, b.payload.size(), fout);
    }
    fclose(fout);

    uint64_t out_size = 20;
    for (const Block& b : blocks) out_size += 16 + b.payload.size();
    fprintf(stderr, "Dimensions: %dx%d  blocks=%d threads=%u\n", W, H, nblocks, nth);
    fprintf(stderr, "Compressed: %ld -> %llu bytes (%.4f b/B)\n",
            fsize, (unsigned long long)out_size, (double)out_size * 8.0 / fsize);
    return 0;
}
