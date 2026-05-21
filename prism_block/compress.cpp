// PRISM Block compressor - independent threaded row blocks.
//
// Usage: ./compress <input_raw> <output.prismb>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>

static constexpr int N_LO_CTX = 8;
static constexpr int DEFAULT_BLOCK_ROWS = 150;
static constexpr uint8_t MAGIC[4] = {'P','R','B','1'};

static constexpr uint32_t MAX_TOTAL = 1u << 16;
static constexpr int N_SYMS = 256;

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

struct RangeEncoder {
    std::vector<uint8_t> out;
    uint64_t low = 0;
    uint32_t range = 0xFFFFFFFFu;
    uint8_t cache = 0;
    uint32_t pending = 0;

    void shift() {
        bool carry = (low >> 32) != 0;
        uint8_t top = (uint8_t)((uint32_t)low >> 24);
        if (top < 0xFF || carry) {
            out.push_back((uint8_t)(cache + (carry ? 1u : 0u)));
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
        uint32_t r     = range / m.total;
        uint32_t cumul = m.prefix_less(sym);
        low  += (uint64_t)cumul * r;
        range = (sym < 255) ? m.freq[sym] * r : range - cumul * r;
        while (range < (1u << 24)) { range <<= 8; shift(); }
        m.update(sym);
    }

    void finish() {
        for (int i = 0; i < 5; i++) shift();
    }
};

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

// Lookup table: ctx = LO_CTX_TAB[min(mW+mN, 1023)]
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

static inline void encode_pixel(RangeEncoder& enc, AdaptModel& m_hi,
                                 AdaptModel* m_lo,
                                 uint16_t pix, uint16_t pred,
                                 int16_t* curr_res, const int16_t* prev_res,
                                 int gx, bool has_above) {
    uint16_t u  = (uint16_t)(pix - pred);
    uint16_t zz = zigzag_enc(u);
    uint8_t  hi = (uint8_t)(zz >> 8);
    uint8_t  lo = (uint8_t)(zz & 0xFF);
    int16_t  res = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
    curr_res[gx] = res;
    int mW = (gx > 0) ? std::abs((int)curr_res[gx - 1]) : 0;
    int mN = has_above ? std::abs((int)prev_res[gx]) : 0;
    enc.encode(m_hi, hi);
    enc.encode(m_lo[lo_context(mW, mN)], lo);
}

static Block encode_block(const std::vector<uint16_t>& image,
                          int width, int row0, int rows) {
    Block block;
    block.row0 = (uint32_t)row0;
    block.rows = (uint32_t)rows;

    AdaptModel m_hi;
    AdaptModel m_lo[N_LO_CTX];
    m_hi.init();
    for (int c = 0; c < N_LO_CTX; c++) m_lo[c].init();

    RangeEncoder enc;
    enc.out.reserve((size_t)rows * width * 5 / 4);

    std::vector<int16_t> prev_res(width, 0), curr_res(width, 0);

    for (int by = 0; by < rows; by++) {
        int gy = row0 + by;
        const uint16_t* row  = image.data() + gy * width;
        const uint16_t* prow = (by > 0) ? (image.data() + (gy - 1) * width) : nullptr;
        bool has_above = (by > 0);

        // First pixel of row: no left neighbor
        {
            int N = has_above ? (int)prow[0] : 0;
            uint16_t pred = (by == 0) ? 0u : (uint16_t)N;
            encode_pixel(enc, m_hi, m_lo, row[0], pred,
                         curr_res.data(), prev_res.data(), 0, has_above);
        }
        // Rest of row: left neighbor always available
        if (has_above) {
            for (int gx = 1; gx < width; gx++) {
                int W  = (int)row[gx - 1];
                int N  = (int)prow[gx];
                int NW = (int)prow[gx - 1];
                encode_pixel(enc, m_hi, m_lo, row[gx],
                             gap_predict(W, N, NW),
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

    enc.finish();
    block.payload = std::move(enc.out);
    return block;
}

// Infer dimensions from file size: try square root, then common widths, then Wx1.
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

int main(int argc, char* argv[]) {
    build_lo_ctx_tab();
    if (argc != 3 && argc != 5) {
        fprintf(stderr, "Usage: %s <input_raw> <output.prismb> [width height]\n", argv[0]);
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
        if (W <= 0 || H <= 0) {
            fprintf(stderr, "Invalid dimensions\n"); fclose(fin); return 1;
        }
        if (fsize != (long)W * H * 2) {
            fprintf(stderr, "File size mismatch: got %ld, expected %d\n", fsize, W * H * 2);
            fclose(fin); return 1;
        }
    } else {
        if (!infer_dims(fsize, W, H)) {
            fprintf(stderr, "Cannot infer dimensions from file size %ld\n", fsize);
            fclose(fin); return 1;
        }
    }

    const int npix = W * H;
    std::vector<uint8_t> input((size_t)fsize);
    if (!read_exact(fin, input.data(), input.size())) {
        fprintf(stderr, "Input read failed\n"); fclose(fin); return 1;
    }
    fclose(fin);

    std::vector<uint16_t> image(npix);
    for (int i = 0; i < npix; i++)
        image[i] = (uint16_t)((input[(size_t)i * 2] << 8) | input[(size_t)i * 2 + 1]);

    const int block_rows = DEFAULT_BLOCK_ROWS;
    const int nblocks = (H + block_rows - 1) / block_rows;
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
                int row0 = b * block_rows;
                int rows = std::min(block_rows, H - row0);
                blocks[b] = encode_block(image, W, row0, rows);
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
