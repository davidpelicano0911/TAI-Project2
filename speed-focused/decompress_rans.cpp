// PRISM Block decompressor - rANS variant.
//
// Usage: ./decompress_rans <input.prmb> <output_raw>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <vector>

static constexpr int     N_LO_CTX  = 8;
static constexpr uint8_t MAGIC[4]  = {'P','R','M','1'};
static constexpr int     N_SYMS    = 256;
static constexpr uint32_t MAX_TOTAL = 1u << 16;
static constexpr uint32_t RANS_L   = 1u << 23;
static constexpr uint32_t RANS_M_BITS = 16;
static constexpr uint32_t RANS_M   = 1u << RANS_M_BITS;

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

static inline uint16_t zigzag_dec(uint16_t z) {
    return (z & 1u) ? (uint16_t)(65536u - (z + 1u) / 2u) : (uint16_t)(z / 2u);
}

// ---------------------------------------------------------------------------
// Adaptive frequency model
// ---------------------------------------------------------------------------

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

    // Find symbol whose cumulative range contains slot; also return its cumul.
    uint8_t find(uint32_t slot, uint32_t& cumul) const {
        uint32_t idx = 0;
        uint32_t sum = 0;
        for (uint32_t bit = N_SYMS >> 1; bit != 0; bit >>= 1) {
            uint32_t next = idx + bit;
            if (next <= N_SYMS && sum + tree[next] <= slot) {
                idx = next;
                sum += tree[next];
            }
        }
        cumul = sum;
        return (uint8_t)idx;
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
// rANS decoder
// ---------------------------------------------------------------------------

struct RansDecoder {
    const uint16_t* ptr = nullptr;
    const uint16_t* end = nullptr;
    uint32_t x = 0;

    void init(const uint8_t* data, size_t size) {
        ptr = reinterpret_cast<const uint16_t*>(data);
        end = ptr + size / 2;
        // Read initial state from first two 16-bit words.
        uint32_t hi = (ptr < end) ? *ptr++ : 0;
        uint32_t lo = (ptr < end) ? *ptr++ : 0;
        x = (hi << 16) | lo;
    }

    uint8_t decode(AdaptModel& m) {
        uint32_t M    = m.total;
        uint32_t slot = x % M;
        uint32_t cumul = 0;
        uint8_t  sym  = m.find(slot, cumul);
        uint32_t f    = m.freq[sym];
        x = f * (x / M) + slot - cumul;
        // Renormalise: refill until x >= L
        while (x < RANS_L) {
            uint32_t word = (ptr < end) ? *ptr++ : 0;
            x = (x << 16) | word;
        }
        m.update(sym);
        return sym;
    }
};

// ---------------------------------------------------------------------------
// Block / I/O helpers
// ---------------------------------------------------------------------------

struct Block {
    uint32_t row0 = 0;
    uint32_t rows = 0;
    std::vector<uint8_t> payload;
};

static bool read_exact(FILE* f, void* data, size_t n) {
    return fread(data, 1, n, f) == n;
}

static bool read_u32le(FILE* f, uint32_t& v) {
    uint8_t b[4];
    if (!read_exact(f, b, 4)) return false;
    v = 0; for (int i = 3; i >= 0; i--) v = (v << 8) | b[i];
    return true;
}

static bool read_u64le(FILE* f, uint64_t& v) {
    uint8_t b[8];
    if (!read_exact(f, b, 8)) return false;
    v = 0; for (int i = 7; i >= 0; i--) v = (v << 8) | b[i];
    return true;
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
// Block decoder
// ---------------------------------------------------------------------------

static void decode_block(const Block& block, int width, std::vector<uint8_t>& output) {
    // Pad payload so decoder never reads past the end.
    std::vector<uint8_t> stream = block.payload;
    stream.resize(stream.size() + 16, 0);

    AdaptModel m_hi;
    AdaptModel m_lo[N_LO_CTX];
    m_hi.init();
    for (int c = 0; c < N_LO_CTX; c++) m_lo[c].init();

    RansDecoder dec;
    dec.init(stream.data(), stream.size());

    std::vector<uint16_t> prev_img(width, 0), curr_img(width, 0);
    std::vector<int16_t>  prev_res(width, 0), curr_res(width, 0);

    for (uint32_t by = 0; by < block.rows; by++) {
        uint32_t gy = block.row0 + by;
        bool has_above = (by > 0);
        size_t row_off = (size_t)gy * width * 2;

        // First pixel
        {
            int mW = 0;
            int mN = has_above ? std::abs((int)prev_res[0]) : 0;
            uint8_t hi  = dec.decode(m_hi);
            uint8_t lo  = dec.decode(m_lo[lo_context(mW, mN)]);
            uint16_t zz = ((uint16_t)hi << 8) | lo;
            uint16_t u  = zigzag_dec(zz);
            uint16_t pred = has_above ? prev_img[0] : 0u;
            curr_img[0] = (uint16_t)(pred + u);
            curr_res[0] = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
        }

        if (has_above) {
            for (int gx = 1; gx < width; gx++) {
                int mW = std::abs((int)curr_res[gx - 1]);
                int mN = std::abs((int)prev_res[gx]);
                uint8_t hi  = dec.decode(m_hi);
                uint8_t lo  = dec.decode(m_lo[lo_context(mW, mN)]);
                uint16_t zz = ((uint16_t)hi << 8) | lo;
                uint16_t u  = zigzag_dec(zz);
                uint16_t pred = gap_predict(curr_img[gx-1], prev_img[gx], prev_img[gx-1]);
                curr_img[gx] = (uint16_t)(pred + u);
                curr_res[gx] = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
            }
        } else {
            for (int gx = 1; gx < width; gx++) {
                int mW = std::abs((int)curr_res[gx - 1]);
                uint8_t hi  = dec.decode(m_hi);
                uint8_t lo  = dec.decode(m_lo[lo_context(mW, 0)]);
                uint16_t zz = ((uint16_t)hi << 8) | lo;
                uint16_t u  = zigzag_dec(zz);
                curr_img[gx] = (uint16_t)(curr_img[gx - 1] + u);
                curr_res[gx] = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
            }
        }

        for (int gx = 0; gx < width; gx++) {
            uint16_t pix = curr_img[gx];
            output[row_off + gx * 2]     = (uint8_t)(pix >> 8);
            output[row_off + gx * 2 + 1] = (uint8_t)(pix & 0xFF);
        }
        prev_res.swap(curr_res);
        prev_img.swap(curr_img);
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    build_lo_ctx_tab();
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.prmb> <output_raw>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", argv[1]); return 1; }

    uint8_t magic[4];
    if (!read_exact(fin, magic, 4) || memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n"); fclose(fin); return 1;
    }

    uint32_t width = 0, height = 0, block_rows = 0, nblocks = 0;
    if (!read_u32le(fin, width) || !read_u32le(fin, height) ||
        !read_u32le(fin, block_rows) || !read_u32le(fin, nblocks)) {
        fprintf(stderr, "Truncated header\n"); fclose(fin); return 1;
    }
    if (width == 0 || height == 0 || block_rows == 0 || nblocks == 0) {
        fprintf(stderr, "Invalid header\n"); fclose(fin); return 1;
    }

    std::vector<Block> blocks(nblocks);
    for (uint32_t i = 0; i < nblocks; i++) {
        uint64_t size = 0;
        if (!read_u32le(fin, blocks[i].row0) || !read_u32le(fin, blocks[i].rows) ||
            !read_u64le(fin, size)) {
            fprintf(stderr, "Truncated block header\n"); fclose(fin); return 1;
        }
        if (blocks[i].row0 + blocks[i].rows > height || size < 4) {
            fprintf(stderr, "Invalid block\n"); fclose(fin); return 1;
        }
        blocks[i].payload.resize((size_t)size);
        if (!read_exact(fin, blocks[i].payload.data(), blocks[i].payload.size())) {
            fprintf(stderr, "Truncated payload\n"); fclose(fin); return 1;
        }
    }
    fclose(fin);

    std::vector<uint8_t> output((size_t)width * height * 2);

    unsigned nth = std::thread::hardware_concurrency();
    if (nth == 0) nth = 1;
    if (nth > nblocks) nth = nblocks;
    std::atomic<uint32_t> next(0);
    std::vector<std::thread> workers;
    workers.reserve(nth);
    for (unsigned t = 0; t < nth; t++) {
        workers.emplace_back([&]() {
            while (true) {
                uint32_t b = next.fetch_add(1, std::memory_order_relaxed);
                if (b >= nblocks) break;
                decode_block(blocks[b], (int)width, output);
            }
        });
    }
    for (auto& th : workers) th.join();

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", argv[2]); return 1; }
    fwrite(output.data(), 1, output.size(), fout);
    fclose(fout);
    return 0;
}
