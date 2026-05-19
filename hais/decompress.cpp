// HAIS — Hybrid Astronomical Image decompressor
// Usage: ./decompress <input.hais> <output>

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>

static constexpr uint8_t  MAGIC[4]   = {'H','A','I','S'};
static constexpr uint32_t SCALE_BITS = 16;
static constexpr uint32_t SCALE      = 1u << SCALE_BITS;

// ---------------------------------------------------------------------------
// FSE / tANS decoder
// ---------------------------------------------------------------------------

struct FseDecodeEntry {
    uint8_t  sym;
    uint8_t  nb_bits;
    uint16_t base;
};

struct FseTable {
    uint32_t freq[256];
    std::vector<FseDecodeEntry> dec;

    void build() {
        dec.assign(SCALE, {});

        uint32_t pos = 0;
        const uint32_t step = (SCALE >> 1) + (SCALE >> 3) + 3;
        std::vector<uint8_t> spread(SCALE);
        for (int s = 0; s < 256; s++) {
            for (uint32_t n = 0; n < freq[s]; n++) {
                spread[pos] = (uint8_t)s;
                pos = (pos + step) & (SCALE - 1);
            }
        }

        uint32_t next[256];
        for (int s = 0; s < 256; s++) next[s] = freq[s];

        for (uint32_t state = 0; state < SCALE; state++) {
            uint8_t sym = spread[state];
            uint32_t x = next[sym]++;
            uint8_t nb = (uint8_t)(SCALE_BITS - (31u - __builtin_clz(x)));
            uint32_t base = (x << nb) - SCALE;
            dec[state] = {sym, nb, (uint16_t)base};
        }
    }
};

struct BitReader {
    const uint8_t* ptr;
    const uint8_t* end;
    uint64_t bitbuf = 0;
    int bitcnt = 0;

    BitReader(const uint8_t* begin, const uint8_t* finish) : ptr(begin), end(finish) {}

    uint32_t get(int n) {
        while (bitcnt < n && ptr < end) {
            bitbuf |= (uint64_t)(*ptr++) << bitcnt;
            bitcnt += 8;
        }
        uint32_t v = (uint32_t)(bitbuf & ((1u << n) - 1u));
        bitbuf >>= n;
        bitcnt -= n;
        return v;
    }
};

struct FseDecoder {
    uint16_t state;
    BitReader br;

    FseDecoder(const uint8_t* bits, const uint8_t* end, uint16_t initial_state)
        : state(initial_state), br(bits, end) {}

    uint8_t next(const FseTable& t) {
        const FseDecodeEntry& e = t.dec[state];
        uint8_t sym = e.sym;
        state = (uint16_t)(e.base + br.get(e.nb_bits));
        return sym;
    }
};

// ---------------------------------------------------------------------------
// Memory-based I/O helpers
// ---------------------------------------------------------------------------

static uint8_t  pget_u8  (const uint8_t*& p) { return *p++; }
static uint32_t pget_u32le(const uint8_t*& p) {
    uint32_t v = (uint32_t)p[0] | ((uint32_t)p[1]<<8) | ((uint32_t)p[2]<<16) | ((uint32_t)p[3]<<24);
    p += 4; return v;
}
static uint64_t pget_u64le(const uint8_t*& p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= (uint64_t)p[i] << (8 * i);
    p += 8; return v;
}

// ---------------------------------------------------------------------------
// Decode one FSE stream from memory pointer (pointer advances past the stream)
// ---------------------------------------------------------------------------

static void decode_stream(const uint8_t*& ptr, std::vector<uint8_t>& out, int n) {
    FseTable tab;
    memset(tab.freq, 0, sizeof(tab.freq));

    uint8_t flag = pget_u8(ptr);
    if (flag == 0xFF) {
        for (int i = 0; i < 256; i++) tab.freq[i] = SCALE / 256;
    } else if (flag == 0) {
        for (int i = 0; i < 256; i++) tab.freq[i] = pget_u32le(ptr);
    } else {
        int nnz = (int)flag;
        for (int k = 0; k < nnz; k++) {
            uint8_t  sym  = pget_u8(ptr);
            uint32_t freq = pget_u32le(ptr);
            tab.freq[sym] = freq;
        }
    }
    tab.build();

    uint64_t nbytes = pget_u64le(ptr);
    const uint8_t* bits = ptr;
    ptr += nbytes;

    out.resize(n);
    uint16_t initial_state = (uint16_t)(bits[0] | ((uint16_t)bits[1] << 8));
    FseDecoder dec(bits + 2, bits + nbytes, initial_state);
    for (int i = 0; i < n; i++) out[i] = dec.next(tab);
}

// ---------------------------------------------------------------------------
// Predictors (inverse)
// ---------------------------------------------------------------------------

static inline int16_t zagzig(uint16_t u) {
    return (u & 1) ? -(int16_t)((u + 1) / 2) : (int16_t)(u / 2);
}

// Mode 1 — smooth average
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
// Per-block decompression task
// ---------------------------------------------------------------------------

struct BlockTask {
    const uint8_t* ptr;  // start of block payload (after mode byte)
    int mode, bw, bh;
};

static void decompress_block(const BlockTask& task, std::vector<uint16_t>& blk) {
    const uint8_t* ptr = task.ptr;
    int npix = task.bw * task.bh;

    float ls_w[3] = {};
    if (task.mode == 2) {
        for (int i = 0; i < 3; i++) {
            uint32_t bits = pget_u32le(ptr);
            memcpy(&ls_w[i], &bits, 4);
        }
    }

    std::vector<uint8_t> hi8, lo8;
    decode_stream(ptr, hi8, npix);
    decode_stream(ptr, lo8, npix);

    blk.resize(npix);
    for (int y = 0; y < task.bh; y++) {
        for (int x = 0; x < task.bw; x++) {
            int idx = y * task.bw + x;
            uint16_t sym = ((uint16_t)hi8[idx] << 8) | lo8[idx];
            uint16_t pixel;
            if (task.mode == 0) {
                pixel = sym;
            } else if (task.mode == 1) {
                pixel = (uint16_t)((int)avg_pred(blk, x, y, task.bw) + zagzig(sym));
            } else {
                uint16_t pred;
                if (y == 0 && x == 0) pred = 0;
                else if (y == 0)      pred = blk[x-1];
                else if (x == 0)      pred = blk[(y-1)*task.bw+x];
                else {
                    float p = ls_w[0]*blk[y*task.bw+(x-1)]
                            + ls_w[1]*blk[(y-1)*task.bw+x]
                            + ls_w[2]*blk[(y-1)*task.bw+(x-1)];
                    pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
                }
                pixel = (uint16_t)((int)pred + zagzig(sym));
            }
            blk[idx] = pixel;
        }
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.hais> <output>\n", argv[0]);
        return 1;
    }

    // Read entire file into memory
    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    rewind(fin);
    std::vector<uint8_t> buf(fsize);
    [[maybe_unused]] auto rb = fread(buf.data(), 1, fsize, fin);
    fclose(fin);

    const uint8_t* p = buf.data();

    if (memcmp(p, MAGIC, 4)) { fprintf(stderr, "Bad magic\n"); return 1; }
    p += 4;
    uint32_t width  = pget_u32le(p);
    uint32_t height = pget_u32le(p);
    uint32_t bs     = pget_u32le(p);

    int blocks_x = ((int)width  + (int)bs - 1) / (int)bs;
    int blocks_y = ((int)height + (int)bs - 1) / (int)bs;
    int total    = blocks_x * blocks_y;

    // Sequential pass: collect block tasks (ptr + metadata per block)
    std::vector<BlockTask> tasks(total);
    std::vector<std::pair<int,int>> block_pos(total); // (bx, by) per block index

    for (int by = 0; by < blocks_y; by++) {
        for (int bx = 0; bx < blocks_x; bx++) {
            int idx = by * blocks_x + bx;
            int bw  = std::min((int)bs, (int)width  - bx * (int)bs);
            int bh  = std::min((int)bs, (int)height - by * (int)bs);
            int mode = (int)pget_u8(p);
            tasks[idx] = {p, mode, bw, bh};
            block_pos[idx] = {bx, by};

            // Skip LS weights if present
            if (mode == 2) p += 12;

            // Skip past the two streams without decoding
            for (int s = 0; s < 2; s++) {
                uint8_t flag = pget_u8(p);
                if (flag == 0xFF) {
                    // no table bytes
                } else if (flag == 0) {
                    p += 256 * 4;
                } else {
                    p += (int)flag * 5;  // nnz × (1 sym + 4 freq)
                }
                uint64_t nbytes = pget_u64le(p);
                p += nbytes;
            }
        }
    }

    // Parallel decompression
    std::vector<std::vector<uint16_t>> blocks(total);
    std::atomic<int> next{0};
    int nthreads = std::max(1, (int)std::thread::hardware_concurrency());

    auto worker = [&]() {
        int idx;
        while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < total)
            decompress_block(tasks[idx], blocks[idx]);
    };
    std::vector<std::thread> pool(nthreads);
    for (auto& t : pool) t = std::thread(worker);
    for (auto& t : pool) t.join();

    // Assemble image
    std::vector<uint16_t> image(width * height);
    for (int idx = 0; idx < total; idx++) {
        auto [bx, by] = block_pos[idx];
        int bw = tasks[idx].bw;
        int bh = tasks[idx].bh;
        const auto& blk = blocks[idx];
        for (int y = 0; y < bh; y++) {
            int gy = by * bs + y;
            for (int x = 0; x < bw; x++)
                image[gy * width + bx * bs + x] = blk[y * bw + x];
        }
    }

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }
    for (uint16_t px : image) {
        uint8_t b[2] = {(uint8_t)(px >> 8), (uint8_t)(px & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
