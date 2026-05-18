// HAIS — Hybrid Astronomical Image compressor
//
// Pipeline por bloco 150×150:
//   1. Gerar símbolos com 3 modos: raw / left / MED
//   2. Escolher modo com menor entropia
//   3. Comprimir com ANS (rANS com split hi8+lo8)
//
// Formato do bloco:
//   [1]      modo (0=raw, 1=left, 2=MED)
//   [256*4]  tabela freq hi8
//   [8]      tamanho stream hi8
//   [N]      stream hi8
//   [256*4]  tabela freq lo8
//   [8]      tamanho stream lo8
//   [M]      stream lo8
//
// Usage: ./compress <input> <output.hais>

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>

static constexpr int     WIDTH      = 1500;
static constexpr int     HEIGHT     = 1500;
static constexpr int     BLOCK_SIZE = 150;
static constexpr int     BLOCKS_X   = WIDTH  / BLOCK_SIZE;
static constexpr int     BLOCKS_Y   = HEIGHT / BLOCK_SIZE;
static constexpr int     BLOCK_NPIX = BLOCK_SIZE * BLOCK_SIZE;

static constexpr uint8_t  MAGIC[4] = {'H','A','I','S'};
static constexpr uint16_t VERSION  = 1;

// ---------------------------------------------------------------------------
// rANS (ANS com emissão de bytes)
// ---------------------------------------------------------------------------

static constexpr uint32_t RANS_SCALE_BITS = 16;
static constexpr uint32_t RANS_SCALE      = 1u << RANS_SCALE_BITS;
static constexpr uint32_t RANS_L          = 1u << 23;

struct RansTable {
    uint32_t freq[256];
    uint32_t cumul[257];

    void build_cumul() {
        cumul[0] = 0;
        for (int i = 0; i < 256; i++) cumul[i+1] = cumul[i] + freq[i];
    }
};

static void normalise_freqs(const uint64_t* counts, uint32_t* freq) {
    uint64_t total = 0;
    for (int i = 0; i < 256; i++) total += counts[i];

    uint32_t assigned = 0;
    for (int i = 0; i < 256; i++) {
        if (counts[i] == 0) { freq[i] = 0; continue; }
        freq[i] = (uint32_t)std::max((uint64_t)1,
                    (counts[i] * (uint64_t)RANS_SCALE) / total);
        assigned += freq[i];
    }
    int best = 0;
    for (int i = 1; i < 256; i++)
        if (counts[i] > counts[best]) best = i;
    if (assigned < RANS_SCALE)      freq[best] += RANS_SCALE - assigned;
    else if (assigned > RANS_SCALE) freq[best] -= assigned - RANS_SCALE;
}

struct RansEncoder {
    uint32_t state = RANS_L;
    std::vector<uint8_t> buf;

    void encode(uint8_t sym, const RansTable& tab) {
        uint32_t f = tab.freq[sym];
        uint32_t c = tab.cumul[sym];
        uint32_t x_max = ((RANS_L / RANS_SCALE) * 256) * f;
        while (state >= x_max) { buf.push_back(state & 0xFF); state >>= 8; }
        state = (state / f) * RANS_SCALE + c + (state % f);
    }

    void flush() {
        buf.push_back( state        & 0xFF);
        buf.push_back((state >>  8) & 0xFF);
        buf.push_back((state >> 16) & 0xFF);
        buf.push_back((state >> 24) & 0xFF);
        std::reverse(buf.begin(), buf.end());
    }
};

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static void write_u8   (FILE* f, uint8_t  v) { fwrite(&v, 1, 1, f); }
static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4]; for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 4, f);
}

// ---------------------------------------------------------------------------
// Zigzag (usa int16 — residuos fora de [-32768,32767] são raros e só ocorrem
// em blocos onde mode 0 vence, logo nunca chegam à codificação)
// ---------------------------------------------------------------------------

static inline uint16_t zigzag(int16_t v) {
    return (uint16_t)((v >= 0) ? (v * 2) : ((-v) * 2 - 1));
}

// ---------------------------------------------------------------------------
// Preditor MED (JPEG-LS)
// ---------------------------------------------------------------------------

static inline uint16_t med_predict(const std::vector<uint16_t>& block,
                                    int x, int y) {
    if (y == 0 && x == 0) return 0;
    if (y == 0)           return block[x - 1];
    if (x == 0)           return block[(y - 1) * BLOCK_SIZE + x];
    uint16_t A = block[y       * BLOCK_SIZE + (x - 1)];
    uint16_t B = block[(y - 1) * BLOCK_SIZE +  x     ];
    uint16_t C = block[(y - 1) * BLOCK_SIZE + (x - 1)];
    int pred = (int)A + (int)B - (int)C;
    int lo   = std::min((int)A, (int)B);
    int hi   = std::max((int)A, (int)B);
    if (pred < lo) pred = lo;
    if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

// ---------------------------------------------------------------------------
// Geração de símbolos (uint16) para cada modo
// ---------------------------------------------------------------------------

static std::vector<uint16_t> make_raw(const std::vector<uint16_t>& block) {
    return block;
}

static std::vector<uint16_t> make_left(const std::vector<uint16_t>& block) {
    std::vector<uint16_t> syms(BLOCK_NPIX);
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            uint16_t pred;
            if      (x > 0) pred = block[y * BLOCK_SIZE + (x - 1)];
            else if (y > 0) pred = block[(y - 1) * BLOCK_SIZE + x];
            else            pred = 0;
            syms[y * BLOCK_SIZE + x] = zigzag((int16_t)(block[y * BLOCK_SIZE + x] - pred));
        }
    }
    return syms;
}

static std::vector<uint16_t> make_med(const std::vector<uint16_t>& block) {
    std::vector<uint16_t> syms(BLOCK_NPIX);
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            uint16_t pred = med_predict(block, x, y);
            syms[y * BLOCK_SIZE + x] = zigzag((int16_t)(block[y * BLOCK_SIZE + x] - pred));
        }
    }
    return syms;
}

// Mode 3: smooth causal predictor — uniform average of 3 causal neighbours
// (Quintas-Torra et al. 2026, PASP 138:044505, Section 3.1 "2×2 uniform")
// pred = round((A + B + C) / 3)  with A=left, B=above, C=above-left
static std::vector<uint16_t> make_smooth(const std::vector<uint16_t>& block) {
    std::vector<uint16_t> syms(BLOCK_NPIX);
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            uint16_t pred;
            if (y == 0 && x == 0)   pred = 0;
            else if (y == 0)         pred = block[x - 1];
            else if (x == 0)         pred = block[(y - 1) * BLOCK_SIZE + x];
            else {
                uint16_t A = block[y       * BLOCK_SIZE + (x - 1)];
                uint16_t B = block[(y - 1) * BLOCK_SIZE +  x     ];
                uint16_t C = block[(y - 1) * BLOCK_SIZE + (x - 1)];
                pred = (uint16_t)(((int)A + (int)B + (int)C + 1) / 3);
            }
            syms[y * BLOCK_SIZE + x] = zigzag((int16_t)(block[y * BLOCK_SIZE + x] - pred));
        }
    }
    return syms;
}

// ---------------------------------------------------------------------------
// Entropia empírica (bits/símbolo, sobre uint16)
// ---------------------------------------------------------------------------

static double entropy(const std::vector<uint16_t>& syms) {
    std::vector<uint64_t> freq(65536, 0);
    for (uint16_t s : syms) freq[s]++;
    double H = 0.0, N = (double)syms.size();
    for (int i = 0; i < 65536; i++) {
        if (freq[i] == 0) continue;
        double p = freq[i] / N;
        H -= p * std::log2(p);
    }
    return H;
}

// ---------------------------------------------------------------------------
// encode_stream into a memory buffer (thread-safe, no FILE* needed)
// ---------------------------------------------------------------------------

static void buf_u8   (std::vector<uint8_t>& b, uint8_t  v) { b.push_back(v); }
static void buf_u32le(std::vector<uint8_t>& b, uint32_t v) {
    for (int i = 0; i < 4; i++) { b.push_back(v & 0xFF); v >>= 8; }
}
static void buf_u64le(std::vector<uint8_t>& b, uint64_t v) {
    for (int i = 0; i < 8; i++) { b.push_back(v & 0xFF); v >>= 8; }
}

static void encode_stream_buf(std::vector<uint8_t>& out, const std::vector<uint8_t>& bytes) {
    uint64_t cnt[256] = {};
    for (uint8_t byte : bytes) cnt[byte]++;

    int nnz = 0;
    for (int i = 0; i < 256; i++) if (cnt[i] > 0) nnz++;

    RansTable tab;
    normalise_freqs(cnt, tab.freq);
    tab.build_cumul();

    RansEncoder enc;
    for (int i = (int)bytes.size() - 1; i >= 0; i--)
        enc.encode(bytes[i], tab);
    enc.flush();

    if (nnz < 256) {
        buf_u8(out, (uint8_t)nnz);
        for (int i = 0; i < 256; i++) {
            if (tab.freq[i] == 0) continue;
            buf_u8(out, (uint8_t)i);
            buf_u32le(out, tab.freq[i]);
        }
    } else {
        buf_u8(out, 0);
        for (int i = 0; i < 256; i++) buf_u32le(out, tab.freq[i]);
    }

    buf_u64le(out, (uint64_t)enc.buf.size());
    out.insert(out.end(), enc.buf.begin(), enc.buf.end());
}

// ---------------------------------------------------------------------------
// Compress one block into a byte buffer — pure function, safe to run in threads
// ---------------------------------------------------------------------------

struct BlockResult {
    uint8_t mode;
    std::vector<uint8_t> payload;
};

static BlockResult compress_block(const std::vector<uint16_t>& image, int bx, int by) {
    std::vector<uint16_t> block(BLOCK_NPIX);
    for (int y = 0; y < BLOCK_SIZE; y++) {
        int gy = by * BLOCK_SIZE + y;
        for (int x = 0; x < BLOCK_SIZE; x++)
            block[y * BLOCK_SIZE + x] = image[gy * WIDTH + bx * BLOCK_SIZE + x];
    }

    auto s0 = make_raw   (block);
    auto s1 = make_left  (block);
    auto s2 = make_med   (block);
    auto s3 = make_smooth(block);

    double H[4] = { entropy(s0), entropy(s1), entropy(s2), entropy(s3) };
    int best = 0;
    if (H[1] < H[best]) best = 1;
    if (H[2] < H[best]) best = 2;
    if (H[3] < H[best]) best = 3;

    const auto& syms = (best == 0) ? s0 : (best == 1) ? s1 : (best == 2) ? s2 : s3;
    std::vector<uint8_t> hi8(BLOCK_NPIX), lo8(BLOCK_NPIX);
    for (int i = 0; i < BLOCK_NPIX; i++) {
        hi8[i] = (uint8_t)(syms[i] >> 8);
        lo8[i] = (uint8_t)(syms[i] & 0xFF);
    }

    BlockResult res;
    res.mode = (uint8_t)best;
    encode_stream_buf(res.payload, hi8);
    encode_stream_buf(res.payload, lo8);
    return res;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input> <output.hais>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    std::vector<uint16_t> image(WIDTH * HEIGHT);
    for (int i = 0; i < WIDTH * HEIGHT; i++) {
        uint8_t b[2];
        if (fread(b, 1, 2, fin) != 2) { fprintf(stderr, "Short read\n"); return 1; }
        image[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u16le(fout, VERSION);
    write_u32le(fout, (uint32_t)WIDTH);
    write_u32le(fout, (uint32_t)HEIGHT);
    write_u32le(fout, (uint32_t)BLOCK_SIZE);

    int total_blocks = BLOCKS_X * BLOCKS_Y;
    std::vector<BlockResult> results(total_blocks);

    // Compress blocks in parallel; each thread grabs the next unprocessed block
    int nthreads = (int)std::thread::hardware_concurrency();
    if (nthreads < 1) nthreads = 1;
    std::atomic<int> next{0};

    auto worker = [&]() {
        int idx;
        while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < total_blocks) {
            int by = idx / BLOCKS_X;
            int bx = idx % BLOCKS_X;
            results[idx] = compress_block(image, bx, by);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(nthreads);
    for (int t = 0; t < nthreads; t++) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    // Write blocks in order (single-threaded, sequential)
    int mode_count[4] = {};
    for (int idx = 0; idx < total_blocks; idx++) {
        const auto& r = results[idx];
        mode_count[r.mode]++;
        write_u8(fout, r.mode);
        fwrite(r.payload.data(), 1, r.payload.size(), fout);
    }

    fclose(fout);

    // Relatório
    long in_bytes  = WIDTH * HEIGHT * 2;
    FILE* ftmp = fopen(argv[2], "rb");
    fseek(ftmp, 0, SEEK_END);
    long out_bytes = ftell(ftmp);
    fclose(ftmp);

    fprintf(stderr, "Modos: raw=%d  left=%d  MED=%d  smooth=%d\n",
            mode_count[0], mode_count[1], mode_count[2], mode_count[3]);
    fprintf(stderr, "Comprimido: %ld → %ld bytes  (%.4f bits/byte)\n",
            in_bytes, out_bytes, out_bytes * 8.0 / in_bytes);
    return 0;
}
