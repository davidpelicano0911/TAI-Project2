// HAIS — Hybrid Astronomical Image Compressor
// Usage: ./compress <input> <output.hais> [width height]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>

#include "coder/FseTable.hpp"
#include "coder/FseEncoder.hpp"
#include "model/CompressorModel.hpp"

static constexpr uint8_t MAGIC[4]   = {'H', 'A', 'I', 'S'};
static constexpr int     BLOCK_SIZE = 500;

// ---------------------------------------------------------------------------
// File write helpers
// ---------------------------------------------------------------------------

static void write_u8   (FILE* f, uint8_t  v) { fwrite(&v, 1, 1, f); }
static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 4, f);
}

// ---------------------------------------------------------------------------
// Per-block compression result
// ---------------------------------------------------------------------------

struct BlockResult {
    uint8_t              mode;
    std::vector<uint8_t> payload;
};

static BlockResult compress_block(const std::vector<uint16_t>& image,
                                  int W, int H, int bs,
                                  int bx, int by, uint16_t gmean) {
    int bw   = std::min(bs, W - bx * bs);
    int bh   = std::min(bs, H - by * bs);
    int npix = bw * bh;

    // Extract block pixels.
    std::vector<uint16_t> blk(npix);
    for (int y = 0; y < bh; y++) {
        int gy = by * bs + y;
        for (int x = 0; x < bw; x++)
            blk[y * bw + x] = image[gy * W + bx * bs + x];
    }

    // Compute residuals for all 5 modes.
    std::vector<uint16_t> s0, s1, s2, s3, s4;
    float ls_w[3], ls4_w[4];

    make_syms(blk, 0, bw, bh, s0, gmean);   // raw
    make_syms(blk, 1, bw, bh, s1, gmean);   // avg
    ls_compute (blk, bw, bh, ls_w,  s2);    // LS(W,N,NW)
    make_syms(blk, 3, bw, bh, s3, gmean);   // global_mean
    ls4_compute(blk, bw, bh, ls4_w, s4);    // LS+bias(W,N,NW,1)

    // Pick mode with lowest entropy cost (overhead for LS modes counted in bits).
    const double ls_overhead  = (12.0 * 8.0) / npix;
    const double ls4_overhead = (16.0 * 8.0) / npix;
    const std::vector<uint16_t>* sp[5] = {&s0, &s1, &s2, &s3, &s4};
    double costs[5] = {
        byte_cost(s0),
        byte_cost(s1),
        byte_cost(s2) + ls_overhead,
        byte_cost(s3),
        byte_cost(s4) + ls4_overhead,
    };
    int best = 0;
    for (int m = 1; m < 5; m++)
        if (costs[m] < costs[best]) best = m;

    const auto& sbest = *sp[best];
    std::vector<uint8_t> hi(npix), lo(npix);
    for (int i = 0; i < npix; i++) { hi[i] = sbest[i] >> 8; lo[i] = sbest[i] & 0xFF; }

    // Build payload: [weights if LS] [hi stream] [lo stream]
    BlockResult r;
    r.mode = (uint8_t)best;
    if (best == 2) {
        uint32_t bits;
        for (int i = 0; i < 3; i++) { memcpy(&bits, &ls_w[i], 4); bput32(r.payload, bits); }
    } else if (best == 4) {
        uint32_t bits;
        for (int i = 0; i < 4; i++) { memcpy(&bits, &ls4_w[i], 4); bput32(r.payload, bits); }
    }
    encode_stream(r.payload, hi);
    encode_stream(r.payload, lo);
    return r;
}

// ---------------------------------------------------------------------------
// Infer image dimensions from file size (npix = fsize/2 pixels).
// Tries square root first, then common widths, then 1×npix fallback.
// ---------------------------------------------------------------------------

static bool infer_dims(long fsize, int& W, int& H) {
    if (fsize <= 0 || fsize % 2 != 0) return false;
    long npix = fsize / 2;
    long sq   = (long)std::sqrt((double)npix);
    if (sq * sq == npix) { W = H = (int)sq; return true; }
    for (int w : {1500, 2048, 4096, 3000, 2000, 1920, 1024, 512, 256}) {
        if (npix % w == 0) { W = w; H = (int)(npix / w); return true; }
    }
    W = (int)npix; H = 1;
    return true;
}

// ---------------------------------------------------------------------------
// Choose block size: largest BS ≤ 512 that divides both W and H exactly,
// yielding at least 4 blocks. Falls back to 256 (with partial blocks).
// ---------------------------------------------------------------------------

static int choose_block_size(int W, int H) {
    for (int bs = 512; bs >= 32; bs--) {
        if (W % bs == 0 && H % bs == 0) {
            if ((W / bs) * (H / bs) >= 4) return bs;
        }
    }
    return 256;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 3 || argc == 4 || argc > 5) {
        fprintf(stderr, "Usage: %s <input> <output.hais> [width height]\n", argv[0]);
        return 1;
    }

    // Read input.
    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    fseek(fin, 0, SEEK_END); long fsize = ftell(fin); rewind(fin);

    int W, H;
    if (argc == 5) {
        W = std::atoi(argv[3]); H = std::atoi(argv[4]);
        if (W <= 0 || H <= 0) { fprintf(stderr, "Invalid dimensions\n"); fclose(fin); return 1; }
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

    std::vector<uint16_t> image(W * H);
    for (int i = 0; i < W * H; i++) {
        uint8_t b[2];
        if (fread(b, 1, 2, fin) != 2) { fprintf(stderr, "Short read\n"); return 1; }
        image[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    // Compute global mean for mode 3.
    uint64_t sum = 0;
    for (auto px : image) sum += px;
    uint16_t gmean = (uint16_t)(sum / (W * H));

    int bs       = choose_block_size(W, H);
    int blocks_x = (W + bs - 1) / bs;
    int blocks_y = (H + bs - 1) / bs;
    int total    = blocks_x * blocks_y;

    // Compress blocks in parallel.
    std::vector<BlockResult> results(total);
    std::atomic<int> next{0};
    int nthreads = std::max(1, (int)std::thread::hardware_concurrency());
    auto worker = [&]() {
        int idx;
        while ((idx = next.fetch_add(1, std::memory_order_relaxed)) < total)
            results[idx] = compress_block(image, W, H, bs,
                                          idx % blocks_x, idx / blocks_x, gmean);
    };
    std::vector<std::thread> pool(nthreads);
    for (auto& t : pool) t = std::thread(worker);
    for (auto& t : pool) t.join();

    // Write output.
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u32le(fout, (uint32_t)W);
    write_u32le(fout, (uint32_t)H);
    write_u32le(fout, (uint32_t)bs);
    write_u16le(fout, gmean);

    int mode_count[5] = {};
    for (int idx = 0; idx < total; idx++) {
        const auto& r = results[idx];
        mode_count[r.mode]++;
        write_u8(fout, r.mode);
        fwrite(r.payload.data(), 1, r.payload.size(), fout);
    }
    fclose(fout);

    long in_bytes  = (long)W * H * 2;
    FILE* ft = fopen(argv[2], "rb"); fseek(ft, 0, SEEK_END); long out_bytes = ftell(ft); fclose(ft);
    fprintf(stderr, "Dimensions: %dx%d  blocks: %dx%d (%d total)\n",
            W, H, blocks_x, blocks_y, total);
    fprintf(stderr, "Modes: raw=%d  avg=%d  ls=%d  mean=%d  ls4=%d\n",
            mode_count[0], mode_count[1], mode_count[2], mode_count[3], mode_count[4]);
    fprintf(stderr, "Compressed: %ld -> %ld bytes  (%.4f bits/byte)\n",
            in_bytes, out_bytes, out_bytes * 8.0 / in_bytes);
    return 0;
}
