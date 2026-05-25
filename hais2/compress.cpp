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

static constexpr uint8_t MAGIC[4] = {'H', 'A', 'I', '2'};

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

    std::vector<uint16_t> s0, s1, s3, s4, s6, s9;
    float ls4_w[4], ls5_w[5], ls6_w[6];

    make_syms(blk, 0, bw, bh, s0, gmean);    // raw
    make_syms(blk, 1, bw, bh, s1, gmean);    // avg
    make_syms(blk, 3, bw, bh, s3, gmean);    // global_mean
    ls4_compute(blk, bw, bh, ls4_w, s4);     // LS+bias(W,N,NW,1)
    ls5_compute(blk, bw, bh, ls5_w, s6);     // LS(W,N,NW,NE,NN)
    ls6_compute(blk, bw, bh, ls6_w, s9);     // LS+bias(W,N,NW,WW,NN,1)

    // Use ctx_cost for all candidates — better matches actual coded size.
    const double ls4_overhead = (16.0 * 8.0) / npix;
    const double ls5_overhead = (20.0 * 8.0) / npix;
    const double ls6_overhead = (24.0 * 8.0) / npix;
    struct Candidate { int mode; double cost; const std::vector<uint16_t>* syms; };
    Candidate cands[6] = {
        {0, ctx_cost(s0),                  &s0},
        {1, ctx_cost(s1),                  &s1},
        {3, ctx_cost(s3),                  &s3},
        {4, ctx_cost(s4) + ls4_overhead,   &s4},
        {6, ctx_cost(s6) + ls5_overhead,   &s6},
        {9, ctx_cost(s9) + ls6_overhead,   &s9},
    };
    int best = 0;
    for (int m = 1; m < 6; m++)
        if (cands[m].cost < cands[best].cost) best = m;

    const auto& sbest = *cands[best].syms;
    int   mode_id     =  cands[best].mode;
    std::vector<uint8_t> hi(npix), lo(npix);
    for (int i = 0; i < npix; i++) { hi[i] = sbest[i] >> 8; lo[i] = sbest[i] & 0xFF; }

    // Build weights header.
    std::vector<uint8_t> whdr;
    if (mode_id == 4) {
        uint32_t bits;
        for (int i = 0; i < 4; i++) { memcpy(&bits, &ls4_w[i], 4); bput32(whdr, bits); }
    } else if (mode_id == 6) {
        uint32_t bits;
        for (int i = 0; i < 5; i++) { memcpy(&bits, &ls5_w[i], 4); bput32(whdr, bits); }
    } else if (mode_id == 9) {
        uint32_t bits;
        for (int i = 0; i < 6; i++) { memcpy(&bits, &ls6_w[i], 4); bput32(whdr, bits); }
    }

    // Encode hi once — shared between plain and ctx payloads.
    std::vector<uint8_t> hi_stream;
    encode_stream(hi_stream, hi);

    // Build plain payload.
    std::vector<uint8_t> plain;
    plain.insert(plain.end(), whdr.begin(), whdr.end());
    plain.insert(plain.end(), hi_stream.begin(), hi_stream.end());
    encode_stream(plain, lo);

    // Build context payload (hi==0 and hi>0 lo-byte streams separately).
    std::vector<uint8_t> lo_zero, lo_nonzero;
    lo_zero.reserve(npix); lo_nonzero.reserve(npix);
    for (int i = 0; i < npix; i++) {
        if (hi[i] == 0) lo_zero.push_back(lo[i]);
        else            lo_nonzero.push_back(lo[i]);
    }
    int nzero = (int)lo_zero.size();
    BlockResult r;
    if (nzero >= npix / 50 && nzero <= npix - npix / 50) {
        std::vector<uint8_t> ctx;
        ctx.insert(ctx.end(), whdr.begin(), whdr.end());
        ctx.insert(ctx.end(), hi_stream.begin(), hi_stream.end());
        encode_stream(ctx, lo_zero);
        encode_stream(ctx, lo_nonzero);
        if (ctx.size() < plain.size()) {
            r.mode    = (uint8_t)(mode_id | 0x80);
            r.payload = std::move(ctx);
            return r;
        }
    }
    r.mode    = (uint8_t)mode_id;
    r.payload = std::move(plain);
    return r;
}

// ---------------------------------------------------------------------------
// Infer image dimensions from file size (npix = fsize/2 pixels).
// Tries square root first, then common widths, then 1×npix fallback.
// ---------------------------------------------------------------------------

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
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <n_rows> <n_cols> <input> <output.hais>\n", argv[0]);
        return 1;
    }

    int H = std::atoi(argv[1]);
    int W = std::atoi(argv[2]);
    if (W <= 0 || H <= 0) { fprintf(stderr, "Invalid dimensions\n"); return 1; }

    FILE* fin = fopen(argv[3], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[3]); return 1; }
    fseek(fin, 0, SEEK_END); long fsize = ftell(fin); rewind(fin);

    if (fsize != (long)W * H * 2) {
        fprintf(stderr, "File size mismatch: got %ld, expected %ld\n", fsize, (long)W * H * 2);
        fclose(fin); return 1;
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
    FILE* fout = fopen(argv[4], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[4]); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u32le(fout, (uint32_t)W);
    write_u32le(fout, (uint32_t)H);
    write_u32le(fout, (uint32_t)bs);
    write_u16le(fout, gmean);

    int mode_count[10] = {}, ctx_count = 0;
    for (int idx = 0; idx < total; idx++) {
        const auto& r = results[idx];
        mode_count[r.mode & 0x7F]++;
        if (r.mode & 0x80) ctx_count++;
        write_u8(fout, r.mode);
        fwrite(r.payload.data(), 1, r.payload.size(), fout);
    }
    fclose(fout);

    long in_bytes  = (long)W * H * 2;
    FILE* ft = fopen(argv[4], "rb"); fseek(ft, 0, SEEK_END); long out_bytes = ftell(ft); fclose(ft);
    fprintf(stderr, "Dimensions: %dx%d  blocks: %dx%d (%d total)\n",
            W, H, blocks_x, blocks_y, total);
    fprintf(stderr, "Modes: raw=%d avg=%d mean=%d ls4=%d ls5=%d ls6=%d  ctx=%d/%d\n",
            mode_count[0], mode_count[1], mode_count[3],
            mode_count[4], mode_count[6], mode_count[9],
            ctx_count, total);
    fprintf(stderr, "Compressed: %ld -> %ld bytes  (%.4f bits/byte)\n",
            in_bytes, out_bytes, out_bytes * 8.0 / in_bytes);
    return 0;
}
