// HAIS — Hybrid Astronomical Image Decompressor
// Usage: ./decompress <input.hais> <output>

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>

#include "coder/FseTable.hpp"
#include "coder/FseDecoder.hpp"
#include "model/Predictors.hpp"

static constexpr uint8_t MAGIC[4] = {'H', 'A', 'I', 'S'};

// ---------------------------------------------------------------------------
// Per-block decompression task (read-only pointer into the file buffer)
// ---------------------------------------------------------------------------

struct BlockTask {
    const uint8_t* ptr;
    int     mode, bw, bh;
    uint16_t gmean;
};

static void decompress_block(const BlockTask& task, std::vector<uint16_t>& blk) {
    const uint8_t* ptr  = task.ptr;
    int             npix = task.bw * task.bh;

    // Read per-block LS weights if present.
    float ls_w[3] = {}, ls4_w[4] = {};
    if (task.mode == 2) {
        for (int i = 0; i < 3; i++) {
            uint32_t bits = pget_u32le(ptr);
            memcpy(&ls_w[i], &bits, 4);
        }
    } else if (task.mode == 4) {
        for (int i = 0; i < 4; i++) {
            uint32_t bits = pget_u32le(ptr);
            memcpy(&ls4_w[i], &bits, 4);
        }
    }

    // Decode hi8 and lo8 streams.
    std::vector<uint8_t> hi8, lo8;
    decode_stream(ptr, hi8, npix);
    decode_stream(ptr, lo8, npix);

    // Reconstruct pixels.
    blk.resize(npix);
    for (int y = 0; y < task.bh; y++) {
        for (int x = 0; x < task.bw; x++) {
            int      idx = y * task.bw + x;
            uint16_t sym = ((uint16_t)hi8[idx] << 8) | lo8[idx];
            uint16_t pixel;

            switch (task.mode) {
            case 0:  // raw
                pixel = sym;
                break;
            case 1:  // avg
                pixel = (uint16_t)((int)avg_pred(blk, x, y, task.bw) + zagzig(sym));
                break;
            case 3:  // global_mean
                pixel = (uint16_t)((int)task.gmean + zagzig(sym));
                break;
            case 2: {  // LS(W, N, NW)
                uint16_t pred;
                if (y == 0 && x == 0)
                    pred = 0;
                else if (y == 0)
                    pred = blk[x - 1];
                else if (x == 0)
                    pred = blk[(y-1)*task.bw+x];
                else {
                    float p = ls_w[0]*blk[y*task.bw+(x-1)]
                            + ls_w[1]*blk[(y-1)*task.bw+x]
                            + ls_w[2]*blk[(y-1)*task.bw+(x-1)];
                    pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
                }
                pixel = (uint16_t)((int)pred + zagzig(sym));
                break;
            }
            case 4: {  // LS+bias(W, N, NW, 1)
                uint16_t pred;
                if (y == 0 && x == 0)
                    pred = (uint16_t)std::max(0.0f, std::min(65535.0f, ls4_w[3] + 0.5f));
                else if (y == 0)
                    pred = (uint16_t)std::max(0.0f, std::min(65535.0f,
                               ls4_w[0]*blk[x-1] + ls4_w[3] + 0.5f));
                else if (x == 0)
                    pred = (uint16_t)std::max(0.0f, std::min(65535.0f,
                               ls4_w[1]*blk[(y-1)*task.bw+x] + ls4_w[3] + 0.5f));
                else {
                    float p = ls4_w[0]*blk[y*task.bw+(x-1)]
                            + ls4_w[1]*blk[(y-1)*task.bw+x]
                            + ls4_w[2]*blk[(y-1)*task.bw+(x-1)]
                            + ls4_w[3];
                    pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
                }
                pixel = (uint16_t)((int)pred + zagzig(sym));
                break;
            }
            default:
                pixel = sym;
                break;
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

    // Read entire file into memory.
    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    fseek(fin, 0, SEEK_END); long fsize = ftell(fin); rewind(fin);
    std::vector<uint8_t> buf(fsize);
    [[maybe_unused]] auto rb = fread(buf.data(), 1, fsize, fin);
    fclose(fin);

    // Parse header.
    const uint8_t* p = buf.data();
    if (memcmp(p, MAGIC, 4)) { fprintf(stderr, "Bad magic\n"); return 1; }
    p += 4;
    uint32_t width  = pget_u32le(p);
    uint32_t height = pget_u32le(p);
    uint32_t bs     = pget_u32le(p);
    uint16_t gmean  = pget_u16le(p);

    int blocks_x = ((int)width  + (int)bs - 1) / (int)bs;
    int blocks_y = ((int)height + (int)bs - 1) / (int)bs;
    int total    = blocks_x * blocks_y;

    // Sequential pass: record block pointers and skip over payloads.
    std::vector<BlockTask>        tasks(total);
    std::vector<std::pair<int,int>> block_pos(total);

    for (int by = 0; by < blocks_y; by++) {
        for (int bx = 0; bx < blocks_x; bx++) {
            int idx = by * blocks_x + bx;
            int bw  = std::min((int)bs, (int)width  - bx * (int)bs);
            int bh  = std::min((int)bs, (int)height - by * (int)bs);
            int mode = (int)pget_u8(p);

            tasks[idx]     = {p, mode, bw, bh, gmean};
            block_pos[idx] = {bx, by};

            // Skip weights.
            if      (mode == 2) p += 12;
            else if (mode == 4) p += 16;

            // Skip 2 FSE streams.
            for (int s = 0; s < 2; s++) {
                uint8_t flag = pget_u8(p);
                if      (flag == 0xFF) { /* uniform — no table bytes */ }
                else if (flag == 0)    { p += 256 * 4; }
                else                   { p += (int)flag * 5; }
                uint64_t nbytes = pget_u64le(p);
                p += nbytes;
            }
        }
    }

    // Parallel decompression.
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

    // Assemble full image.
    std::vector<uint16_t> image(width * height);
    for (int idx = 0; idx < total; idx++) {
        auto [bx, by] = block_pos[idx];
        int bw = tasks[idx].bw, bh = tasks[idx].bh;
        const auto& blk = blocks[idx];
        for (int y = 0; y < bh; y++) {
            int gy = by * bs + y;
            for (int x = 0; x < bw; x++)
                image[gy * width + bx * bs + x] = blk[y * bw + x];
        }
    }

    // Write big-endian u16 output.
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }
    for (uint16_t px : image) {
        uint8_t b[2] = {(uint8_t)(px >> 8), (uint8_t)(px & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
