/*
 * TAI Project 2 — Astronomical Image Decompressor (ox-astro-balanced)
 *
 * Usage:  decompress_astro_balanced <input_file> <output_file>
 *         decompress_astro_balanced          (stdin → stdout)
 */

#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

#include "coder/RangeCoder.hpp"
#include "coder/FrequencyTable.hpp"
#include "coder/FenwickFrequencyTable.hpp"
#include "model/ImagePredictor.hpp"

static constexpr uint8_t MAGIC[4] = {'T','A','2','B'};

static uint32_t read_u32le(std::istream& in) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i)
        v |= static_cast<uint32_t>(static_cast<uint8_t>(in.get())) << (8*i);
    return v;
}

static uint16_t read_u16le(std::istream& in) {
    uint16_t lo = static_cast<uint8_t>(in.get());
    uint16_t hi = static_cast<uint8_t>(in.get());
    return static_cast<uint16_t>(lo | (hi << 8));
}

static FenwickFrequencyTable make_fenwick_from_hist(const std::array<uint8_t,256>& hist) {
    FenwickFrequencyTable t(256);
    for (int i = 0; i < 256; ++i) t.set(i, hist[i]);
    return t;
}

int main(int argc, char* argv[]) {
    std::istream* in_ptr  = &std::cin;
    std::ostream* out_ptr = &std::cout;
    std::ifstream fin;
    std::ofstream fout;

    if (argc >= 3) {
        fin.open(argv[1], std::ios::binary);
        if (!fin)  { std::cerr << "Cannot open input: "  << argv[1] << '\n'; return 1; }
        fout.open(argv[2], std::ios::binary);
        if (!fout) { std::cerr << "Cannot open output: " << argv[2] << '\n'; return 1; }
        in_ptr  = &fin;
        out_ptr = &fout;
    } else if (argc == 1) {
        std::cin.sync_with_stdio(false);
    } else {
        std::cerr << "Usage: decompress_astro_balanced <input> <output>\n"; return 1;
    }

    uint8_t magic[4];
    in_ptr->read(reinterpret_cast<char*>(magic), 4);
    for (int i = 0; i < 4; ++i)
        if (magic[i] != MAGIC[i]) { std::cerr << "Bad magic\n"; return 1; }

    const uint32_t width  = read_u32le(*in_ptr);
    const uint32_t height = read_u32le(*in_ptr);
    const uint64_t npix   = static_cast<uint64_t>(width) * height;
    const uint16_t global_mean = read_u16le(*in_ptr);

    // Read both sets of warm-model histograms.
    std::array<uint8_t,256> hi_hist_gap{},    lo_hi0hist_gap{},    lo_hiphist_gap{};
    std::array<uint8_t,256> hi_hist_mean{},   lo_hi0hist_mean{},   lo_hiphist_mean{};
    in_ptr->read(reinterpret_cast<char*>(hi_hist_gap.data()),    256);
    in_ptr->read(reinterpret_cast<char*>(lo_hi0hist_gap.data()), 256);
    in_ptr->read(reinterpret_cast<char*>(lo_hiphist_gap.data()), 256);
    in_ptr->read(reinterpret_cast<char*>(hi_hist_mean.data()),    256);
    in_ptr->read(reinterpret_cast<char*>(lo_hi0hist_mean.data()), 256);
    in_ptr->read(reinterpret_cast<char*>(lo_hiphist_mean.data()), 256);

    const unsigned ntiles = static_cast<uint8_t>(in_ptr->get());
    std::vector<uint8_t> tile_pred_modes(ntiles);
    in_ptr->read(reinterpret_cast<char*>(tile_pred_modes.data()), ntiles);

    std::vector<uint32_t> tile_sizes(ntiles);
    for (unsigned t = 0; t < ntiles; ++t) tile_sizes[t] = read_u32le(*in_ptr);

    std::vector<std::string> tile_data(ntiles);
    for (unsigned t = 0; t < ntiles; ++t) {
        tile_data[t].resize(tile_sizes[t]);
        in_ptr->read(tile_data[t].data(), tile_sizes[t]);
    }

    auto quantize = [](uint8_t h) -> int {
        if (h == 0) return 0; if (h == 1) return 1; if (h <= 7) return 2; return 3;
    };

    auto quantize_prev_lo = [](uint8_t lo) -> int {
        if (lo == 0) return 0; if (lo < 4) return 1; if (lo < 8) return 2;
        if (lo < 16) return 3; if (lo < 32) return 4; if (lo < 64) return 5;
        if (lo < 128) return 6; return 7;
    };

    std::vector<uint16_t> pixels(npix);
    std::vector<std::thread> threads(ntiles);

    for (unsigned t = 0; t < ntiles; ++t) {
        threads[t] = std::thread([&, t]() {
            const uint32_t row_start  = static_cast<uint32_t>(height * t / ntiles);
            const uint32_t row_end    = static_cast<uint32_t>(height * (t+1) / ntiles);
            const bool     use_mean   = (tile_pred_modes[t] == 1);

            const auto& wh    = use_mean ? hi_hist_mean    : hi_hist_gap;
            const auto& wlh0  = use_mean ? lo_hi0hist_mean : lo_hi0hist_gap;
            const auto& wlhip = use_mean ? lo_hiphist_mean : lo_hiphist_gap;

            std::istringstream iss(tile_data[t], std::ios::binary);
            std::vector<FenwickFrequencyTable> hi_models(16,  make_fenwick_from_hist(wh));
            std::vector<FenwickFrequencyTable> lo_hi0_models(32,  make_fenwick_from_hist(wlh0));
            std::vector<FenwickFrequencyTable> lo_hip_models(255, make_fenwick_from_hist(wlhip));

            std::vector<int> C(16,0), B(16,0), Nc(16,0);
            std::vector<uint8_t> hi_N(width, 0);
            uint8_t prev_lo_hi0 = 0;

            RangeDecoder dec(iss);

            for (uint32_t row = row_start; row < row_end; ++row) {
                const uint32_t local_row = row - row_start;
                uint8_t hi_W = 0;
                for (uint32_t col = 0; col < width; ++col) {
                    uint16_t W  = (col > 0)                    ? pixels[row*width+col-1]         : 0u;
                    uint16_t WW = (col > 1)                    ? pixels[row*width+col-2]         : W;
                    uint16_t N  = (local_row > 0)              ? pixels[(row-1)*width+col]       : W;
                    uint16_t NN = (local_row > 1)              ? pixels[(row-2)*width+col]       : N;
                    uint16_t NW = (local_row>0 && col>0)       ? pixels[(row-1)*width+col-1]     : W;
                    uint16_t NE = (local_row>0 && col<width-1) ? pixels[(row-1)*width+col+1]     : N;

                    uint8_t ctx_hi = static_cast<uint8_t>(quantize(hi_W)*4 + quantize(hi_N[col]));

                    uint8_t hi = static_cast<uint8_t>(dec.read(hi_models[ctx_hi]));
                    hi_models[ctx_hi].increment(hi);

                    uint8_t lo;
                    if (hi == 0) {
                        int gc = std::min(quantize(hi_W) + quantize(hi_N[col]), 3);
                        int lo_idx = gc * 8 + quantize_prev_lo(prev_lo_hi0);
                        lo = static_cast<uint8_t>(dec.read(lo_hi0_models[lo_idx]));
                        lo_hi0_models[lo_idx].increment(lo);
                        prev_lo_hi0 = lo;
                    } else {
                        lo = static_cast<uint8_t>(dec.read(lo_hip_models[hi-1]));
                        lo_hip_models[hi-1].increment(lo);
                    }

                    uint16_t u = zigzag_decode((static_cast<uint16_t>(hi) << 8) | lo);

                    uint16_t pred;
                    if (use_mean) {
                        pred = global_mean;
                    } else {
                        pred = (local_row==0 && col==0)
                            ? 0u
                            : robust_gap_predict(W,WW,N,NN,NW,NE,hi_W,hi_N[col]);
                    }
                    uint16_t pred_adj = static_cast<uint16_t>((int)pred + C[ctx_hi]);
                    pixels[row * width + col] = static_cast<uint16_t>(pred_adj + u);

                    int r_s = (u <= 32767u) ? (int)u : (int)u - 65536;
                    B[ctx_hi] += r_s; Nc[ctx_hi]++;
                    if      (B[ctx_hi] >  Nc[ctx_hi]) { C[ctx_hi]++; B[ctx_hi] -= Nc[ctx_hi]; if (B[ctx_hi]>Nc[ctx_hi])  B[ctx_hi]=Nc[ctx_hi]; }
                    else if (B[ctx_hi] < -Nc[ctx_hi]) { C[ctx_hi]--; B[ctx_hi] += Nc[ctx_hi]; if (B[ctx_hi]<-Nc[ctx_hi]) B[ctx_hi]=-Nc[ctx_hi]; }
                    if (Nc[ctx_hi] == 64) { Nc[ctx_hi] >>= 1; B[ctx_hi] >>= 1; }

                    hi_W = hi; hi_N[col] = hi;
                }
            }
        });
    }

    for (auto& th : threads) th.join();

    std::vector<uint8_t> out_bytes(npix * 2);
    for (uint64_t i = 0; i < npix; ++i) {
        out_bytes[2*i]   = static_cast<uint8_t>(pixels[i] >> 8);
        out_bytes[2*i+1] = static_cast<uint8_t>(pixels[i] & 0xFF);
    }
    out_ptr->write(reinterpret_cast<const char*>(out_bytes.data()),
                   static_cast<std::streamsize>(out_bytes.size()));
    return 0;
}
