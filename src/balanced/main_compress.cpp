/*
 * TAI Project 2 — Astronomical Image Compressor (ox-astro-balanced)
 *
 * Pipeline:
 *   1. Read raw bytes; interpret as 2D big-endian uint16_t array.
 *   2. Global entropy scan: accumulate 6 × 256-bucket histograms for both
 *        GAP and mean predictors. Scale both sets for warm-model storage.
 *   3. Tile-parallel encoding (ntiles = min(8, hardware_concurrency)):
 *        Each tile runs its own entropy scan over its row range and
 *        independently selects GAP or mean predictor. Models are warm-started
 *        from the global histogram set matching the tile's selected predictor.
 *   4. Wrapping residual → modular zigzag → hi/lo byte split.
 *   5. Per-context bias correction (16 contexts, LOCO-I style).
 *   6. Encode hi with spatial context: quantize(hi_W)*4 + quantize(hi_N).
 *   7. Encode lo:
 *        hi == 0 : model[(grad_class)*8 + quantize_prev_lo]  (32 models)
 *        hi  > 0 : model[hi - 1]                             (255 models)
 *
 * Compressed file format
 * ──────────────────────
 *  Bytes  0–  3 : magic "TA2B"
 *  Bytes  4–  7 : width           uint32_t LE
 *  Bytes  8– 11 : height          uint32_t LE
 *  Bytes 12– 13 : global_mean     uint16_t LE
 *  Byte  14   : ntiles              uint8_t
 *  Bytes 15.. : pred_mode[ntiles]  uint8_t each (0=GAP, 1=mean, per tile)
 *  Then       : tile_hists[ntiles × 3 × 256] uint8_t (hi + lo_hi0 + lo_hip per tile)
 *  Then       : tile_sizes[ntiles] uint32_t LE each
 *  Then       : tile bitstreams concatenated
 *
 * Usage:  compress_astro_balanced <input_file> <output_file>
 *         compress_astro_balanced          (stdin → stdout)
 */

#include <algorithm>
#include <array>
#include <cmath>
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

static void write_u32le(std::ostream& out, uint32_t v) {
    out.put(static_cast<char>(v         & 0xFF));
    out.put(static_cast<char>((v >>  8) & 0xFF));
    out.put(static_cast<char>((v >> 16) & 0xFF));
    out.put(static_cast<char>((v >> 24) & 0xFF));
}

static void write_u16le(std::ostream& out, uint16_t v) {
    out.put(static_cast<char>(v       & 0xFF));
    out.put(static_cast<char>((v >> 8) & 0xFF));
}

static std::array<uint8_t,256> scale_hist(const std::array<uint64_t,256>& cnt) {
    std::array<uint8_t,256> out{};
    uint64_t max_cnt = *std::max_element(cnt.begin(), cnt.end());
    for (int i = 0; i < 256; ++i) {
        uint64_t s = (max_cnt > 0) ? (cnt[i] * 127u / max_cnt) : 0u;
        out[i] = static_cast<uint8_t>(s + 1u);
    }
    return out;
}

static FenwickFrequencyTable make_fenwick_from_hist(const std::array<uint8_t,256>& hist) {
    FenwickFrequencyTable t(256);
    for (int i = 0; i < 256; ++i)
        t.set(i, hist[i]);
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
        std::cerr << "Usage: compress_astro_balanced <input> <output>\n"; return 1;
    }

    std::vector<uint8_t> raw(std::istreambuf_iterator<char>(*in_ptr), {});
    if (raw.size() % 2 != 0) { std::cerr << "Input size not even\n"; return 1; }

    const uint64_t npix = raw.size() / 2;
    uint32_t width = 0, height = 0;
    uint32_t sq = static_cast<uint32_t>(std::sqrt(static_cast<double>(npix)));
    if ((uint64_t)sq * sq == npix) {
        width = height = sq;
    } else {
        width = static_cast<uint32_t>(npix); height = 1;
        for (uint32_t w : {1500u, 2048u, 1024u, 512u, 256u})
            if (npix % w == 0) { width = w; height = static_cast<uint32_t>(npix / w); break; }
    }

    std::vector<uint16_t> pixels(npix);
    for (uint64_t i = 0; i < npix; ++i)
        pixels[i] = (static_cast<uint16_t>(raw[2*i]) << 8) | raw[2*i+1];

    // ── Global mean ────────────────────────────────────────────────────────
    uint64_t px_sum = 0;
    for (uint64_t i = 0; i < npix; ++i) px_sum += pixels[i];
    const uint16_t global_mean = static_cast<uint16_t>(px_sum / npix);

    // ── Shared helpers ─────────────────────────────────────────────────────
    auto entropy = [](const std::array<uint64_t,256>& cnt) -> double {
        uint64_t total = 0; for (auto c : cnt) total += c;
        if (total == 0) return 0.0;
        double h = 0.0; const double inv = 1.0 / static_cast<double>(total);
        for (auto c : cnt) if (c > 0) { double p = c * inv; h -= p * std::log2(p); }
        return h;
    };

    auto quantize = [](uint8_t h) -> int {
        if (h == 0) return 0; if (h == 1) return 1; if (h <= 7) return 2; return 3;
    };

    auto quantize_prev_lo = [](uint8_t lo) -> int {
        if (lo == 0) return 0; if (lo < 4) return 1; if (lo < 8) return 2;
        if (lo < 16) return 3; if (lo < 32) return 4; if (lo < 64) return 5;
        if (lo < 128) return 6; return 7;
    };

    // ── Tile-parallel Phase 1: scan (parallel) ────────────────────────────
    const unsigned ntiles = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
    std::vector<std::string>  tile_bufs(ntiles);
    std::vector<uint8_t>      tile_pred_modes(ntiles);
    std::vector<std::thread>  threads(ntiles);

    std::vector<std::array<uint64_t,256>> t_hi_gap(ntiles),     t_hi_mean(ntiles);
    std::vector<std::array<uint64_t,256>> t_lo_hi0_gap(ntiles), t_lo_hi0_mean(ntiles);
    std::vector<std::array<uint64_t,256>> t_lo_hip_gap(ntiles), t_lo_hip_mean(ntiles);

    for (unsigned t = 0; t < ntiles; ++t) {
        threads[t] = std::thread([&, t]() {
            const uint32_t row_start = static_cast<uint32_t>(height * t / ntiles);
            const uint32_t row_end   = static_cast<uint32_t>(height * (t + 1) / ntiles);
            const uint64_t tile_npix = static_cast<uint64_t>(row_end - row_start) * width;

            std::array<uint64_t,256> ti_hi_gap{},     ti_hi_mean{};
            std::array<uint64_t,256> ti_lo_hi0_gap{}, ti_lo_hi0_mean{};
            std::array<uint64_t,256> ti_lo_hip_gap{}, ti_lo_hip_mean{};

            std::vector<uint8_t> scan_hi_N(width, 0);
            for (uint32_t row = row_start; row < row_end; ++row) {
                const uint32_t lr = row - row_start;
                uint8_t scan_hi_W = 0;
                for (uint32_t col = 0; col < width; ++col) {
                    uint16_t px = pixels[row * width + col];
                    uint16_t W  = (col > 0)               ? pixels[row*width+col-1]           : 0u;
                    uint16_t WW = (col > 1)               ? pixels[row*width+col-2]           : W;
                    uint16_t N  = (lr > 0)                ? pixels[(row-1)*width+col]         : W;
                    uint16_t NN = (lr > 1)                ? pixels[(row-2)*width+col]         : N;
                    uint16_t NW = (lr > 0 && col > 0)     ? pixels[(row-1)*width+col-1]       : W;
                    uint16_t NE = (lr > 0 && col<width-1) ? pixels[(row-1)*width+col+1]       : N;

                    uint16_t pg = (lr == 0 && col == 0) ? 0u
                                : robust_gap_predict(W,WW,N,NN,NW,NE,scan_hi_W,scan_hi_N[col]);
                    uint16_t zg = zigzag_encode(static_cast<uint16_t>((int)px-(int)pg));
                    uint8_t hig = zg>>8; uint8_t log_ = zg&0xFF;
                    ti_hi_gap[hig]++;
                    if (hig==0) ti_lo_hi0_gap[log_]++; else ti_lo_hip_gap[log_]++;

                    uint16_t zm = zigzag_encode(static_cast<uint16_t>((int)px-(int)global_mean));
                    uint8_t him = zm>>8; uint8_t lom = zm&0xFF;
                    ti_hi_mean[him]++;
                    if (him==0) ti_lo_hi0_mean[lom]++; else ti_lo_hip_mean[lom]++;

                    scan_hi_W = hig;
                    scan_hi_N[col] = hig;
                }
            }

            auto tile_est = [&](const std::array<uint64_t,256>& hc,
                                const std::array<uint64_t,256>& lc0,
                                const std::array<uint64_t,256>& lcp) -> double {
                uint64_t n0 = 0; for (auto c : lc0) n0 += c;
                double f0 = static_cast<double>(n0) / static_cast<double>(tile_npix);
                return entropy(hc) + f0*entropy(lc0) + (1.0-f0)*entropy(lcp);
            };

            const bool use_mean = tile_est(ti_hi_mean, ti_lo_hi0_mean, ti_lo_hip_mean) <
                                  tile_est(ti_hi_gap,  ti_lo_hi0_gap,  ti_lo_hip_gap);
            tile_pred_modes[t] = use_mean ? 1u : 0u;

            t_hi_gap[t]     = ti_hi_gap;     t_hi_mean[t]    = ti_hi_mean;
            t_lo_hi0_gap[t] = ti_lo_hi0_gap; t_lo_hi0_mean[t]= ti_lo_hi0_mean;
            t_lo_hip_gap[t] = ti_lo_hip_gap; t_lo_hip_mean[t]= ti_lo_hip_mean;
        });
    }
    for (auto& th : threads) th.join();

    // Merge per-tile histograms → global (full-image quality warm-start).
    std::array<uint64_t,256> g_hi_gap{},     g_hi_mean{};
    std::array<uint64_t,256> g_lo_hi0_gap{}, g_lo_hi0_mean{};
    std::array<uint64_t,256> g_lo_hip_gap{}, g_lo_hip_mean{};
    for (unsigned t = 0; t < ntiles; ++t)
        for (int i = 0; i < 256; ++i) {
            g_hi_gap[i]     += t_hi_gap[t][i];     g_hi_mean[i]    += t_hi_mean[t][i];
            g_lo_hi0_gap[i] += t_lo_hi0_gap[t][i]; g_lo_hi0_mean[i]+= t_lo_hi0_mean[t][i];
            g_lo_hip_gap[i] += t_lo_hip_gap[t][i]; g_lo_hip_mean[i]+= t_lo_hip_mean[t][i];
        }

    const auto hi_hist_gap     = scale_hist(g_hi_gap);
    const auto lo_hi0hist_gap  = scale_hist(g_lo_hi0_gap);
    const auto lo_hiphist_gap  = scale_hist(g_lo_hip_gap);
    const auto hi_hist_mean    = scale_hist(g_hi_mean);
    const auto lo_hi0hist_mean = scale_hist(g_lo_hi0_mean);
    const auto lo_hiphist_mean = scale_hist(g_lo_hip_mean);

    // ── Tile-parallel Phase 2: encode ──────────────────────────────────────
    for (unsigned t = 0; t < ntiles; ++t) {
        threads[t] = std::thread([&, t]() {
            const uint32_t row_start = static_cast<uint32_t>(height * t / ntiles);
            const uint32_t row_end   = static_cast<uint32_t>(height * (t + 1) / ntiles);
            const bool use_mean = (tile_pred_modes[t] == 1);

            const auto& wh    = use_mean ? hi_hist_mean    : hi_hist_gap;
            const auto& wlh0  = use_mean ? lo_hi0hist_mean : lo_hi0hist_gap;
            const auto& wlhip = use_mean ? lo_hiphist_mean : lo_hiphist_gap;

            std::ostringstream oss;
            std::vector<FenwickFrequencyTable> hi_models(16,  make_fenwick_from_hist(wh));
            std::vector<FenwickFrequencyTable> lo_hi0_models(32,  make_fenwick_from_hist(wlh0));
            std::vector<FenwickFrequencyTable> lo_hip_models(255, make_fenwick_from_hist(wlhip));

            std::vector<int> C(16,0), B(16,0), Nc(16,0);
            std::vector<uint8_t> hi_N(width, 0);
            uint8_t prev_lo_hi0 = 0;

            RangeEncoder enc(oss);

            for (uint32_t row = row_start; row < row_end; ++row) {
                const uint32_t local_row = row - row_start;
                uint8_t hi_W = 0;
                for (uint32_t col = 0; col < width; ++col) {
                    uint16_t px = pixels[row * width + col];
                    uint16_t W  = (col > 0)                    ? pixels[row*width+col-1]         : 0u;
                    uint16_t WW = (col > 1)                    ? pixels[row*width+col-2]         : W;
                    uint16_t N  = (local_row > 0)              ? pixels[(row-1)*width+col]       : W;
                    uint16_t NN = (local_row > 1)              ? pixels[(row-2)*width+col]       : N;
                    uint16_t NW = (local_row>0 && col>0)       ? pixels[(row-1)*width+col-1]     : W;
                    uint16_t NE = (local_row>0 && col<width-1) ? pixels[(row-1)*width+col+1]     : N;

                    uint8_t ctx_hi = static_cast<uint8_t>(quantize(hi_W)*4 + quantize(hi_N[col]));

                    uint16_t pred;
                    if (use_mean) {
                        pred = global_mean;
                    } else {
                        pred = (local_row==0 && col==0)
                            ? 0u
                            : robust_gap_predict(W,WW,N,NN,NW,NE,hi_W,hi_N[col]);
                    }
                    uint16_t pred_adj = static_cast<uint16_t>((int)pred + C[ctx_hi]);

                    uint16_t u  = static_cast<uint16_t>(px - pred_adj);
                    uint16_t z  = zigzag_encode(u);
                    uint8_t  hi = z >> 8;
                    uint8_t  lo = z & 0xFF;

                    enc.write(hi_models[ctx_hi], hi);
                    hi_models[ctx_hi].increment(hi);

                    if (hi == 0) {
                        int gc = std::min(quantize(hi_W) + quantize(hi_N[col]), 3);
                        int lo_idx = gc * 8 + quantize_prev_lo(prev_lo_hi0);
                        enc.write(lo_hi0_models[lo_idx], lo);
                        lo_hi0_models[lo_idx].increment(lo);
                        prev_lo_hi0 = lo;
                    } else {
                        enc.write(lo_hip_models[hi-1], lo);
                        lo_hip_models[hi-1].increment(lo);
                    }

                    int r_s = (u <= 32767u) ? (int)u : (int)u - 65536;
                    B[ctx_hi] += r_s; Nc[ctx_hi]++;
                    if      (B[ctx_hi] >  Nc[ctx_hi]) { C[ctx_hi]++; B[ctx_hi] -= Nc[ctx_hi]; if (B[ctx_hi]>Nc[ctx_hi])  B[ctx_hi]=Nc[ctx_hi]; }
                    else if (B[ctx_hi] < -Nc[ctx_hi]) { C[ctx_hi]--; B[ctx_hi] += Nc[ctx_hi]; if (B[ctx_hi]<-Nc[ctx_hi]) B[ctx_hi]=-Nc[ctx_hi]; }
                    if (Nc[ctx_hi] == 64) { Nc[ctx_hi] >>= 1; B[ctx_hi] >>= 1; }

                    hi_W = hi; hi_N[col] = hi;
                }
            }

            enc.finish();
            tile_bufs[t] = oss.str();
        });
    }
    for (auto& th : threads) th.join();

    // ── Header ─────────────────────────────────────────────────────────────
    out_ptr->write(reinterpret_cast<const char*>(MAGIC), 4);
    write_u32le(*out_ptr, width);
    write_u32le(*out_ptr, height);
    write_u16le(*out_ptr, global_mean);
    out_ptr->write(reinterpret_cast<const char*>(hi_hist_gap.data()),    256);
    out_ptr->write(reinterpret_cast<const char*>(lo_hi0hist_gap.data()), 256);
    out_ptr->write(reinterpret_cast<const char*>(lo_hiphist_gap.data()), 256);
    out_ptr->write(reinterpret_cast<const char*>(hi_hist_mean.data()),    256);
    out_ptr->write(reinterpret_cast<const char*>(lo_hi0hist_mean.data()), 256);
    out_ptr->write(reinterpret_cast<const char*>(lo_hiphist_mean.data()), 256);
    out_ptr->put(static_cast<char>(static_cast<uint8_t>(ntiles)));
    out_ptr->write(reinterpret_cast<const char*>(tile_pred_modes.data()), ntiles);
    for (unsigned t = 0; t < ntiles; ++t)
        write_u32le(*out_ptr, static_cast<uint32_t>(tile_bufs[t].size()));
    for (unsigned t = 0; t < ntiles; ++t)
        out_ptr->write(tile_bufs[t].data(), static_cast<std::streamsize>(tile_bufs[t].size()));

    return 0;
}
