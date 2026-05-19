/*
 * TAI Project 2 — Astronomical Image Compressor (ox-astro-balanced)
 *
 * Specialized for raw 16-bit big-endian raster images (e.g. 1500×1500).
 *
 * Pipeline:
 *   1. Read raw bytes; interpret as 2D array of big-endian uint16_t.
 *   2. Apply JPEG-LS MED spatial predictor → signed residuals.
 *   3. Zigzag-encode residuals to uint16_t (concentrates mass near 0).
 *   4. Split each uint16_t into high byte (bits[15:8]) and low byte (bits[7:0]).
 *   5. Encode interleaved per pixel:
 *        hi  with hi_model[prev_hi]  (ORDER-1 on previous high byte)
 *        lo  with lo_model[hi]       (conditioned on current high byte)
 *   6. Single adaptive range encoder across all pixels.
 *
 * Compressed file format
 * ──────────────────────
 *  Bytes  0–3  : magic "TA2A"
 *  Bytes  4–7  : width   uint32_t LE
 *  Bytes  8–11 : height  uint32_t LE
 *  Bytes 12+   : range-coded bitstream (flush with finish())
 *
 * Usage:  compress_astro <input_file> <output_file>
 *         compress_astro          (stdin → stdout)
 */

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "coder/RangeCoder.hpp"
#include "coder/FrequencyTable.hpp"
#include "model/ImagePredictor.hpp"

static constexpr uint8_t MAGIC[4] = {'T','A','2','A'};

static void write_u32le(std::ostream& out, uint32_t v) {
    out.put(static_cast<char>(v        & 0xFF));
    out.put(static_cast<char>((v >>  8) & 0xFF));
    out.put(static_cast<char>((v >> 16) & 0xFF));
    out.put(static_cast<char>((v >> 24) & 0xFF));
}

static SimpleFrequencyTable make_flat256() {
    return SimpleFrequencyTable(std::vector<uint32_t>(256, 1u));
}

int main(int argc, char* argv[]) {
    std::istream* in_ptr  = &std::cin;
    std::ostream* out_ptr = &std::cout;
    std::ifstream fin;
    std::ofstream fout;

    if (argc >= 3) {
        fin.open(argv[1], std::ios::binary);
        if (!fin) { std::cerr << "Cannot open input: " << argv[1] << '\n'; return 1; }
        fout.open(argv[2], std::ios::binary);
        if (!fout) { std::cerr << "Cannot open output: " << argv[2] << '\n'; return 1; }
        in_ptr  = &fin;
        out_ptr = &fout;
    } else if (argc == 1) {
        std::cin.sync_with_stdio(false);
    } else {
        std::cerr << "Usage: compress_astro <input> <output>\n"; return 1;
    }

    // Read entire input
    std::vector<uint8_t> raw(std::istreambuf_iterator<char>(*in_ptr), {});
    if (raw.size() % 2 != 0) {
        std::cerr << "Input size not even — not a valid 16-bit raster\n"; return 1;
    }

    // Infer dimensions: expect square image, or rectangular if exact match
    const uint64_t npix = raw.size() / 2;
    uint32_t width = 0, height = 0;
    // Try square first
    uint32_t sq = static_cast<uint32_t>(std::sqrt(static_cast<double>(npix)));
    if ((uint64_t)sq * sq == npix) {
        width = height = sq;
    } else {
        // Common rectangular sizes: try 1500×1500, 1000×1000, 2048×1024, etc.
        // Fall back to 1×npix (1D mode)
        width  = static_cast<uint32_t>(npix);
        height = 1;
        // Try common widths
        for (uint32_t w : {1500u, 2048u, 1024u, 512u, 256u}) {
            if (npix % w == 0) { width = w; height = static_cast<uint32_t>(npix / w); break; }
        }
    }

    // Decode big-endian uint16_t pixels into host order
    std::vector<uint16_t> pixels(npix);
    for (uint64_t i = 0; i < npix; ++i)
        pixels[i] = (static_cast<uint16_t>(raw[2*i]) << 8) | raw[2*i+1];

    // Write header
    out_ptr->write(reinterpret_cast<const char*>(MAGIC), 4);
    write_u32le(*out_ptr, width);
    write_u32le(*out_ptr, height);

    // 256 adaptive models for hi bytes (ORDER-1: context = prev hi byte)
    // 256 adaptive models for lo bytes (conditioned on current hi byte)
    std::vector<SimpleFrequencyTable> hi_models, lo_models;
    hi_models.reserve(256);
    lo_models.reserve(256);
    for (int i = 0; i < 256; ++i) {
        hi_models.push_back(make_flat256());
        lo_models.push_back(make_flat256());
    }

    RangeEncoder enc(*out_ptr);

    uint8_t prev_hi = 0;
    for (uint32_t row = 0; row < height; ++row) {
        for (uint32_t col = 0; col < width; ++col) {
            uint16_t px = pixels[row * width + col];

            // Prediction
            uint16_t W  = (col  > 0)              ? pixels[row * width + col - 1]           : 0;
            uint16_t N  = (row  > 0)              ? pixels[(row-1) * width + col]            : W;
            uint16_t NW = (row  > 0 && col  > 0)  ? pixels[(row-1) * width + col - 1]       : W;

            uint16_t pred = (row == 0 && col == 0) ? 0u : med_predict(W, N, NW);
            uint16_t u    = zigzag_encode(static_cast<uint16_t>(px - pred));

            uint8_t hi = static_cast<uint8_t>(u >> 8);
            uint8_t lo = static_cast<uint8_t>(u & 0xFF);

            enc.write(hi_models[prev_hi], hi);
            hi_models[prev_hi].increment(hi);

            enc.write(lo_models[hi], lo);
            lo_models[hi].increment(lo);

            prev_hi = hi;
        }
    }

    enc.finish();
    return 0;
}
