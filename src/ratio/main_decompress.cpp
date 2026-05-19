/*
 * TAI Project 2 — Astronomical Image Decompressor (ox-astro-balanced)
 *
 * Inverse of main_compress.cpp:
 *   1. Read header: magic "TA2A", width, height.
 *   2. Decode hi/lo bytes per pixel via ORDER-1 range decoding.
 *   3. Reconstruct zigzag uint16_t residuals from byte pairs.
 *   4. Inverse-zigzag → signed residual; add JPEG-LS MED prediction.
 *   5. Write pixels as big-endian uint16_t.
 *
 * Usage:  decompress_astro <input_file> <output_file>
 *         decompress_astro          (stdin → stdout)
 */

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "coder/RangeCoder.hpp"
#include "coder/FrequencyTable.hpp"
#include "model/ImagePredictor.hpp"

static constexpr uint8_t MAGIC[4] = {'T','A','2','A'};

static uint32_t read_u32le(std::istream& in) {
    uint32_t v = 0;
    for (int i = 0; i < 4; ++i)
        v |= static_cast<uint32_t>(static_cast<uint8_t>(in.get())) << (8 * i);
    return v;
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
        std::cerr << "Usage: decompress_astro <input> <output>\n"; return 1;
    }

    // Verify magic
    uint8_t magic[4];
    in_ptr->read(reinterpret_cast<char*>(magic), 4);
    for (int i = 0; i < 4; ++i) {
        if (magic[i] != MAGIC[i]) {
            std::cerr << "Bad magic number\n"; return 1;
        }
    }

    const uint32_t width  = read_u32le(*in_ptr);
    const uint32_t height = read_u32le(*in_ptr);
    const uint64_t npix   = static_cast<uint64_t>(width) * height;

    std::vector<uint16_t> pixels(npix);

    std::vector<SimpleFrequencyTable> hi_models, lo_models;
    hi_models.reserve(256);
    lo_models.reserve(256);
    for (int i = 0; i < 256; ++i) {
        hi_models.push_back(make_flat256());
        lo_models.push_back(make_flat256());
    }

    RangeDecoder dec(*in_ptr);

    uint8_t prev_hi = 0;
    for (uint32_t row = 0; row < height; ++row) {
        for (uint32_t col = 0; col < width; ++col) {
            uint8_t hi = static_cast<uint8_t>(dec.read(hi_models[prev_hi]));
            hi_models[prev_hi].increment(hi);

            uint8_t lo = static_cast<uint8_t>(dec.read(lo_models[hi]));
            lo_models[hi].increment(lo);

            uint16_t z    = (static_cast<uint16_t>(hi) << 8) | lo;
            uint16_t u    = zigzag_decode(z);

            uint16_t W  = (col  > 0)             ? pixels[row * width + col - 1]       : 0u;
            uint16_t N  = (row  > 0)             ? pixels[(row-1) * width + col]        : W;
            uint16_t NW = (row  > 0 && col  > 0) ? pixels[(row-1) * width + col - 1]   : W;

            uint16_t pred = (row == 0 && col == 0) ? 0u : med_predict(W, N, NW);
            pixels[row * width + col] = static_cast<uint16_t>(pred + u);

            prev_hi = hi;
        }
    }

    // Write as big-endian uint16_t
    std::vector<uint8_t> out_bytes(npix * 2);
    for (uint64_t i = 0; i < npix; ++i) {
        out_bytes[2*i]   = static_cast<uint8_t>(pixels[i] >> 8);
        out_bytes[2*i+1] = static_cast<uint8_t>(pixels[i] & 0xFF);
    }
    out_ptr->write(reinterpret_cast<const char*>(out_bytes.data()),
                   static_cast<std::streamsize>(out_bytes.size()));
    return 0;
}
