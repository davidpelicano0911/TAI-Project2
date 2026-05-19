/*
 * TAI Project 2 — Astronomical Image Compressor (ox-astro-fast)
 *
 * Optimised for minimum latency.
 *
 * Pipeline:
 *   1. Read raw bytes; interpret as 2D big-endian uint16_t array.
 *   2. Horizontal delta prediction: pred = left neighbour (or 0 at row start).
 *   3. Wrapping residual → modular zigzag.
 *   4. Split into hi and lo bytes.
 *   5. Encode hi with a single ORDER-0 adaptive model.
 *      Encode lo with a single ORDER-0 adaptive model.
 *
 * Two models total (vs 257 for balanced, 512 for ratio) → excellent cache
 * behaviour and minimal per-pixel overhead.
 *
 * Compressed file format
 * ──────────────────────
 *  Bytes  0–3  : magic "TA2F"
 *  Bytes  4–7  : width   uint32_t LE
 *  Bytes  8–11 : height  uint32_t LE
 *  Bytes 12+   : range-coded bitstream
 *
 * Usage:  compress_astro_fast <input_file> <output_file>
 *         compress_astro_fast          (stdin → stdout)
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

static constexpr uint8_t MAGIC[4] = {'T','A','2','F'};

static void write_u32le(std::ostream& out, uint32_t v) {
    out.put(static_cast<char>(v         & 0xFF));
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
        if (!fin)  { std::cerr << "Cannot open input: "  << argv[1] << '\n'; return 1; }
        fout.open(argv[2], std::ios::binary);
        if (!fout) { std::cerr << "Cannot open output: " << argv[2] << '\n'; return 1; }
        in_ptr  = &fin;
        out_ptr = &fout;
    } else if (argc == 1) {
        std::cin.sync_with_stdio(false);
    } else {
        std::cerr << "Usage: compress_astro_fast <input> <output>\n"; return 1;
    }

    std::vector<uint8_t> raw(std::istreambuf_iterator<char>(*in_ptr), {});
    if (raw.size() % 2 != 0) {
        std::cerr << "Input size not even\n"; return 1;
    }

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

    out_ptr->write(reinterpret_cast<const char*>(MAGIC), 4);
    write_u32le(*out_ptr, width);
    write_u32le(*out_ptr, height);

    // Two ORDER-0 models: one for hi bytes, one for lo bytes
    SimpleFrequencyTable hi_model = make_flat256();
    SimpleFrequencyTable lo_model = make_flat256();

    RangeEncoder enc(*out_ptr);

    for (uint32_t row = 0; row < height; ++row) {
        for (uint32_t col = 0; col < width; ++col) {
            uint16_t px   = pixels[row * width + col];
            // Horizontal delta: predict from left neighbour; restart at each row
            uint16_t pred = (col > 0) ? pixels[row * width + col - 1] : 0u;
            uint16_t u    = zigzag_encode(static_cast<uint16_t>(px - pred));

            uint8_t hi = static_cast<uint8_t>(u >> 8);
            uint8_t lo = static_cast<uint8_t>(u & 0xFF);

            enc.write(hi_model, hi);
            hi_model.increment(hi);

            enc.write(lo_model, lo);
            lo_model.increment(lo);
        }
    }

    enc.finish();
    return 0;
}
