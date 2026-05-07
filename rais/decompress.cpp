//
// TAI Project 2 — Decompressor (matches compress.cpp format)
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cassert>

static constexpr uint8_t  MAGIC[4]       = {'R','A','I','S'};
static constexpr uint16_t VERSION        = 1;
static constexpr uint32_t RANS_SCALE_BITS = 16;
static constexpr uint32_t RANS_SCALE      = 1u << RANS_SCALE_BITS;
static constexpr uint32_t RANS_L          = 1u << 23;

// ---------------------------------------------------------------------------
// Predictor (identical to compressor)
// ---------------------------------------------------------------------------

static inline uint16_t predict(const uint16_t* row, const uint16_t* prev_row,
                                int x, int width) {
    if (prev_row == nullptr && x == 0) return 0;
    if (prev_row == nullptr)           return row[x - 1];
    if (x == 0)                        return prev_row[x];

    uint16_t a = row[x - 1];
    uint16_t b = prev_row[x];
    uint16_t c = prev_row[x - 1];

    int pred = (int)a + (int)b - (int)c;
    int lo   = std::min((int)a, (int)b);
    int hi   = std::max((int)a, (int)b);
    if      (pred < lo) pred = lo;
    else if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

static inline int16_t zigzag_dec(uint16_t v) {
    return (v & 1) ? -(int16_t)((v + 1) / 2) : (int16_t)(v / 2);
}

// ---------------------------------------------------------------------------
// rANS table + decoder
// ---------------------------------------------------------------------------

struct RansTable {
    uint32_t freq[256];
    uint32_t cumul[257];
    uint8_t  sym_of_slot[RANS_SCALE];

    void build_cumul() {
        cumul[0] = 0;
        for (int i = 0; i < 256; i++) cumul[i+1] = cumul[i] + freq[i];
    }
    void build_lookup() {
        for (int s = 0; s < 256; s++)
            for (uint32_t j = cumul[s]; j < cumul[s+1]; j++)
                sym_of_slot[j] = (uint8_t)s;
    }
};

struct RansDecoder {
    uint32_t state;
    const uint8_t* ptr;

    void init(const uint8_t* data) {
        ptr    = data;
        state  = (uint32_t)(*ptr++) << 24;
        state |= (uint32_t)(*ptr++) << 16;
        state |= (uint32_t)(*ptr++) <<  8;
        state |= (uint32_t)(*ptr++);
    }

    uint8_t decode(const RansTable& tab) {
        uint32_t slot = state % RANS_SCALE;
        uint8_t  sym  = tab.sym_of_slot[slot];
        uint32_t f    = tab.freq[sym];
        uint32_t c    = tab.cumul[sym];
        state = f * (state / RANS_SCALE) + slot - c;
        while (state < RANS_L) state = (state << 8) | (*ptr++);
        return sym;
    }
};

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static uint16_t read_u16le(FILE* f) {
    uint8_t b[2]; fread(b, 1, 2, f);
    return (uint16_t)(b[0] | ((uint16_t)b[1] << 8));
}
static uint32_t read_u32le(FILE* f) {
    uint8_t b[4]; fread(b, 1, 4, f);
    uint32_t v = 0;
    for (int i = 3; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}
static uint64_t read_u64le(FILE* f) {
    uint8_t b[8]; fread(b, 1, 8, f);
    uint64_t v = 0;
    for (int i = 7; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.rais> <output>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open: %s\n", argv[1]); return 1; }

    uint8_t magic[4]; fread(magic, 1, 4, fin);
    if (memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n"); fclose(fin); return 1;
    }
    uint16_t version = read_u16le(fin);
    if (version != VERSION) {
        fprintf(stderr, "Unknown version %d\n", version); fclose(fin); return 1;
    }
    uint32_t width    = read_u32le(fin);
    uint32_t height   = read_u32le(fin);
    uint32_t nstreams = read_u32le(fin);
    if (nstreams != 2) {
        fprintf(stderr, "Expected 2 streams\n"); fclose(fin); return 1;
    }

    int npix = (int)(width * height);

    auto read_stream = [&](RansTable& tab, std::vector<uint8_t>& bits) {
        uint32_t nsyms = read_u32le(fin); (void)nsyms;
        for (int i = 0; i < 256; i++) tab.freq[i] = read_u32le(fin);
        tab.build_cumul();
        tab.build_lookup();
        uint64_t nbytes = read_u64le(fin);
        bits.resize(nbytes);
        fread(bits.data(), 1, nbytes, fin);
    };

    RansTable tab_hi, tab_lo;
    std::vector<uint8_t> bits_hi, bits_lo;
    read_stream(tab_hi, bits_hi);
    read_stream(tab_lo, bits_lo);
    fclose(fin);

    std::vector<uint8_t> hi(npix), lo(npix);
    RansDecoder dec_hi, dec_lo;
    dec_hi.init(bits_hi.data());
    dec_lo.init(bits_lo.data());
    for (int i = 0; i < npix; i++) hi[i] = dec_hi.decode(tab_hi);
    for (int i = 0; i < npix; i++) lo[i] = dec_lo.decode(tab_lo);

    std::vector<uint16_t> pixels(npix);
    for (int i = 0; i < npix; i++) {
        int x = i % (int)width, y = i / (int)width;
        const uint16_t* row      = pixels.data() + y * width;
        const uint16_t* prev_row = (y > 0) ? (pixels.data() + (y-1) * width) : nullptr;
        uint16_t pred = predict(row, prev_row, x, (int)width);
        uint16_t zz   = ((uint16_t)hi[i] << 8) | lo[i];
        int16_t  res  = zigzag_dec(zz);
        pixels[i]     = (uint16_t)((int)pred + res);
    }

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", argv[2]); return 1; }
    for (int i = 0; i < npix; i++) {
        uint8_t b[2] = {(uint8_t)(pixels[i] >> 8), (uint8_t)(pixels[i] & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
