//
// ASTRA decompressor
//
// Reverses compress.cpp exactly:
//   1. Read header: background, T, anomaly list
//   2. Decode rANS delta stream
//   3. Walk pixels in order: anomaly positions get their stored value,
//      all other positions get background + delta
//   4. Write big-endian uint16 output
//
// Usage: ./decompress <input.astr> <output>
//

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cassert>

static constexpr uint8_t  MAGIC[4] = {'A','S','T','R'};
static constexpr uint16_t VERSION  = 1;

static constexpr uint32_t RANS_SCALE_BITS = 16;
static constexpr uint32_t RANS_SCALE      = 1u << RANS_SCALE_BITS;
static constexpr uint32_t RANS_L          = 1u << 23;

// ---------------------------------------------------------------------------
// rANS table + decoder (identical to compressor)
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
        ptr = data;
        state  = (uint32_t)(*ptr++) << 24;
        state |= (uint32_t)(*ptr++) << 16;
        state |= (uint32_t)(*ptr++) <<  8;
        state |= (uint32_t)(*ptr++);
    }

    uint8_t decode(const RansTable& tab) {
        uint32_t slot = state % RANS_SCALE;
        uint8_t  sym  = tab.sym_of_slot[slot];
        state = tab.freq[sym] * (state / RANS_SCALE) + slot - tab.cumul[sym];
        while (state < RANS_L) state = (state << 8) | (*ptr++);
        return sym;
    }
};

// ---------------------------------------------------------------------------
// Zigzag decode
// ---------------------------------------------------------------------------

static inline int32_t zagzig16(uint16_t u) {
    return (u & 1) ? -(int32_t)((u + 1) / 2) : (int32_t)(u / 2);
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static uint16_t read_u16le(FILE* f) {
    uint8_t b[2]; fread(b, 1, 2, f);
    return (uint16_t)(b[0] | ((uint16_t)b[1] << 8));
}
static uint32_t read_u32le(FILE* f) {
    uint8_t b[4]; fread(b, 1, 4, f);
    uint32_t v = 0; for (int i=3;i>=0;i--) v=(v<<8)|b[i]; return v;
}
static uint64_t read_u64le(FILE* f) {
    uint8_t b[8]; fread(b, 1, 8, f);
    uint64_t v = 0; for (int i=7;i>=0;i--) v=(v<<8)|b[i]; return v;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.astr> <output>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    // --- Read header ---
    uint8_t magic[4]; fread(magic, 1, 4, fin);
    if (memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n"); return 1;
    }
    uint16_t version = read_u16le(fin);
    if (version != VERSION) {
        fprintf(stderr, "Unknown version %d\n", version); return 1;
    }
    uint32_t width      = read_u32le(fin);
    uint32_t height     = read_u32le(fin);
    uint16_t background = read_u16le(fin);
    read_u16le(fin); // T — used by compressor only, not needed for reconstruction
    uint32_t n_anom     = read_u32le(fin);

    int npix = (int)(width * height);

    // --- Read anomaly list ---
    // Store as a map: pixel index → value
    // We use a sorted list and walk it in order during reconstruction
    struct Anomaly { uint32_t idx; uint16_t val; };
    std::vector<Anomaly> anomalies(n_anom);
    for (uint32_t i = 0; i < n_anom; i++) {
        anomalies[i].idx = read_u32le(fin);
        anomalies[i].val = read_u16le(fin);
    }
    // anomalies are already in pixel order (written left-to-right, top-to-bottom)

    // --- Read rANS streams ---
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

    // --- Decode delta symbols ---
    int n_bg = npix - (int)n_anom;
    std::vector<uint8_t> sym_hi(n_bg), sym_lo(n_bg);
    RansDecoder dec_hi, dec_lo;
    dec_hi.init(bits_hi.data());
    dec_lo.init(bits_lo.data());
    for (int i = 0; i < n_bg; i++) sym_hi[i] = dec_hi.decode(tab_hi);
    for (int i = 0; i < n_bg; i++) sym_lo[i] = dec_lo.decode(tab_lo);

    // --- Reconstruct pixels ---
    // Walk pixels in order. Use anomaly list as a queue.
    std::vector<uint16_t> pixels(npix);
    int bg_idx   = 0;  // index into decoded delta stream
    int anom_idx = 0;  // index into anomaly list

    for (int i = 0; i < npix; i++) {
        if (anom_idx < (int)n_anom && anomalies[anom_idx].idx == (uint32_t)i) {
            // anomaly: use stored value directly
            pixels[i] = anomalies[anom_idx].val;
            anom_idx++;
        } else {
            // background: recover delta, add to background
            uint16_t zz    = ((uint16_t)sym_hi[bg_idx] << 8) | sym_lo[bg_idx];
            int32_t  delta = zagzig16(zz);
            pixels[i]      = (uint16_t)((int32_t)background + delta);
            bg_idx++;
        }
    }

    // --- Write output (big-endian uint16) ---
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output %s\n", argv[2]); return 1; }
    for (int i = 0; i < npix; i++) {
        uint8_t b[2] = {(uint8_t)(pixels[i] >> 8), (uint8_t)(pixels[i] & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
