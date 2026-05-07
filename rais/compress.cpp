//
// TAI Project 2 — Speed-focused compressor
//
// Pipeline:
//   1. MED predictor  (exploits 2-D spatial smoothness of astro images)
//   2. Zigzag map     signed residual → unsigned (0,1,2,3,... for 0,-1,1,-2,2,...)
//   3. Split high/low bytes of the mapped residual
//   4. rANS entropy coding of each byte stream independently
//
// Format (binary, little-endian integers):
//   [4]  magic  "RAIS"
//   [2]  version = 1
//   [4]  width  (pixels)
//   [4]  height (pixels)
//   [4]  n_streams  (= 2: high-byte stream, low-byte stream)
//   per stream:
//     [4]  nsyms  (= 256)
//     [256*4] freq table (uint32, sum = total pixels)
//     [8]  compressed_bytes (uint64)
//     [N]  rANS bitstream (written MSB→LSB, decoded LSB→MSB)
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

static constexpr int     WIDTH_DEFAULT  = 1500;
static constexpr int     HEIGHT_DEFAULT = 1500;
static constexpr uint8_t MAGIC[4]       = {'R','A','I','S'};
static constexpr uint16_t VERSION       = 1;

// rANS parameters
static constexpr uint32_t RANS_SCALE_BITS = 16;
static constexpr uint32_t RANS_SCALE      = 1u << RANS_SCALE_BITS;  // M = 65536
static constexpr uint32_t RANS_L          = 1u << 23;               // lower bound

// ---------------------------------------------------------------------------
// Predictor
// ---------------------------------------------------------------------------

// MED (Median Edge Detection) predictor — same as used in LOCO-I/JPEG-LS
static inline uint16_t predict(const uint16_t* row, const uint16_t* prev_row,
                                int x, int width) {
    if (prev_row == nullptr && x == 0) return 0;
    if (prev_row == nullptr)           return row[x - 1];
    if (x == 0)                        return prev_row[x];

    uint16_t a = row[x - 1];        // left
    uint16_t b = prev_row[x];       // up
    uint16_t c = prev_row[x - 1];   // up-left

    int pred = (int)a + (int)b - (int)c;
    // clamp to [min(a,b), max(a,b)]
    int lo = std::min((int)a, (int)b);
    int hi = std::max((int)a, (int)b);
    if      (pred < lo) pred = lo;
    else if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

// ---------------------------------------------------------------------------
// Zigzag: signed int16 residual → uint16
// ---------------------------------------------------------------------------

static inline uint16_t zigzag_enc(int16_t v) {
    return (uint16_t)((v >= 0) ? (v * 2) : ((-v) * 2 - 1));
}

static inline int16_t zigzag_dec(uint16_t v) {
    return (v & 1) ? -(int16_t)((v + 1) / 2) : (int16_t)(v / 2);
}

// ---------------------------------------------------------------------------
// rANS — table-based, alias-free
//
// We use a standard static rANS with a cumulative frequency table.
// The freq table is normalised to RANS_SCALE (power of 2) for fast decoding.
// ---------------------------------------------------------------------------

struct RansTable {
    uint32_t freq[256];   // normalised frequencies, sum = RANS_SCALE
    uint32_t cumul[257];  // cumulative: cumul[s] = sum freq[0..s-1]

    // symbol lookup by slot: slot -> symbol (for decoder)
    uint8_t  sym_of_slot[RANS_SCALE]; // indexed by (state % RANS_SCALE)

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

// Normalise raw counts to sum = RANS_SCALE, ensuring freq >= 1 for used symbols
static void normalise_freqs(const uint64_t* counts, uint32_t* freq, int nsyms) {
    uint64_t total = 0;
    for (int i = 0; i < nsyms; i++) total += counts[i];

    // Assign at least 1 to every symbol that appears
    uint32_t assigned = 0;
    for (int i = 0; i < nsyms; i++) {
        if (counts[i] == 0) { freq[i] = 0; continue; }
        freq[i] = (uint32_t)std::max((uint64_t)1,
                      (counts[i] * (uint64_t)RANS_SCALE) / total);
        assigned += freq[i];
    }

    // Fix rounding so sum == RANS_SCALE
    // Add/remove from the most frequent symbol
    int best = 0;
    for (int i = 1; i < nsyms; i++)
        if (counts[i] > counts[best]) best = i;

    if (assigned < RANS_SCALE)      freq[best] += RANS_SCALE - assigned;
    else if (assigned > RANS_SCALE) freq[best] -= assigned - RANS_SCALE;
}

// rANS encoder — emits bytes to a vector (reversed at the end)
struct RansEncoder {
    uint32_t state;
    std::vector<uint8_t> buf; // bytes in reverse order

    RansEncoder() : state(RANS_L) {}

    void encode(uint8_t sym, const RansTable& tab) {
        uint32_t f = tab.freq[sym];
        uint32_t c = tab.cumul[sym];

        // Renormalise: push bytes until state is in [L*f/M, L*f/M * 256)
        uint32_t x_max = ((RANS_L / RANS_SCALE) * 256) * f;
        while (state >= x_max) {
            buf.push_back((uint8_t)(state & 0xFF));
            state >>= 8;
        }
        // Encode
        state = (state / f) * RANS_SCALE + c + (state % f);
    }

    // Flush final state (4 bytes, big-endian for stream)
    void flush() {
        buf.push_back((uint8_t)( state        & 0xFF));
        buf.push_back((uint8_t)((state >>  8) & 0xFF));
        buf.push_back((uint8_t)((state >> 16) & 0xFF));
        buf.push_back((uint8_t)((state >> 24) & 0xFF));
        // reverse so the decoder can read forward
        std::reverse(buf.begin(), buf.end());
    }
};

// rANS decoder — reads forward from a byte pointer
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
        uint32_t f    = tab.freq[sym];
        uint32_t c    = tab.cumul[sym];

        // Update state
        state = f * (state / RANS_SCALE) + slot - c;

        // Renormalise: pull bytes
        while (state < RANS_L) {
            state = (state << 8) | (*ptr++);
        }
        return sym;
    }
};

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 4, f);
}
static void write_u64le(FILE* f, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) { b[i] = v & 0xFF; v >>= 8; }
    fwrite(b, 1, 8, f);
}

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
// Compress
// ---------------------------------------------------------------------------

static int compress(const char* in_path, const char* out_path) {
    FILE* fin = fopen(in_path, "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", in_path); return 1; }

    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    rewind(fin);

    int width  = WIDTH_DEFAULT;
    int height = HEIGHT_DEFAULT;
    int npix   = width * height;

    if (fsize != (long)npix * 2) {
        fprintf(stderr, "Unexpected file size %ld (expected %d)\n", fsize, npix * 2);
        fclose(fin);
        return 1;
    }

    // Read raw big-endian uint16
    std::vector<uint16_t> pixels(npix);
    for (int i = 0; i < npix; i++) {
        uint8_t b[2]; fread(b, 1, 2, fin);
        pixels[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    // --- Predict + zigzag ---
    std::vector<uint8_t> hi(npix), lo(npix);
    for (int i = 0; i < npix; i++) {
        int x = i % width, y = i / width;
        const uint16_t* row      = pixels.data() + y * width;
        const uint16_t* prev_row = (y > 0) ? (pixels.data() + (y-1) * width) : nullptr;
        uint16_t pred = predict(row, prev_row, x, width);
        int16_t  res  = (int16_t)(pixels[i] - pred);
        uint16_t zz   = zigzag_enc(res);
        hi[i] = (uint8_t)(zz >> 8);
        lo[i] = (uint8_t)(zz & 0xFF);
    }

    // --- Build frequency tables ---
    uint64_t cnt_hi[256] = {}, cnt_lo[256] = {};
    for (int i = 0; i < npix; i++) { cnt_hi[hi[i]]++; cnt_lo[lo[i]]++; }

    RansTable tab_hi, tab_lo;
    normalise_freqs(cnt_hi, tab_hi.freq, 256);
    normalise_freqs(cnt_lo, tab_lo.freq, 256);
    tab_hi.build_cumul(); tab_hi.build_lookup();
    tab_lo.build_cumul(); tab_lo.build_lookup();

    // --- Encode (symbols must be fed in REVERSE order for rANS) ---
    RansEncoder enc_hi, enc_lo;
    for (int i = npix - 1; i >= 0; i--) enc_hi.encode(hi[i], tab_hi);
    for (int i = npix - 1; i >= 0; i--) enc_lo.encode(lo[i], tab_lo);
    enc_hi.flush();
    enc_lo.flush();

    // --- Write output ---
    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u16le(fout, VERSION);
    write_u32le(fout, (uint32_t)width);
    write_u32le(fout, (uint32_t)height);
    write_u32le(fout, 2); // n_streams

    // Stream 0: high bytes
    write_u32le(fout, 256);
    for (int i = 0; i < 256; i++) write_u32le(fout, tab_hi.freq[i]);
    write_u64le(fout, (uint64_t)enc_hi.buf.size());
    fwrite(enc_hi.buf.data(), 1, enc_hi.buf.size(), fout);

    // Stream 1: low bytes
    write_u32le(fout, 256);
    for (int i = 0; i < 256; i++) write_u32le(fout, tab_lo.freq[i]);
    write_u64le(fout, (uint64_t)enc_lo.buf.size());
    fwrite(enc_lo.buf.data(), 1, enc_lo.buf.size(), fout);

    fclose(fout);

    long out_size = (long)(4 + 2 + 4 + 4 + 4
                        + (4 + 256*4 + 8 + enc_hi.buf.size())
                        + (4 + 256*4 + 8 + enc_lo.buf.size()));
    fprintf(stderr, "Compressed: %ld -> %ld bytes (%.4f bpp)\n",
            fsize, out_size, (double)out_size * 8.0 / npix);
    return 0;
}

// ---------------------------------------------------------------------------
// Decompress
// ---------------------------------------------------------------------------

static int decompress(const char* in_path, const char* out_path) {
    FILE* fin = fopen(in_path, "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", in_path); return 1; }

    // Check magic
    uint8_t magic[4]; fread(magic, 1, 4, fin);
    if (memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic bytes\n"); fclose(fin); return 1;
    }
    uint16_t version = read_u16le(fin);
    if (version != VERSION) {
        fprintf(stderr, "Unknown version %d\n", version); fclose(fin); return 1;
    }
    uint32_t width    = read_u32le(fin);
    uint32_t height   = read_u32le(fin);
    uint32_t nstreams = read_u32le(fin);
    if (nstreams != 2) {
        fprintf(stderr, "Expected 2 streams, got %d\n", nstreams); fclose(fin); return 1;
    }

    int npix = (int)(width * height);

    auto read_stream = [&](RansTable& tab, std::vector<uint8_t>& bits) {
        uint32_t nsyms = read_u32le(fin);
        assert(nsyms == 256);
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

    // --- Decode ---
    std::vector<uint8_t> hi(npix), lo(npix);
    RansDecoder dec_hi, dec_lo;
    dec_hi.init(bits_hi.data());
    dec_lo.init(bits_lo.data());
    for (int i = 0; i < npix; i++) hi[i] = dec_hi.decode(tab_hi);
    for (int i = 0; i < npix; i++) lo[i] = dec_lo.decode(tab_lo);

    // --- Reconstruct pixels ---
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

    // --- Write big-endian output ---
    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }
    for (int i = 0; i < npix; i++) {
        uint8_t b[2] = {(uint8_t)(pixels[i] >> 8), (uint8_t)(pixels[i] & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input> <output>\n", argv[0]);
        return 1;
    }
    return compress(argv[1], argv[2]);
}
