//
// ASTRA compressor
//
// Model: background + exceptions
//
//   Every pixel is classified as either:
//     - background: |pixel - background| <= T  → store small signed delta
//     - anomaly:    |pixel - background| >  T  → store (position, value) explicitly
//
//   Background and T are computed automatically from the image histogram:
//     - background = mode (most frequent pixel value)
//     - T          = 4 * MAD (median absolute deviation), rounded to next power of 2
//
//   The deltas for background pixels are entropy-coded with rANS.
//   Anomalies are stored as raw (uint32 index, uint16 value) pairs — there are
//   very few of them in astronomical images.
//
// Format:
//   [4]  magic "ASTR"
//   [2]  version = 1
//   [4]  width
//   [4]  height
//   [2]  background value
//   [2]  threshold T
//   [4]  n_anomalies
//   [n_anomalies * 6]  anomaly list: each entry is [uint32 pixel_index][uint16 value]
//   [4]  nsyms = 256
//   [256*4]  rANS frequency table for deltas (mapped to 0..255 via zigzag)
//   [8]  compressed delta stream size in bytes
//   [N]  rANS compressed delta stream
//
// Usage: ./compress <input> <output.astr>
//

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>

// ---------------------------------------------------------------------------
// Image constants
// ---------------------------------------------------------------------------

static constexpr int WIDTH  = 1500;
static constexpr int HEIGHT = 1500;
static constexpr int NPIX   = WIDTH * HEIGHT;

static constexpr uint8_t  MAGIC[4] = {'A','S','T','R'};
static constexpr uint16_t VERSION  = 1;

// ---------------------------------------------------------------------------
// rANS constants
// ---------------------------------------------------------------------------

static constexpr uint32_t RANS_SCALE_BITS = 16;
static constexpr uint32_t RANS_SCALE      = 1u << RANS_SCALE_BITS; // 65536
static constexpr uint32_t RANS_L          = 1u << 23;

// ---------------------------------------------------------------------------
// Step 1 — find background (mode) from histogram
// ---------------------------------------------------------------------------

static uint16_t find_background(const uint32_t* hist) {
    uint16_t bg = 0;
    uint32_t best = 0;
    for (int v = 0; v < 65536; v++) {
        if (hist[v] > best) { best = hist[v]; bg = (uint16_t)v; }
    }
    return bg;
}

// ---------------------------------------------------------------------------
// Step 2 — compute MAD and auto-select T
// ---------------------------------------------------------------------------

static uint32_t find_threshold(const uint32_t* hist, uint16_t background) {
    // Build histogram of |pixel - background|
    std::vector<uint32_t> dev_hist(65536, 0);
    for (int v = 0; v < 65536; v++) {
        if (hist[v] == 0) continue;
        int dev = std::abs(v - (int)background);
        dev_hist[dev] += hist[v];
    }

    // Median of absolute deviations
    uint32_t mad = 0;
    uint64_t cumsum = 0;
    for (int d = 0; d < 65536; d++) {
        cumsum += dev_hist[d];
        if (cumsum >= (uint64_t)NPIX / 2) { mad = (uint32_t)d; break; }
    }

    // T = 4 * MAD, rounded up to next power of two, minimum 1
    uint32_t T = std::max(1u, 4u * mad);
    uint32_t pow2 = 1;
    while (pow2 < T) pow2 <<= 1;
    return pow2;
}

// ---------------------------------------------------------------------------
// Zigzag: map signed delta → unsigned symbol for rANS
//   0 → 0,  -1 → 1,  +1 → 2,  -2 → 3,  +2 → 4, ...
//
// Deltas are in range [-T, +T]. We map them to [0, 2T] via zigzag.
// Then we split into low 8 bits (lo) for rANS and keep hi separate if T > 127.
// For simplicity: if T <= 127 (delta fits in 8-bit zigzag), use 1 byte.
// Otherwise fall back to 2 bytes split into hi+lo streams.
// ---------------------------------------------------------------------------

static inline uint16_t zigzag16(int32_t v) {
    return (uint16_t)((v >= 0) ? (v * 2) : ((-v) * 2 - 1));
}

static inline int32_t zagzig16(uint16_t u) {
    return (u & 1) ? -(int32_t)((u + 1) / 2) : (int32_t)(u / 2);
}

// ---------------------------------------------------------------------------
// rANS table
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

static void normalise_freqs(const uint64_t* counts, uint32_t* freq, int nsyms) {
    uint64_t total = 0;
    for (int i = 0; i < nsyms; i++) total += counts[i];

    uint32_t assigned = 0;
    for (int i = 0; i < nsyms; i++) {
        if (counts[i] == 0) { freq[i] = 0; continue; }
        freq[i] = (uint32_t)std::max((uint64_t)1,
                    (counts[i] * (uint64_t)RANS_SCALE) / total);
        assigned += freq[i];
    }
    // fix rounding on the most frequent symbol
    int best = 0;
    for (int i = 1; i < nsyms; i++)
        if (counts[i] > counts[best]) best = i;
    if (assigned < RANS_SCALE)      freq[best] += RANS_SCALE - assigned;
    else if (assigned > RANS_SCALE) freq[best] -= assigned - RANS_SCALE;
}

// ---------------------------------------------------------------------------
// rANS encoder
// ---------------------------------------------------------------------------

struct RansEncoder {
    uint32_t state = RANS_L;
    std::vector<uint8_t> buf;

    void encode(uint8_t sym, const RansTable& tab) {
        uint32_t f = tab.freq[sym];
        uint32_t c = tab.cumul[sym];
        uint32_t x_max = ((RANS_L / RANS_SCALE) * 256) * f;
        while (state >= x_max) { buf.push_back(state & 0xFF); state >>= 8; }
        state = (state / f) * RANS_SCALE + c + (state % f);
    }

    void flush() {
        buf.push_back( state        & 0xFF);
        buf.push_back((state >>  8) & 0xFF);
        buf.push_back((state >> 16) & 0xFF);
        buf.push_back((state >> 24) & 0xFF);
        std::reverse(buf.begin(), buf.end());
    }
};

// ---------------------------------------------------------------------------
// rANS decoder
// ---------------------------------------------------------------------------

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
// I/O helpers
// ---------------------------------------------------------------------------

static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}
static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4]; for (int i=0;i<4;i++){b[i]=v&0xFF;v>>=8;} fwrite(b,1,4,f);
}
static void write_u64le(FILE* f, uint64_t v) {
    uint8_t b[8]; for (int i=0;i<8;i++){b[i]=v&0xFF;v>>=8;} fwrite(b,1,8,f);
}

// ---------------------------------------------------------------------------
// Compress
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input> <output.astr>\n", argv[0]);
        return 1;
    }

    // --- Read image ---
    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    std::vector<uint16_t> px(NPIX);
    for (int i = 0; i < NPIX; i++) {
        uint8_t b[2];
        if (fread(b, 1, 2, fin) != 2) { fprintf(stderr, "Short read\n"); return 1; }
        px[i] = (uint16_t)((b[0] << 8) | b[1]);
    }
    fclose(fin);

    // --- Step 1: histogram ---
    uint32_t hist[65536] = {};
    for (int i = 0; i < NPIX; i++) hist[px[i]]++;

    // --- Step 2: background + threshold ---
    uint16_t background = find_background(hist);
    uint32_t T          = find_threshold(hist, background);

    // --- Step 3: classify pixels ---
    // Background pixels: collect their deltas (signed, range [-T, T])
    // Anomaly pixels:    collect (index, value) pairs
    struct Anomaly { uint32_t idx; uint16_t val; };
    std::vector<Anomaly>  anomalies;
    std::vector<int32_t>  deltas;   // only for background pixels, in order

    for (int i = 0; i < NPIX; i++) {
        int32_t delta = (int32_t)px[i] - (int32_t)background;
        if (std::abs(delta) <= (int32_t)T) {
            deltas.push_back(delta);
        } else {
            anomalies.push_back({(uint32_t)i, px[i]});
        }
    }

    // --- Step 4: map deltas to symbols for rANS ---
    // Deltas are in [-T, +T]. Zigzag maps them to [0, 2T].
    // We then take symbol = zigzag(delta) % 256  (low byte)
    // and a separate high-byte stream if T > 127.
    //
    // For simplicity: we always use two streams (hi + lo byte of zigzag value),
    // same as RAIS. Most deltas will have hi=0, making it highly compressible.
    int n_bg = (int)deltas.size();
    std::vector<uint8_t> sym_hi(n_bg), sym_lo(n_bg);
    for (int i = 0; i < n_bg; i++) {
        uint16_t zz = zigzag16(deltas[i]);
        sym_hi[i] = (uint8_t)(zz >> 8);
        sym_lo[i] = (uint8_t)(zz & 0xFF);
    }

    // --- Step 5: build rANS frequency tables ---
    uint64_t cnt_hi[256] = {}, cnt_lo[256] = {};
    for (int i = 0; i < n_bg; i++) { cnt_hi[sym_hi[i]]++; cnt_lo[sym_lo[i]]++; }

    RansTable tab_hi, tab_lo;
    normalise_freqs(cnt_hi, tab_hi.freq, 256);
    normalise_freqs(cnt_lo, tab_lo.freq, 256);
    tab_hi.build_cumul(); tab_hi.build_lookup();
    tab_lo.build_cumul(); tab_lo.build_lookup();

    // --- Step 6: rANS encode (reverse order) ---
    RansEncoder enc_hi, enc_lo;
    for (int i = n_bg - 1; i >= 0; i--) enc_hi.encode(sym_hi[i], tab_hi);
    for (int i = n_bg - 1; i >= 0; i--) enc_lo.encode(sym_lo[i], tab_lo);
    enc_hi.flush();
    enc_lo.flush();

    // --- Step 7: write compressed file ---
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output %s\n", argv[2]); return 1; }

    // Header
    fwrite(MAGIC, 1, 4, fout);
    write_u16le(fout, VERSION);
    write_u32le(fout, (uint32_t)WIDTH);
    write_u32le(fout, (uint32_t)HEIGHT);
    write_u16le(fout, background);
    write_u16le(fout, (uint16_t)T);
    write_u32le(fout, (uint32_t)anomalies.size());

    // Anomaly list: each is 4 bytes (index) + 2 bytes (value) = 6 bytes
    for (const auto& a : anomalies) {
        write_u32le(fout, a.idx);
        write_u16le(fout, a.val);
    }

    // rANS stream 0 (hi bytes)
    write_u32le(fout, 256);
    for (int i = 0; i < 256; i++) write_u32le(fout, tab_hi.freq[i]);
    write_u64le(fout, (uint64_t)enc_hi.buf.size());
    fwrite(enc_hi.buf.data(), 1, enc_hi.buf.size(), fout);

    // rANS stream 1 (lo bytes)
    write_u32le(fout, 256);
    for (int i = 0; i < 256; i++) write_u32le(fout, tab_lo.freq[i]);
    write_u64le(fout, (uint64_t)enc_lo.buf.size());
    fwrite(enc_lo.buf.data(), 1, enc_lo.buf.size(), fout);

    fclose(fout);

    // --- Report ---
    long in_bytes  = NPIX * 2;
    long out_bytes = (long)(4+2+4+4+2+2+4
                   + anomalies.size()*6
                   + (4+256*4+8+enc_hi.buf.size())
                   + (4+256*4+8+enc_lo.buf.size()));
    fprintf(stderr, "background=%u  T=%u  anomalies=%zu (%.2f%%)\n",
            background, T, anomalies.size(),
            anomalies.size() * 100.0 / NPIX);
    fprintf(stderr, "compressed: %ld -> %ld bytes  (%.4f bits/byte)\n",
            in_bytes, out_bytes, out_bytes * 8.0 / in_bytes);
    return 0;
}
