// PRISM v3 decompressor — adaptive range coder edition

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>

static constexpr uint8_t MAGIC[4] = {'P','R','M','3'};
static constexpr int     N_LO_CTX = 8;
static constexpr int     HEADER_SIZE = 16;

// ---------------------------------------------------------------------------
// GAP predictor
// ---------------------------------------------------------------------------

static inline uint16_t gap_predict(int W, int N, int NW) {
    int dh = std::abs(W - NW);
    int dv = std::abs(N - NW);
    if (dh > dv * 2) return (uint16_t)N;
    if (dv > dh * 2) return (uint16_t)W;
    int wW = dv + 1, wN = dh + 1;
    int pred = (wW * W + wN * N) / (wW + wN);
    int lo = std::min(W, N), hi = std::max(W, N);
    if (pred < lo) pred = lo;
    if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

// ---------------------------------------------------------------------------
// Modular zigzag decode
// ---------------------------------------------------------------------------

static inline uint16_t zigzag_dec(uint16_t z) {
    return (z & 1u) ? (uint16_t)(65536u - (z + 1u) / 2u) : (uint16_t)(z / 2u);
}

// ---------------------------------------------------------------------------
// Adaptive frequency model
// ---------------------------------------------------------------------------

static constexpr uint32_t MAX_TOTAL = 1u << 16;

struct AdaptModel {
    uint32_t freq[256];
    uint32_t cumul[257];
    uint32_t total;

    void init() {
        for (int i = 0; i < 256; i++) freq[i] = 1;
        total = 256;
        build_cumul();
    }

    void build_cumul() {
        cumul[0] = 0;
        for (int i = 0; i < 256; i++) cumul[i+1] = cumul[i] + freq[i];
    }

    void update(uint8_t sym) {
        freq[sym]++;
        total++;
        if (total >= MAX_TOTAL) {
            total = 0;
            for (int i = 0; i < 256; i++) {
                freq[i] = (freq[i] + 1) >> 1;
                total += freq[i];
            }
            build_cumul();
            return;
        }
        for (int i = (int)sym + 1; i <= 256; i++) cumul[i]++;
    }

    uint8_t find(uint32_t slot) const {
        int lo = 0, hi = 255;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (cumul[mid+1] <= slot) lo = mid + 1;
            else                      hi = mid;
        }
        return (uint8_t)lo;
    }
};

// ---------------------------------------------------------------------------
// Range decoder
// ---------------------------------------------------------------------------

struct RangeDecoder {
    const uint8_t* ptr;
    uint32_t range = 0xFFFFFFFFu;
    uint32_t code  = 0;

    void init(const uint8_t* data) {
        ptr = data;
        for (int i = 0; i < 5; i++) code = (code << 8) | (*ptr++);
    }

    uint8_t decode(AdaptModel& m) {
        uint32_t r    = range / m.total;
        uint32_t slot = std::min(code / r, m.total - 1);
        uint8_t  sym  = m.find(slot);
        code  -= m.cumul[sym] * r;
        range  = (sym < 255) ? (m.cumul[sym+1] - m.cumul[sym]) * r
                             : range - m.cumul[sym] * r;
        while (range < (1u << 24)) { code = (code << 8) | (*ptr++); range <<= 8; }
        m.update(sym);
        return sym;
    }
};

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static bool read_exact(FILE* f, void* data, size_t n) {
    return fread(data, 1, n, f) == n;
}

static bool read_u32le(FILE* f, uint32_t& v) {
    uint8_t b[4];
    if (!read_exact(f, b, 4)) return false;
    v = 0;
    for (int i = 3; i >= 0; i--) v = (v << 8) | b[i];
    return true;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) { fprintf(stderr, "Usage: %s <input.prism> <output>\n", argv[0]); return 1; }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open: %s\n", argv[1]); return 1; }

    uint8_t magic[4];
    if (!read_exact(fin, magic, 4)) { fprintf(stderr, "Truncated header\n"); fclose(fin); return 1; }
    if (memcmp(magic, MAGIC, 4) != 0) { fprintf(stderr, "Bad magic\n"); fclose(fin); return 1; }

    uint32_t width = 0, height = 0, n_lo_ctx = 0;
    if (!read_u32le(fin, width) || !read_u32le(fin, height) || !read_u32le(fin, n_lo_ctx)) {
        fprintf(stderr, "Truncated header\n");
        fclose(fin);
        return 1;
    }
    if (n_lo_ctx != N_LO_CTX) { fprintf(stderr, "Unexpected n_lo_ctx\n"); fclose(fin); return 1; }

    // Read remaining bytes as the range-coded stream
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    long stream_sz = fsize - HEADER_SIZE;
    if (stream_sz < 5) { fprintf(stderr, "Truncated range stream\n"); fclose(fin); return 1; }
    fseek(fin, HEADER_SIZE, SEEK_SET);
    std::vector<uint8_t> stream(stream_sz + 16, 0);  // pad for overread safety
    if (!read_exact(fin, stream.data(), (size_t)stream_sz)) {
        fprintf(stderr, "Range stream read failed\n");
        fclose(fin);
        return 1;
    }
    fclose(fin);

    int npix = (int)(width * height);

    AdaptModel m_hi;
    m_hi.init();
    AdaptModel m_lo[N_LO_CTX];
    for (int c = 0; c < N_LO_CTX; c++) m_lo[c].init();

    RangeDecoder dec;
    dec.init(stream.data());

    std::vector<uint16_t> image(npix);
    std::vector<int16_t>  raw_res(npix, 0);

    for (int gy = 0; gy < (int)height; gy++) {
        for (int gx = 0; gx < (int)width; gx++) {
            int W  = (gx > 0)           ? image[gy*(int)width+gx-1]     : 0;
            int N  = (gy > 0)           ? image[(gy-1)*(int)width+gx]   : W;
            int NW = (gy > 0 && gx > 0) ? image[(gy-1)*(int)width+gx-1] : W;

            uint16_t pred = (gy == 0 && gx == 0) ? 0u : gap_predict(W, N, NW);

            int mW = (gx > 0) ? std::abs((int)raw_res[gy*(int)width+gx-1]) : 0;
            int mN = (gy > 0) ? std::abs((int)raw_res[(gy-1)*(int)width+gx]) : 0;
            int mg = mW + mN;
            int ctx;
            if      (mg ==  0) ctx = 0;
            else if (mg <=  2) ctx = 1;
            else if (mg <=  8) ctx = 2;
            else if (mg <= 32) ctx = 3;
            else if (mg <= 64) ctx = 4;
            else if (mg <= 128) ctx = 5;
            else if (mg <= 512) ctx = 6;
            else                ctx = 7;

            uint8_t  hi  = dec.decode(m_hi);
            uint8_t  lo  = dec.decode(m_lo[ctx]);
            uint16_t zz  = ((uint16_t)hi << 8) | lo;
            uint16_t u   = zigzag_dec(zz);
            int16_t  res = (u <= 32767u) ? (int16_t)u : (int16_t)((int)u - 65536);
            raw_res[gy*(int)width+gx] = res;
            image[gy*(int)width+gx]   = (uint16_t)(pred + u);
        }
    }

    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", argv[2]); return 1; }
    for (int i = 0; i < npix; i++) {
        uint8_t b[2] = { (uint8_t)(image[i] >> 8), (uint8_t)(image[i] & 0xFF) };
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
