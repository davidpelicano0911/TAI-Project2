// HAIS — Hybrid Astronomical Image decompressor
// Usage: ./decompress <input.hais> <output>

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>

static constexpr uint8_t MAGIC[4] = {'H','A','I','S'};

static constexpr uint32_t RANS_SCALE_BITS = 16;
static constexpr uint32_t RANS_SCALE      = 1u << RANS_SCALE_BITS;
static constexpr uint32_t RANS_L          = 1u << 23;

// ---------------------------------------------------------------------------
// rANS decoder
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
        state = tab.freq[sym] * (state / RANS_SCALE) + slot - tab.cumul[sym];
        while (state < RANS_L) state = (state << 8) | (*ptr++);
        return sym;
    }
};

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static uint8_t  read_u8   (FILE* f) { uint8_t v; fread(&v, 1, 1, f); return v; }
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
// Zagzig: inverso do zigzag (uint16 → int16)
// ---------------------------------------------------------------------------

static inline int16_t zagzig(uint16_t u) {
    return (u & 1) ? -(int16_t)((u + 1) / 2) : (int16_t)(u / 2);
}

// ---------------------------------------------------------------------------
// Preditor MED — igual ao compressor
// ---------------------------------------------------------------------------

static inline uint16_t med_predict(const std::vector<uint16_t>& block,
                                    int x, int y, int bs) {
    if (y == 0 && x == 0) return 0;
    if (y == 0)           return block[x - 1];
    if (x == 0)           return block[(y - 1) * bs + x];
    uint16_t A = block[y       * bs + (x - 1)];
    uint16_t B = block[(y - 1) * bs +  x     ];
    uint16_t C = block[(y - 1) * bs + (x - 1)];
    int pred = (int)A + (int)B - (int)C;
    int lo   = std::min((int)A, (int)B);
    int hi   = std::max((int)A, (int)B);
    if (pred < lo) pred = lo;
    if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

// ---------------------------------------------------------------------------
// Ler e descodificar um stream rANS do ficheiro
// ---------------------------------------------------------------------------

static void decode_stream(FILE* fin, std::vector<uint8_t>& out, int n) {
    RansTable tab;
    memset(tab.freq, 0, sizeof(tab.freq));

    // Ler tabela: 0=densa, 1-255=esparsa
    uint8_t flag = read_u8(fin);
    if (flag == 0) {
        for (int i = 0; i < 256; i++) tab.freq[i] = read_u32le(fin);
    } else {
        int nnz = (int)flag;
        for (int k = 0; k < nnz; k++) {
            uint8_t  sym  = read_u8(fin);
            uint32_t freq = read_u32le(fin);
            tab.freq[sym] = freq;
        }
    }

    tab.build_cumul();
    tab.build_lookup();

    uint64_t nbytes = read_u64le(fin);
    std::vector<uint8_t> bits(nbytes);
    fread(bits.data(), 1, nbytes, fin);

    out.resize(n);
    RansDecoder dec;
    dec.init(bits.data());
    for (int i = 0; i < n; i++) out[i] = dec.decode(tab);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input.hais> <output>\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }

    uint8_t magic[4]; fread(magic, 1, 4, fin);
    if (memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n"); fclose(fin); return 1;
    }
    read_u16le(fin);
    uint32_t width  = read_u32le(fin);
    uint32_t height = read_u32le(fin);
    uint32_t bs     = read_u32le(fin);

    int blocks_x   = (int)width  / (int)bs;
    int blocks_y   = (int)height / (int)bs;
    int block_npix = (int)(bs * bs);

    std::vector<uint16_t> image(width * height);
    std::vector<uint16_t> block(block_npix);

    for (int by = 0; by < blocks_y; by++) {
        for (int bx = 0; bx < blocks_x; bx++) {
            int mode = (int)read_u8(fin);

            // Descodificar streams hi8 e lo8
            std::vector<uint8_t> hi8, lo8;
            decode_stream(fin, hi8, block_npix);
            decode_stream(fin, lo8, block_npix);

            // Reconstruir símbolos uint16
            // Inverter modelo para obter pixels
            for (int y = 0; y < (int)bs; y++) {
                for (int x = 0; x < (int)bs; x++) {
                    int idx = y * bs + x;
                    uint16_t sym = ((uint16_t)hi8[idx] << 8) | lo8[idx];
                    uint16_t pixel;

                    if (mode == 0) {
                        pixel = sym;
                    } else if (mode == 1) {
                        uint16_t pred;
                        if      (x > 0) pred = block[y * bs + (x - 1)];
                        else if (y > 0) pred = block[(y - 1) * bs + x];
                        else            pred = 0;
                        pixel = (uint16_t)(pred + zagzig(sym));
                    } else if (mode == 2) {
                        uint16_t pred = med_predict(block, x, y, (int)bs);
                        pixel = (uint16_t)(pred + zagzig(sym));
                    } else {
                        // mode 3: smooth causal predictor (Quintas-Torra et al. 2026)
                        uint16_t pred;
                        if (y == 0 && x == 0)   pred = 0;
                        else if (y == 0)         pred = block[x - 1];
                        else if (x == 0)         pred = block[(y - 1) * bs + x];
                        else {
                            uint16_t A = block[y       * bs + (x - 1)];
                            uint16_t B = block[(y - 1) * bs +  x     ];
                            uint16_t C = block[(y - 1) * bs + (x - 1)];
                            pred = (uint16_t)(((int)A + (int)B + (int)C + 1) / 3);
                        }
                        pixel = (uint16_t)(pred + zagzig(sym));
                    }

                    block[idx] = pixel;
                }
            }

            // Copiar bloco para imagem
            for (int y = 0; y < (int)bs; y++) {
                int gy = by * bs + y;
                for (int x = 0; x < (int)bs; x++)
                    image[gy * width + bx * bs + x] = block[y * bs + x];
            }
        }
    }
    fclose(fin);

    // Escrever RAW output (big-endian uint16)
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open %s\n", argv[2]); return 1; }
    for (uint16_t px : image) {
        uint8_t b[2] = {(uint8_t)(px >> 8), (uint8_t)(px & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}
