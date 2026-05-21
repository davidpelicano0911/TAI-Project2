//
// HELIX v7 decompressor — exact mirror of compress.cpp
//

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <algorithm>

// ---------------------------------------------------------------------------
// FSE / tANS  (inlined)
// ---------------------------------------------------------------------------

static constexpr uint32_t FSE_SCALE_BITS = 16;
static constexpr uint32_t FSE_SCALE      = 1u << FSE_SCALE_BITS;

struct FseDecodeEntry { uint8_t sym; uint8_t nb_bits; uint16_t base; };

struct FseTable {
    uint32_t       freq[256]  = {};
    FseDecodeEntry dec[FSE_SCALE];
    uint16_t       enc_flat[FSE_SCALE];
    uint32_t       enc_offset[256];
    int            nb_base[256];

    void build() {
        const uint32_t step = (FSE_SCALE >> 1) + (FSE_SCALE >> 3) + 3;
        uint8_t spread[FSE_SCALE];
        uint32_t pos = 0;
        for (int s = 0; s < 256; s++)
            for (uint32_t n = 0; n < freq[s]; n++) {
                spread[pos] = (uint8_t)s;
                pos = (pos + step) & (FSE_SCALE - 1);
            }
        uint32_t cumul = 0;
        for (int s = 0; s < 256; s++) {
            enc_offset[s] = cumul;
            cumul += freq[s];
            nb_base[s] = (freq[s] > 0)
                ? std::max(0, (int)(__builtin_clz(freq[s]) + FSE_SCALE_BITS - 32))
                : 0;
        }
        uint32_t next[256];
        for (int s = 0; s < 256; s++) next[s] = freq[s];
        for (uint32_t state = 0; state < FSE_SCALE; state++) {
            uint8_t  sym  = spread[state];
            uint32_t x    = next[sym]++;
            uint8_t  nb   = (uint8_t)(FSE_SCALE_BITS - (31u - __builtin_clz(x)));
            uint32_t base = (x << nb) - FSE_SCALE;
            dec[state]    = {sym, nb, (uint16_t)base};
            enc_flat[enc_offset[sym] + (x - freq[sym])] = (uint16_t)state;
        }
    }
};

struct FseBitReader {
    const uint8_t* ptr;
    const uint8_t* end;
    uint64_t bitbuf = 0;
    int      bitcnt = 0;
    FseBitReader(const uint8_t* b, const uint8_t* e) : ptr(b), end(e) {}
    uint32_t get(int n) {
        while (bitcnt < n && ptr < end) {
            bitbuf |= (uint64_t)(*ptr++) << bitcnt;
            bitcnt += 8;
        }
        uint32_t v = (uint32_t)(bitbuf & ((1u << n) - 1u));
        bitbuf >>= n; bitcnt -= n;
        return v;
    }
};

struct FseDecoder {
    uint16_t    state;
    FseBitReader br;
    FseDecoder(const uint8_t* bits, const uint8_t* end, uint16_t s)
        : state(s), br(bits, end) {}
    uint8_t next(const FseTable& t) {
        const FseDecodeEntry& e = t.dec[state];
        uint8_t sym = e.sym;
        state = (uint16_t)(e.base + br.get(e.nb_bits));
        return sym;
    }
};

static inline uint32_t pget_u32le(const uint8_t*& p) {
    uint32_t v = (uint32_t)p[0] | ((uint32_t)p[1] << 8)
               | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
    p += 4; return v;
}
static inline uint64_t pget_u64le(const uint8_t*& p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= (uint64_t)p[i] << (8*i);
    p += 8; return v;
}

static void decode_stream(const uint8_t*& ptr, std::vector<uint8_t>& out, int n) {
    uint8_t flag = *ptr++;
    if (flag == 0x01) { out.assign(n, *ptr++); return; }
    FseTable tab;
    if (flag == 0) {
        for (int i = 0; i < 256; i++) tab.freq[i] = pget_u32le(ptr);
    } else {
        for (int k = 0; k < (int)flag; k++) {
            uint8_t sym = *ptr++;
            tab.freq[sym] = pget_u32le(ptr);
        }
    }
    tab.build();
    uint64_t nbytes = pget_u64le(ptr);
    const uint8_t* bits = ptr; ptr += nbytes;
    uint16_t init = (uint16_t)(bits[0] | ((uint16_t)bits[1] << 8));
    FseDecoder dec(bits + 2, bits + nbytes, init);
    out.resize(n);
    for (int i = 0; i < n; i++) out[i] = dec.next(tab);
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

static constexpr uint8_t  MAGIC[4] = {'H','E','L','X'};
static constexpr int      N_CTX    = 27;

// ---------------------------------------------------------------------------
// Predictors + context (mirror of encoder)
// ---------------------------------------------------------------------------

static inline int32_t med_predict(int32_t a, int32_t b, int32_t c) {
    int32_t lo = std::min(a, b);
    int32_t hi = std::max(a, b);
    if      (c >= hi) return lo;
    else if (c <= lo) return hi;
    else              return a + b - c;
}

static inline int32_t avg_predict(int32_t a, int32_t b, int32_t c,
                                   int x, int y) {
    if (y == 0 && x == 0) return 0;
    if (y == 0)            return a;
    if (x == 0)            return b;
    return (a + b + c + 1) / 3;
}

static inline int quantize_grad3(int32_t g, int32_t T1) {
    if (g < -T1) return -1;
    if (g <=  T1) return  0;
    return  1;
}

static inline int ctx_index3(int q1, int q2, int q3) {
    return (q1 + 1) * 9 + (q2 + 1) * 3 + (q3 + 1);
}

static inline int32_t unmap_signed(uint32_t u) {
    return (u & 1) ? -(int32_t)((u + 1) >> 1) : (int32_t)(u >> 1);
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
    uint32_t v = 0;
    for (int i = 3; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 3) { fprintf(stderr, "Usage: %s <input.hlx> <output>\n", argv[0]); return 1; }

    FILE* fin = fopen(argv[1], "rb");
    if (!fin) { fprintf(stderr, "Cannot open: %s\n", argv[1]); return 1; }

    uint8_t magic[4]; fread(magic, 1, 4, fin);
    if (memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n"); fclose(fin); return 1;
    }
    uint16_t version = read_u16le(fin);
    if (version != 7) {
        fprintf(stderr, "Unknown version %u (expected 7)\n", version);
        fclose(fin); return 1;
    }
    uint32_t width  = read_u32le(fin);
    uint32_t height = read_u32le(fin);
    (void)read_u16le(fin); // flags
    uint16_t gmean  = read_u16le(fin);
    uint32_t T1     = read_u32le(fin);
    int npix = (int)(width * height);

    uint32_t bs   = read_u32le(fin);
    uint32_t nbx  = read_u32le(fin);
    uint32_t nby  = read_u32le(fin);
    uint32_t nblk = nbx * nby;

    std::vector<uint8_t> blk_mode(nblk);
    fread(blk_mode.data(), 1, nblk, fin);

    std::vector<std::array<float, 6>> blk_w(nblk);
    for (uint32_t b = 0; b < nblk; b++) {
        for (int k = 0; k < 6; k++) {
            uint32_t bits = read_u32le(fin);
            memcpy(&blk_w[b][k], &bits, 4);
        }
    }

    uint32_t nctx = read_u32le(fin);

    long cur_pos = ftell(fin);
    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    fseek(fin, cur_pos, SEEK_SET);
    std::vector<uint8_t> filebuf(fsize - cur_pos);
    fread(filebuf.data(), 1, filebuf.size(), fin);
    fclose(fin);

    const uint8_t* ptr = filebuf.data();

    // Decode per-context FSE streams
    std::vector<std::vector<uint8_t>> ctx_lo(nctx), ctx_hi(nctx), ctx_top(nctx);
    for (uint32_t q = 0; q < nctx; q++) {
        uint32_t raw = (uint32_t)ptr[0] | ((uint32_t)ptr[1] << 8)
                     | ((uint32_t)ptr[2] << 16) | ((uint32_t)ptr[3] << 24);
        ptr += 4;
        bool lo_split = (raw & 0x80000000u) != 0;
        int  n        = (int)(raw & 0x7FFFFFFFu);
        if (n > 0) {
            decode_stream(ptr, ctx_hi[q],  n);
            decode_stream(ptr, ctx_top[q], n);
            if (lo_split) {
                int n_zero = 0;
                for (int i = 0; i < n; i++)
                    if (ctx_hi[q][i] == 0 && ctx_top[q][i] == 0) n_zero++;
                std::vector<uint8_t> lo_zero, lo_nonzero;
                decode_stream(ptr, lo_zero,    n_zero);
                decode_stream(ptr, lo_nonzero, n - n_zero);
                ctx_lo[q].resize(n);
                int iz = 0, in_ = 0;
                for (int i = 0; i < n; i++) {
                    if (ctx_hi[q][i] == 0 && ctx_top[q][i] == 0)
                        ctx_lo[q][i] = lo_zero[iz++];
                    else
                        ctx_lo[q][i] = lo_nonzero[in_++];
                }
            } else {
                decode_stream(ptr, ctx_lo[q], n);
            }
        }
    }

    // Reconstruct image
    std::vector<int> ctx_idx(nctx, 0);
    std::vector<int32_t> img(npix);

    for (uint32_t y = 0; y < height; y++) {
        for (uint32_t x = 0; x < width; x++) {
            int32_t a = (x > 0)                  ? img[y * width + (x - 1)]        : 0;
            int32_t b = (y > 0)                  ? img[(y - 1) * width + x]        : 0;
            int32_t c = (x > 0 && y > 0)         ? img[(y - 1) * width + (x - 1)] : 0;
            int32_t d = (y > 0 && x + 1 < width) ? img[(y - 1) * width + (x + 1)] : b;

            int32_t g1 = d - b, g2 = b - c, g3 = c - a;
            int q1 = quantize_grad3(g1, (int32_t)T1);
            int q2 = quantize_grad3(g2, (int32_t)T1);
            int q3 = quantize_grad3(g3, (int32_t)T1);
            int sign = 1;
            if (q1 < 0 || (q1 == 0 && q2 < 0) || (q1 == 0 && q2 == 0 && q3 < 0)) {
                sign = -1;
                q1 = -q1; q2 = -q2; q3 = -q3;
            }
            int Q = ctx_index3(q1, q2, q3);

            uint32_t bidx_cur = (y / bs) * nbx + (x / bs);
            int mode = blk_mode[bidx_cur];
            int32_t pred;
            if (mode == 0) {
                if (y == 0 && x == 0) pred = 0;
                else if (y == 0)      pred = a;
                else if (x == 0)      pred = b;
                else                  pred = med_predict(a, b, c);
            } else if (mode == 1) {
                pred = avg_predict(a, b, c, (int)x, (int)y);
            } else if (mode == 2) {
                pred = (int32_t)gmean;
            } else {
                if (y < 2 || x < 2) {
                    if (y == 0 && x == 0) pred = 0;
                    else if (y == 0)      pred = a;
                    else if (x == 0)      pred = b;
                    else                  pred = med_predict(a, b, c);
                } else {
                    const auto& w = blk_w[bidx_cur];
                    float p = w[0]*img[y*width+(x-1)]
                            + w[1]*img[(y-1)*width+x]
                            + w[2]*img[(y-1)*width+(x-1)]
                            + w[3]*img[y*width+(x-2)]
                            + w[4]*img[(y-2)*width+x]
                            + w[5];
                    pred = std::max(0, std::min(65535, (int32_t)(p + 0.5f)));
                }
            }

            int idx = ctx_idx[Q]++;
            uint8_t lo  = ctx_lo[Q][idx];
            uint8_t hi  = ctx_hi[Q][idx];
            uint8_t top = ctx_top[Q][idx];
            uint32_t u  = ((uint32_t)top << 16) | ((uint32_t)hi << 8) | lo;

            int32_t e_c   = unmap_signed(u);
            int32_t e_raw = sign * e_c;
            img[y * width + x] = pred + e_raw;
        }
    }

    // Write big-endian output
    FILE* fout = fopen(argv[2], "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", argv[2]); return 1; }

    for (uint32_t y = 0; y < height; y++) {
        for (uint32_t x = 0; x < width; x++) {
            int32_t v = std::max(0, std::min(65535, img[y * width + x]));
            uint16_t px = (uint16_t)v;
            uint8_t b2[2] = {(uint8_t)(px >> 8), (uint8_t)(px & 0xFF)};
            fwrite(b2, 1, 2, fout);
        }
    }
    fclose(fout);
    return 0;
}
