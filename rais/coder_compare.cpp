//
// Experimental RAIS entropy-coder comparison tool.
//
// Usage:
//   ./coder_compare c <rans|huffman|arith> <input> <output>
//   ./coder_compare d <input> <output>
//
// The predictor and residual representation match rais/compress.cpp:
//   raw uint16 pixels -> MED prediction -> zigzag residual -> high/low byte streams.
//

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <vector>

static constexpr int WIDTH_DEFAULT  = 1500;
static constexpr int HEIGHT_DEFAULT = 1500;
static constexpr uint8_t MAGIC[4]   = {'R','A','I','C'};
static constexpr uint16_t VERSION   = 1;

enum class Coder : uint8_t {
    Rans    = 1,
    Huffman = 2,
    Arith   = 3,
};

static constexpr uint32_t ENT_SCALE_BITS = 16;
static constexpr uint32_t ENT_SCALE      = 1u << ENT_SCALE_BITS;
static constexpr uint32_t RANS_L         = 1u << 23;

// ---------------------------------------------------------------------------
// Predictor and residual mapping
// ---------------------------------------------------------------------------

static inline uint16_t predict(const uint16_t* row, const uint16_t* prev_row,
                               int x, int /*width*/) {
    if (prev_row == nullptr && x == 0) return 0;
    if (prev_row == nullptr) return row[x - 1];
    if (x == 0) return prev_row[x];

    uint16_t a = row[x - 1];
    uint16_t b = prev_row[x];
    uint16_t c = prev_row[x - 1];

    int pred = (int)a + (int)b - (int)c;
    int lo = std::min((int)a, (int)b);
    int hi = std::max((int)a, (int)b);
    if (pred < lo) pred = lo;
    if (pred > hi) pred = hi;
    return (uint16_t)pred;
}

static inline uint16_t zigzag_enc(int16_t v) {
    return (uint16_t)((v >= 0) ? (v * 2) : ((-v) * 2 - 1));
}

static inline int16_t zigzag_dec(uint16_t v) {
    return (v & 1) ? -(int16_t)((v + 1) / 2) : (int16_t)(v / 2);
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

static void write_u8(FILE* f, uint8_t v) {
    fwrite(&v, 1, 1, f);
}

static void write_u16le(FILE* f, uint16_t v) {
    uint8_t b[2] = {(uint8_t)(v & 0xFF), (uint8_t)(v >> 8)};
    fwrite(b, 1, 2, f);
}

static void write_u32le(FILE* f, uint32_t v) {
    uint8_t b[4];
    for (int i = 0; i < 4; i++) { b[i] = (uint8_t)(v & 0xFF); v >>= 8; }
    fwrite(b, 1, 4, f);
}

static void write_u64le(FILE* f, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) { b[i] = (uint8_t)(v & 0xFF); v >>= 8; }
    fwrite(b, 1, 8, f);
}

static bool read_exact(FILE* f, void* data, size_t n) {
    return fread(data, 1, n, f) == n;
}

static uint8_t read_u8(FILE* f) {
    uint8_t v = 0;
    read_exact(f, &v, 1);
    return v;
}

static uint16_t read_u16le(FILE* f) {
    uint8_t b[2]; read_exact(f, b, 2);
    return (uint16_t)(b[0] | ((uint16_t)b[1] << 8));
}

static uint32_t read_u32le(FILE* f) {
    uint8_t b[4]; read_exact(f, b, 4);
    uint32_t v = 0;
    for (int i = 3; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}

static uint64_t read_u64le(FILE* f) {
    uint8_t b[8]; read_exact(f, b, 8);
    uint64_t v = 0;
    for (int i = 7; i >= 0; i--) v = (v << 8) | b[i];
    return v;
}

// ---------------------------------------------------------------------------
// Bit I/O
// ---------------------------------------------------------------------------

struct BitWriter {
    std::vector<uint8_t> out;
    uint8_t cur = 0;
    int fill = 0;

    void put_bit(int bit) {
        cur = (uint8_t)((cur << 1) | (bit & 1));
        fill++;
        if (fill == 8) {
            out.push_back(cur);
            cur = 0;
            fill = 0;
        }
    }

    void put_bits(uint32_t code, int nbits) {
        for (int i = nbits - 1; i >= 0; --i) put_bit((code >> i) & 1u);
    }

    void flush_zero() {
        if (fill) {
            cur <<= (8 - fill);
            out.push_back(cur);
            cur = 0;
            fill = 0;
        }
    }
};

struct BitReader {
    const uint8_t* data = nullptr;
    size_t size = 0;
    size_t pos = 0;
    uint8_t cur = 0;
    int left = 0;

    BitReader() = default;
    BitReader(const uint8_t* data_, size_t size_) : data(data_), size(size_) {}

    int get_bit() {
        if (left == 0) {
            if (pos >= size) return 0;
            cur = data[pos++];
            left = 8;
        }
        int bit = (cur >> 7) & 1;
        cur <<= 1;
        left--;
        return bit;
    }
};

// ---------------------------------------------------------------------------
// Shared static model
// ---------------------------------------------------------------------------

struct EntTable {
    uint32_t freq[256] = {};
    uint32_t cumul[257] = {};
    uint8_t sym_of_slot[ENT_SCALE] = {};

    void build_cumul() {
        cumul[0] = 0;
        for (int i = 0; i < 256; i++) cumul[i + 1] = cumul[i] + freq[i];
    }

    void build_lookup() {
        for (int s = 0; s < 256; s++) {
            for (uint32_t j = cumul[s]; j < cumul[s + 1]; j++) {
                sym_of_slot[j] = (uint8_t)s;
            }
        }
    }
};

static void normalise_freqs(const uint64_t* counts, uint32_t* freq) {
    uint64_t total = 0;
    for (int i = 0; i < 256; i++) total += counts[i];
    if (total == 0) return;

    uint32_t assigned = 0;
    for (int i = 0; i < 256; i++) {
        if (counts[i] == 0) {
            freq[i] = 0;
        } else {
            uint64_t scaled = (counts[i] * (uint64_t)ENT_SCALE) / total;
            freq[i] = (uint32_t)std::max<uint64_t>(1, scaled);
            assigned += freq[i];
        }
    }

    while (assigned < ENT_SCALE) {
        int best = 0;
        for (int i = 1; i < 256; i++) {
            if (counts[i] > counts[best]) best = i;
        }
        freq[best]++;
        assigned++;
    }

    while (assigned > ENT_SCALE) {
        int best = -1;
        for (int i = 0; i < 256; i++) {
            if (freq[i] > 1 && (best < 0 || freq[i] > freq[best])) best = i;
        }
        if (best < 0) break;
        freq[best]--;
        assigned--;
    }
}

// ---------------------------------------------------------------------------
// rANS coder
// ---------------------------------------------------------------------------

struct RansEncoder {
    uint32_t state = RANS_L;
    std::vector<uint8_t> buf;

    void encode(uint8_t sym, const EntTable& tab) {
        uint32_t f = tab.freq[sym];
        uint32_t c = tab.cumul[sym];
        uint32_t x_max = ((RANS_L / ENT_SCALE) * 256) * f;
        while (state >= x_max) {
            buf.push_back((uint8_t)(state & 0xFF));
            state >>= 8;
        }
        state = (state / f) * ENT_SCALE + c + (state % f);
    }

    void flush() {
        buf.push_back((uint8_t)(state & 0xFF));
        buf.push_back((uint8_t)((state >> 8) & 0xFF));
        buf.push_back((uint8_t)((state >> 16) & 0xFF));
        buf.push_back((uint8_t)((state >> 24) & 0xFF));
        std::reverse(buf.begin(), buf.end());
    }
};

struct RansDecoder {
    uint32_t state = 0;
    const uint8_t* ptr = nullptr;

    void init(const uint8_t* data) {
        ptr = data;
        state  = (uint32_t)(*ptr++) << 24;
        state |= (uint32_t)(*ptr++) << 16;
        state |= (uint32_t)(*ptr++) << 8;
        state |= (uint32_t)(*ptr++);
    }

    uint8_t decode(const EntTable& tab) {
        uint32_t slot = state % ENT_SCALE;
        uint8_t sym = tab.sym_of_slot[slot];
        state = tab.freq[sym] * (state / ENT_SCALE) + slot - tab.cumul[sym];
        while (state < RANS_L) state = (state << 8) | (*ptr++);
        return sym;
    }
};

static std::vector<uint8_t> encode_rans(const std::vector<uint8_t>& src,
                                        const EntTable& tab) {
    RansEncoder enc;
    for (int i = (int)src.size() - 1; i >= 0; --i) enc.encode(src[i], tab);
    enc.flush();
    return enc.buf;
}

static std::vector<uint8_t> decode_rans(const std::vector<uint8_t>& bits,
                                        const EntTable& tab,
                                        int nsyms) {
    std::vector<uint8_t> out(nsyms);
    RansDecoder dec;
    dec.init(bits.data());
    for (int i = 0; i < nsyms; i++) out[i] = dec.decode(tab);
    return out;
}

// ---------------------------------------------------------------------------
// Canonical Huffman coder
// ---------------------------------------------------------------------------

struct HuffCode {
    uint32_t code = 0;
    uint8_t len = 0;
};

struct HuffNode {
    uint64_t freq = 0;
    int sym = -1;
    int left = -1;
    int right = -1;
};

struct HuffQueueItem {
    uint64_t freq;
    int node;
    bool operator>(const HuffQueueItem& other) const {
        if (freq != other.freq) return freq > other.freq;
        return node > other.node;
    }
};

static void assign_huff_lengths(const std::vector<HuffNode>& nodes,
                                int idx,
                                int depth,
                                uint8_t* lengths) {
    const HuffNode& n = nodes[idx];
    if (n.sym >= 0) {
        lengths[n.sym] = (uint8_t)std::max(1, depth);
        return;
    }
    assign_huff_lengths(nodes, n.left, depth + 1, lengths);
    assign_huff_lengths(nodes, n.right, depth + 1, lengths);
}

static void build_huffman_lengths(const uint64_t* counts, uint8_t* lengths) {
    std::fill(lengths, lengths + 256, 0);
    std::vector<HuffNode> nodes;
    std::priority_queue<HuffQueueItem, std::vector<HuffQueueItem>,
                        std::greater<HuffQueueItem>> pq;

    for (int s = 0; s < 256; s++) {
        if (counts[s] == 0) continue;
        int idx = (int)nodes.size();
        nodes.push_back({counts[s], s, -1, -1});
        pq.push({counts[s], idx});
    }

    if (pq.empty()) return;
    if (pq.size() == 1) {
        lengths[nodes[pq.top().node].sym] = 1;
        return;
    }

    while (pq.size() > 1) {
        auto a = pq.top(); pq.pop();
        auto b = pq.top(); pq.pop();
        int idx = (int)nodes.size();
        nodes.push_back({a.freq + b.freq, -1, a.node, b.node});
        pq.push({a.freq + b.freq, idx});
    }

    assign_huff_lengths(nodes, pq.top().node, 0, lengths);
}

static std::array<HuffCode, 256> build_canonical_codes(const uint8_t* lengths) {
    std::vector<std::pair<int, int>> symbols;
    for (int s = 0; s < 256; s++) {
        if (lengths[s]) symbols.push_back({lengths[s], s});
    }
    std::sort(symbols.begin(), symbols.end());

    std::array<HuffCode, 256> codes{};
    uint32_t code = 0;
    int prev_len = 0;
    for (auto [len, sym] : symbols) {
        assert(len <= 32);
        code <<= (len - prev_len);
        codes[sym] = {code, (uint8_t)len};
        code++;
        prev_len = len;
    }
    return codes;
}

static std::vector<uint8_t> encode_huffman(const std::vector<uint8_t>& src,
                                           const uint8_t* lengths) {
    auto codes = build_canonical_codes(lengths);
    BitWriter bw;
    for (uint8_t sym : src) {
        bw.put_bits(codes[sym].code, codes[sym].len);
    }
    bw.flush_zero();
    return bw.out;
}

struct HuffDecodeNode {
    int child[2] = {-1, -1};
    int sym = -1;
};

static std::vector<uint8_t> decode_huffman(const std::vector<uint8_t>& bits,
                                           const uint8_t* lengths,
                                           int nsyms) {
    auto codes = build_canonical_codes(lengths);
    std::vector<HuffDecodeNode> tree(1);

    for (int s = 0; s < 256; s++) {
        if (!codes[s].len) continue;
        int node = 0;
        for (int bitpos = codes[s].len - 1; bitpos >= 0; --bitpos) {
            int bit = (codes[s].code >> bitpos) & 1u;
            if (tree[node].child[bit] < 0) {
                tree[node].child[bit] = (int)tree.size();
                tree.push_back({});
            }
            node = tree[node].child[bit];
        }
        tree[node].sym = s;
    }

    BitReader br(bits.data(), bits.size());
    std::vector<uint8_t> out(nsyms);
    for (int i = 0; i < nsyms; i++) {
        int node = 0;
        while (tree[node].sym < 0) {
            int bit = br.get_bit();
            node = tree[node].child[bit];
            if (node < 0) {
                fprintf(stderr, "Invalid Huffman bitstream\n");
                std::exit(1);
            }
        }
        out[i] = (uint8_t)tree[node].sym;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Static arithmetic coder
// ---------------------------------------------------------------------------

static constexpr int AC_BITS = 56;
static constexpr uint64_t AC_TOP = (1ULL << AC_BITS) - 1;
static constexpr uint64_t AC_FIRST_QTR = (AC_TOP + 1) / 4;
static constexpr uint64_t AC_HALF = AC_FIRST_QTR * 2;
static constexpr uint64_t AC_THIRD_QTR = AC_FIRST_QTR * 3;

struct ArithmeticEncoder {
    uint64_t low = 0;
    uint64_t high = AC_TOP;
    uint64_t pending = 0;
    BitWriter bw;

    void output_bit_plus_pending(int bit) {
        bw.put_bit(bit);
        while (pending > 0) {
            bw.put_bit(!bit);
            pending--;
        }
    }

    void encode(uint8_t sym, const EntTable& tab) {
        uint64_t range = high - low + 1;
        uint32_t c0 = tab.cumul[sym];
        uint32_t c1 = tab.cumul[sym + 1];
        high = low + (uint64_t)(((__uint128_t)range * c1) / ENT_SCALE) - 1;
        low  = low + (uint64_t)(((__uint128_t)range * c0) / ENT_SCALE);

        for (;;) {
            if (high < AC_HALF) {
                output_bit_plus_pending(0);
            } else if (low >= AC_HALF) {
                output_bit_plus_pending(1);
                low -= AC_HALF;
                high -= AC_HALF;
            } else if (low >= AC_FIRST_QTR && high < AC_THIRD_QTR) {
                pending++;
                low -= AC_FIRST_QTR;
                high -= AC_FIRST_QTR;
            } else {
                break;
            }
            low <<= 1;
            high = (high << 1) | 1u;
        }
    }

    std::vector<uint8_t> finish() {
        pending++;
        if (low < AC_FIRST_QTR) output_bit_plus_pending(0);
        else output_bit_plus_pending(1);
        bw.flush_zero();
        return bw.out;
    }
};

struct ArithmeticDecoder {
    uint64_t low = 0;
    uint64_t high = AC_TOP;
    uint64_t code = 0;
    BitReader br;

    ArithmeticDecoder(const uint8_t* data, size_t size) : br(data, size) {
        for (int i = 0; i < AC_BITS; i++) code = (code << 1) | (uint64_t)br.get_bit();
    }

    uint8_t decode(const EntTable& tab) {
        uint64_t range = high - low + 1;
        uint32_t value = (uint32_t)((((__uint128_t)(code - low + 1) * ENT_SCALE) - 1) / range);
        uint8_t sym = tab.sym_of_slot[value];

        uint32_t c0 = tab.cumul[sym];
        uint32_t c1 = tab.cumul[sym + 1];
        high = low + (uint64_t)(((__uint128_t)range * c1) / ENT_SCALE) - 1;
        low  = low + (uint64_t)(((__uint128_t)range * c0) / ENT_SCALE);

        for (;;) {
            if (high < AC_HALF) {
                // keep interval in lower half
            } else if (low >= AC_HALF) {
                code -= AC_HALF;
                low -= AC_HALF;
                high -= AC_HALF;
            } else if (low >= AC_FIRST_QTR && high < AC_THIRD_QTR) {
                code -= AC_FIRST_QTR;
                low -= AC_FIRST_QTR;
                high -= AC_FIRST_QTR;
            } else {
                break;
            }
            low <<= 1;
            high = (high << 1) | 1u;
            code = (code << 1) | (uint64_t)br.get_bit();
        }
        return sym;
    }
};

static std::vector<uint8_t> encode_arithmetic(const std::vector<uint8_t>& src,
                                              const EntTable& tab) {
    ArithmeticEncoder enc;
    for (uint8_t sym : src) enc.encode(sym, tab);
    return enc.finish();
}

static std::vector<uint8_t> decode_arithmetic(const std::vector<uint8_t>& bits,
                                              const EntTable& tab,
                                              int nsyms) {
    std::vector<uint8_t> out(nsyms);
    ArithmeticDecoder dec(bits.data(), bits.size());
    for (int i = 0; i < nsyms; i++) out[i] = dec.decode(tab);
    return out;
}

// ---------------------------------------------------------------------------
// Image transform
// ---------------------------------------------------------------------------

static bool read_pixels(const char* path, std::vector<uint16_t>& pixels) {
    FILE* f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "Cannot open input: %s\n", path); return false; }
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    rewind(f);
    int npix = WIDTH_DEFAULT * HEIGHT_DEFAULT;
    if (fsize != (long)npix * 2) {
        fprintf(stderr, "Unexpected file size %ld (expected %d)\n", fsize, npix * 2);
        fclose(f);
        return false;
    }
    pixels.resize(npix);
    for (int i = 0; i < npix; i++) {
        uint8_t b[2];
        if (!read_exact(f, b, 2)) {
            fprintf(stderr, "Short read\n");
            fclose(f);
            return false;
        }
        pixels[i] = (uint16_t)(((uint16_t)b[0] << 8) | b[1]);
    }
    fclose(f);
    return true;
}

static void make_streams(const std::vector<uint16_t>& pixels,
                         std::vector<uint8_t>& hi,
                         std::vector<uint8_t>& lo) {
    int width = WIDTH_DEFAULT;
    int npix = (int)pixels.size();
    hi.resize(npix);
    lo.resize(npix);
    for (int i = 0; i < npix; i++) {
        int x = i % width;
        int y = i / width;
        const uint16_t* row = pixels.data() + y * width;
        const uint16_t* prev_row = (y > 0) ? (pixels.data() + (y - 1) * width) : nullptr;
        uint16_t pred = predict(row, prev_row, x, width);
        int16_t res = (int16_t)(pixels[i] - pred);
        uint16_t zz = zigzag_enc(res);
        hi[i] = (uint8_t)(zz >> 8);
        lo[i] = (uint8_t)(zz & 0xFF);
    }
}

static std::vector<uint16_t> reconstruct(const std::vector<uint8_t>& hi,
                                         const std::vector<uint8_t>& lo,
                                         int width,
                                         int height) {
    int npix = width * height;
    std::vector<uint16_t> pixels(npix);
    for (int i = 0; i < npix; i++) {
        int x = i % width;
        int y = i / width;
        const uint16_t* row = pixels.data() + y * width;
        const uint16_t* prev_row = (y > 0) ? (pixels.data() + (y - 1) * width) : nullptr;
        uint16_t pred = predict(row, prev_row, x, width);
        uint16_t zz = ((uint16_t)hi[i] << 8) | lo[i];
        int16_t res = zigzag_dec(zz);
        pixels[i] = (uint16_t)((int)pred + res);
    }
    return pixels;
}

// ---------------------------------------------------------------------------
// File format streams
// ---------------------------------------------------------------------------

struct EncodedStream {
    EntTable ent;
    uint8_t huff_lengths[256] = {};
    std::vector<uint8_t> payload;
};

static EncodedStream encode_stream(const std::vector<uint8_t>& src, Coder coder) {
    uint64_t counts[256] = {};
    for (uint8_t v : src) counts[v]++;

    EncodedStream s;
    if (coder == Coder::Huffman) {
        build_huffman_lengths(counts, s.huff_lengths);
        s.payload = encode_huffman(src, s.huff_lengths);
        return s;
    }

    normalise_freqs(counts, s.ent.freq);
    s.ent.build_cumul();
    s.ent.build_lookup();
    if (coder == Coder::Rans) s.payload = encode_rans(src, s.ent);
    else s.payload = encode_arithmetic(src, s.ent);
    return s;
}

static void write_stream(FILE* f, const EncodedStream& s, Coder coder) {
    if (coder == Coder::Huffman) {
        fwrite(s.huff_lengths, 1, 256, f);
    } else {
        for (int i = 0; i < 256; i++) write_u32le(f, s.ent.freq[i]);
    }
    write_u64le(f, (uint64_t)s.payload.size());
    fwrite(s.payload.data(), 1, s.payload.size(), f);
}

static bool read_stream(FILE* f, EncodedStream& s, Coder coder) {
    if (coder == Coder::Huffman) {
        if (!read_exact(f, s.huff_lengths, 256)) return false;
    } else {
        for (int i = 0; i < 256; i++) s.ent.freq[i] = read_u32le(f);
        s.ent.build_cumul();
        s.ent.build_lookup();
    }
    uint64_t nbytes = read_u64le(f);
    s.payload.resize((size_t)nbytes);
    return read_exact(f, s.payload.data(), s.payload.size());
}

static std::vector<uint8_t> decode_stream(const EncodedStream& s, Coder coder, int nsyms) {
    if (coder == Coder::Huffman) return decode_huffman(s.payload, s.huff_lengths, nsyms);
    if (coder == Coder::Rans) return decode_rans(s.payload, s.ent, nsyms);
    return decode_arithmetic(s.payload, s.ent, nsyms);
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

static bool parse_coder(const char* name, Coder& coder) {
    std::string s(name);
    if (s == "rans") { coder = Coder::Rans; return true; }
    if (s == "huffman" || s == "huff") { coder = Coder::Huffman; return true; }
    if (s == "arith" || s == "arithmetic") { coder = Coder::Arith; return true; }
    return false;
}

static int compress_file(Coder coder, const char* in_path, const char* out_path) {
    std::vector<uint16_t> pixels;
    if (!read_pixels(in_path, pixels)) return 1;

    std::vector<uint8_t> hi, lo;
    make_streams(pixels, hi, lo);
    EncodedStream enc_hi = encode_stream(hi, coder);
    EncodedStream enc_lo = encode_stream(lo, coder);

    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }

    fwrite(MAGIC, 1, 4, fout);
    write_u16le(fout, VERSION);
    write_u8(fout, (uint8_t)coder);
    write_u8(fout, 0);
    write_u32le(fout, WIDTH_DEFAULT);
    write_u32le(fout, HEIGHT_DEFAULT);
    write_u32le(fout, 2);
    write_stream(fout, enc_hi, coder);
    write_stream(fout, enc_lo, coder);
    fclose(fout);

    return 0;
}

static int decompress_file(const char* in_path, const char* out_path) {
    FILE* fin = fopen(in_path, "rb");
    if (!fin) { fprintf(stderr, "Cannot open input: %s\n", in_path); return 1; }

    uint8_t magic[4];
    if (!read_exact(fin, magic, 4) || memcmp(magic, MAGIC, 4) != 0) {
        fprintf(stderr, "Bad magic\n");
        fclose(fin);
        return 1;
    }
    uint16_t version = read_u16le(fin);
    if (version != VERSION) {
        fprintf(stderr, "Unsupported version %u\n", version);
        fclose(fin);
        return 1;
    }
    Coder coder = (Coder)read_u8(fin);
    read_u8(fin);
    uint32_t width = read_u32le(fin);
    uint32_t height = read_u32le(fin);
    uint32_t nstreams = read_u32le(fin);
    if (nstreams != 2) {
        fprintf(stderr, "Expected 2 streams\n");
        fclose(fin);
        return 1;
    }
    if (coder != Coder::Rans && coder != Coder::Huffman && coder != Coder::Arith) {
        fprintf(stderr, "Unknown coder id\n");
        fclose(fin);
        return 1;
    }

    EncodedStream enc_hi, enc_lo;
    if (!read_stream(fin, enc_hi, coder) || !read_stream(fin, enc_lo, coder)) {
        fprintf(stderr, "Short read while reading streams\n");
        fclose(fin);
        return 1;
    }
    fclose(fin);

    int npix = (int)(width * height);
    std::vector<uint8_t> hi = decode_stream(enc_hi, coder, npix);
    std::vector<uint8_t> lo = decode_stream(enc_lo, coder, npix);
    std::vector<uint16_t> pixels = reconstruct(hi, lo, (int)width, (int)height);

    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "Cannot open output: %s\n", out_path); return 1; }
    for (int i = 0; i < npix; i++) {
        uint8_t b[2] = {(uint8_t)(pixels[i] >> 8), (uint8_t)(pixels[i] & 0xFF)};
        fwrite(b, 1, 2, fout);
    }
    fclose(fout);
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr,
                "Usage:\n"
                "  %s c <rans|huffman|arith> <input> <output>\n"
                "  %s d <input> <output>\n",
                argv[0], argv[0]);
        return 1;
    }

    std::string mode(argv[1]);
    if (mode == "c") {
        if (argc != 5) {
            fprintf(stderr, "Usage: %s c <rans|huffman|arith> <input> <output>\n", argv[0]);
            return 1;
        }
        Coder coder;
        if (!parse_coder(argv[2], coder)) {
            fprintf(stderr, "Unknown coder: %s\n", argv[2]);
            return 1;
        }
        return compress_file(coder, argv[3], argv[4]);
    }
    if (mode == "d") {
        if (argc != 4) {
            fprintf(stderr, "Usage: %s d <input> <output>\n", argv[0]);
            return 1;
        }
        return decompress_file(argv[2], argv[3]);
    }

    fprintf(stderr, "Unknown mode: %s\n", argv[1]);
    return 1;
}
