#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>
#include "FseTable.hpp"

// ---------------------------------------------------------------------------
// Normalize raw counts into FSE frequencies that sum to SCALE.
// ---------------------------------------------------------------------------

static void fit_freqs(const uint64_t* cnt, uint32_t* freq) {
    uint64_t total = 0;
    for (int i = 0; i < 256; i++) total += cnt[i];
    uint32_t used = 0;
    for (int i = 0; i < 256; i++) {
        if (!cnt[i]) { freq[i] = 0; continue; }
        freq[i] = (uint32_t)std::max((uint64_t)1, cnt[i] * (uint64_t)SCALE / total);
        used += freq[i];
    }
    int peak = 0;
    for (int i = 1; i < 256; i++) if (cnt[i] > cnt[peak]) peak = i;
    if (used < SCALE) freq[peak] += SCALE - used;
    else              freq[peak] -= used - SCALE;
}

// ---------------------------------------------------------------------------
// Buffer write helpers
// ---------------------------------------------------------------------------

static inline void bput8 (std::vector<uint8_t>& b, uint8_t  v) { b.push_back(v); }
static inline void bput32(std::vector<uint8_t>& b, uint32_t v) {
    for (int i = 0; i < 4; i++) { b.push_back(v & 0xFF); v >>= 8; }
}
static inline void bput64(std::vector<uint8_t>& b, uint64_t v) {
    for (int i = 0; i < 8; i++) { b.push_back(v & 0xFF); v >>= 8; }
}

// ---------------------------------------------------------------------------
// Bit writer — stores (bits, nbits) ops and flushes them in reverse order,
// so the decoder can read symbols in the original forward order.
// ---------------------------------------------------------------------------

struct BitWriter {
    struct Op { uint16_t bits; uint8_t nbits; };
    std::vector<Op> ops;
    std::vector<uint8_t> bytes;

    void put(uint32_t bits, int n) { ops.push_back({(uint16_t)bits, (uint8_t)n}); }

    void flush_reverse() {
        uint64_t bitbuf = 0; int bitcnt = 0;
        for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
            bitbuf |= (uint64_t)it->bits << bitcnt;
            bitcnt += it->nbits;
            while (bitcnt >= 8) {
                bytes.push_back((uint8_t)(bitbuf & 0xFF));
                bitbuf >>= 8; bitcnt -= 8;
            }
        }
        if (bitcnt > 0) bytes.push_back((uint8_t)(bitbuf & 0xFF));
    }
};

// ---------------------------------------------------------------------------
// FSE encoder: pushes symbols in reverse, emits bits via BitWriter.
// Output layout: [state_lo][state_hi][bits...]
// ---------------------------------------------------------------------------

struct FseEncoder {
    uint16_t state = 0;
    BitWriter bw;

    void push(uint8_t sym, const FseTable& t) {
        uint32_t xs = (uint32_t)state + SCALE;
        int nb = t.nb_base[sym];
        if ((xs >> nb) >= 2 * t.freq[sym]) ++nb;
        bw.put(xs & ((1u << nb) - 1u), nb);
        state = t.enc_flat[t.enc_offset[sym] + (xs >> nb) - t.freq[sym]];
    }

    void finish(std::vector<uint8_t>& out) {
        bw.flush_reverse();
        out.push_back((uint8_t)(state & 0xFF));
        out.push_back((uint8_t)(state >> 8));
        out.insert(out.end(), bw.bytes.begin(), bw.bytes.end());
    }
};

// ---------------------------------------------------------------------------
// Encode a byte stream into buf.
// Format: [flag:1B][table_data][size:8B][state:2B + bits]
//
// Flag values:
//   0x00          — dense  (256 × u32le frequencies)
//   0x01          — single-symbol: just one sym byte, no table, no bitstream
//   0x02..0xFE    — sparse (nnz × {sym:u8, freq:u32le}), nnz = flag
// ---------------------------------------------------------------------------

static void encode_stream(std::vector<uint8_t>& buf, const std::vector<uint8_t>& bytes) {
    uint64_t cnt[256] = {};
    for (uint8_t b : bytes) cnt[b]++;
    int nnz = 0;
    for (int i = 0; i < 256; i++) if (cnt[i]) nnz++;

    // Single-symbol: no entropy needed, just store the symbol.
    // Saves 13 bytes vs sparse-with-1-entry (no table, no bitstream, no nbytes).
    if (nnz == 1) {
        int sym = 0; while (!cnt[sym]) sym++;
        bput8(buf, 0x01);
        bput8(buf, (uint8_t)sym);
        return;
    }

    FseTable tab;
    fit_freqs(cnt, tab.freq);
    tab.build();
    if (nnz <= 254) {
        bput8(buf, (uint8_t)nnz);  // flag = nnz = 2..254 = 0x02..0xFE
        for (int i = 0; i < 256; i++) {
            if (!tab.freq[i]) continue;
            bput8(buf, (uint8_t)i);
            bput32(buf, tab.freq[i]);
        }
    } else {
        // Dense format for nnz=255 or nnz=256 — stores all 256 freqs.
        bput8(buf, 0);
        for (int i = 0; i < 256; i++) bput32(buf, tab.freq[i]);
    }

    FseEncoder enc;
    enc.bw.ops.reserve(bytes.size());
    std::vector<uint8_t> stream;
    for (int i = (int)bytes.size() - 1; i >= 0; i--) enc.push(bytes[i], tab);
    enc.finish(stream);
    bput64(buf, (uint64_t)stream.size());
    buf.insert(buf.end(), stream.begin(), stream.end());
}
