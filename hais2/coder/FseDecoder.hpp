#pragma once

#include <cstdint>
#include <vector>
#include "FseTable.hpp"

// ---------------------------------------------------------------------------
// Memory read helpers (advance pointer)
// ---------------------------------------------------------------------------

static inline uint8_t  pget_u8   (const uint8_t*& p) { return *p++; }
static inline uint16_t pget_u16le(const uint8_t*& p) {
    uint16_t v = (uint16_t)p[0] | ((uint16_t)p[1] << 8); p += 2; return v;
}
static inline uint32_t pget_u32le(const uint8_t*& p) {
    uint32_t v = (uint32_t)p[0] | ((uint32_t)p[1] << 8)
               | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
    p += 4; return v;
}
static inline uint64_t pget_u64le(const uint8_t*& p) {
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= (uint64_t)p[i] << (8 * i);
    p += 8; return v;
}

// ---------------------------------------------------------------------------
// Bit reader — feeds bits from a byte buffer in LSB-first order.
// ---------------------------------------------------------------------------

struct BitReader {
    const uint8_t* ptr;
    const uint8_t* end;
    uint64_t bitbuf = 0;
    int      bitcnt = 0;

    BitReader(const uint8_t* b, const uint8_t* e) : ptr(b), end(e) {}

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

// ---------------------------------------------------------------------------
// FSE decoder: reads symbols forward using the state machine.
// ---------------------------------------------------------------------------

struct FseDecoder {
    uint16_t  state;
    BitReader br;

    FseDecoder(const uint8_t* bits, const uint8_t* end, uint16_t s)
        : state(s), br(bits, end) {}

    uint8_t next(const FseTable& t) {
        const FseDecodeEntry& e = t.dec[state];
        uint8_t sym = e.sym;
        state = (uint16_t)(e.base + br.get(e.nb_bits));
        return sym;
    }
};

// ---------------------------------------------------------------------------
// Decode one FSE stream from ptr into out (n symbols).
// Advances ptr past the consumed bytes.
// ---------------------------------------------------------------------------

static void decode_stream(const uint8_t*& ptr, std::vector<uint8_t>& out, int n) {
    uint8_t flag = pget_u8(ptr);

    // Single-symbol: fill output with the stored symbol (flag=0x01, no bitstream).
    if (flag == 0x01) {
        out.assign(n, pget_u8(ptr));
        return;
    }

    FseTable tab;
    if (flag == 0) {
        for (int i = 0; i < 256; i++) tab.freq[i] = pget_u32le(ptr);
    } else {
        int nnz = (int)flag;  // 0x02..0xFE
        for (int k = 0; k < nnz; k++) {
            uint8_t sym = pget_u8(ptr);
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
