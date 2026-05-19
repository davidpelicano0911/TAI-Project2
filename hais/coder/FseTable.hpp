#pragma once

#include <cstdint>
#include <cstring>
#include <algorithm>

static constexpr uint32_t SCALE_BITS = 12;
static constexpr uint32_t SCALE      = 1u << SCALE_BITS;

struct FseDecodeEntry {
    uint8_t  sym;
    uint8_t  nb_bits;
    uint16_t base;
};

// FSE table: holds frequencies, decode array, and encode lookup arrays.
struct FseTable {
    uint32_t       freq[256]    = {};
    FseDecodeEntry dec[SCALE];
    uint16_t       enc_flat[SCALE];
    uint32_t       enc_offset[256];
    int            nb_base[256];

    void build() {
        // Spread symbols across SCALE slots deterministically.
        const uint32_t step = (SCALE >> 1) + (SCALE >> 3) + 3;
        uint8_t spread[SCALE];
        uint32_t pos = 0;
        for (int s = 0; s < 256; s++)
            for (uint32_t n = 0; n < freq[s]; n++) {
                spread[pos] = (uint8_t)s;
                pos = (pos + step) & (SCALE - 1);
            }

        // Build cumulative offsets for the encoder.
        uint32_t cumul = 0;
        for (int s = 0; s < 256; s++) {
            enc_offset[s] = cumul;
            cumul += freq[s];
            nb_base[s] = (freq[s] > 0)
                ? std::max(0, (int)(__builtin_clz(freq[s]) + SCALE_BITS - 32))
                : 0;
        }

        // Build decode and encode tables simultaneously.
        uint32_t next[256];
        for (int s = 0; s < 256; s++) next[s] = freq[s];
        for (uint32_t state = 0; state < SCALE; state++) {
            uint8_t  sym  = spread[state];
            uint32_t x    = next[sym]++;
            uint8_t  nb   = (uint8_t)(SCALE_BITS - (31u - __builtin_clz(x)));
            uint32_t base = (x << nb) - SCALE;
            dec[state]    = {sym, nb, (uint16_t)base};
            enc_flat[enc_offset[sym] + (x - freq[sym])] = (uint16_t)state;
        }
    }
};

