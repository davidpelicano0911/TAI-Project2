#pragma once
// Shared predictor utilities used by both compressor and decompressor.

#include <cstdint>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Zigzag: signed int16 <-> unsigned uint16
//   encode:  v >= 0 -> 2v        v < 0 -> -2v - 1
//   decode:  even  -> u/2        odd   -> -(u+1)/2
// ---------------------------------------------------------------------------

static inline uint16_t zigzag(int16_t v) {
    return (v >= 0) ? (uint16_t)(v * 2) : (uint16_t)((-v) * 2 - 1);
}
static inline int16_t zagzig(uint16_t u) {
    return (u & 1) ? -(int16_t)((u + 1) / 2) : (int16_t)(u / 2);
}

// ---------------------------------------------------------------------------
// avg predictor: (W + N + NW + 1) / 3
// Causal neighbours: A=left, B=above, C=above-left (all already decoded).
// ---------------------------------------------------------------------------

static inline uint16_t avg_pred(const std::vector<uint16_t>& blk,
                                int x, int y, int bw) {
    if (y == 0 && x == 0) return 0;
    if (y == 0)            return blk[x - 1];
    if (x == 0)            return blk[(y - 1) * bw + x];
    int A = blk[y * bw + (x - 1)];
    int B = blk[(y - 1) * bw + x];
    int C = blk[(y - 1) * bw + (x - 1)];
    return (uint16_t)((A + B + C + 1) / 3);
}
