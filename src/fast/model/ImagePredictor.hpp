#pragma once

#include <algorithm>
#include <cstdint>

// JPEG-LS MED (median-edge detector) predictor for 16-bit samples.
inline uint16_t med_predict(uint16_t W, uint16_t N, uint16_t NW) {
    if (NW >= std::max(W, N)) return std::min(W, N);
    if (NW <= std::min(W, N)) return std::max(W, N);
    return static_cast<uint16_t>(W + N - NW);
}

// Modular zigzag: maps wrapping residual u (= pixel - pred, mod 2^16) so that
// small magnitudes (both positive and negative) map to values near 0.
//   0 → 0,  65535 → 1,  1 → 2,  65534 → 3,  2 → 4, ...
inline uint16_t zigzag_encode(uint16_t u) {
    return (u <= 32767u) ? static_cast<uint16_t>(u * 2u)
                         : static_cast<uint16_t>((65536u - u) * 2u - 1u);
}

// Inverse modular zigzag.
inline uint16_t zigzag_decode(uint16_t z) {
    return (z & 1u) ? static_cast<uint16_t>(65536u - (z + 1u) / 2u)
                    : static_cast<uint16_t>(z / 2u);
}
