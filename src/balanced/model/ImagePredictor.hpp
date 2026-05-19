#pragma once

#include <algorithm>
#include <cstdint>

// JPEG-LS MED (median-edge detector) predictor for 16-bit samples.
inline uint16_t med_predict(uint16_t W, uint16_t N, uint16_t NW) {
    if (NW >= std::max(W, N)) return std::min(W, N);
    if (NW <= std::min(W, N)) return std::max(W, N);
    return static_cast<uint16_t>(W + N - NW);
}

// GAP (Gradient-Adjusted Predictor) for 16-bit samples.
// Measures horizontal activity dh and vertical activity dv using a 6-pixel
// window. Picks the predictor from the more stable axis; falls back to MED.
static constexpr int GAP_THRESH = 128;

inline uint16_t gap_predict(uint16_t W, uint16_t WW,
                             uint16_t N, uint16_t NN,
                             uint16_t NW, uint16_t NE) {
    int dh = std::abs((int)W - (int)WW)
           + std::abs((int)N - (int)NW)
           + std::abs((int)N - (int)NE);
    int dv = std::abs((int)N - (int)NN)
           + std::abs((int)W - (int)NW)
           + std::abs((int)NE - (int)N);
    if (dv - dh > GAP_THRESH) return W;
    if (dh - dv > GAP_THRESH) return N;
    return med_predict(W, N, NW);
}

// Robust GAP: if a spatial neighbour had a large residual (hi byte >= threshold)
// it was likely a star or bright outlier. Using it in the prediction creates a
// strong negative lag-1 correlation in residuals (GAP over-predicts the next
// background pixel). Fall back to the clean neighbour instead.
// For most images hi_W == hi_N == 0 always, so this never fires and behaviour
// is identical to plain gap_predict.
static constexpr uint8_t ROBUST_HI_THRESH = 8;

inline uint16_t robust_gap_predict(uint16_t W, uint16_t WW,
                                    uint16_t N, uint16_t NN,
                                    uint16_t NW, uint16_t NE,
                                    uint8_t hi_W, uint8_t hi_N) {
    const bool w_bad = (hi_W >= ROBUST_HI_THRESH);
    const bool n_bad = (hi_N >= ROBUST_HI_THRESH);
    if (w_bad && !n_bad) return N;
    if (!w_bad && n_bad) return W;
    if (w_bad && n_bad)  return NW;
    return gap_predict(W, WW, N, NN, NW, NE);
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
