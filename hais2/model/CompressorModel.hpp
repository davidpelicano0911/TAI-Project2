#pragma once
// HAIS2 compressor model — adds GAP, MED, LS6, LOCO-I bias correction,
// and finer context-aware cost estimation.

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Predictors.hpp"

// ---------------------------------------------------------------------------
// Generic N×N Gauss-Jordan solver
// ---------------------------------------------------------------------------

template<int N>
static void ls_solve(double XtX[N][N], double Xty[N], float w[N]) {
    double M[N][N + 1];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) M[i][j] = XtX[i][j];
        M[i][N] = Xty[i];
    }
    for (int col = 0; col < N; col++) {
        int piv = col;
        for (int r = col + 1; r < N; r++)
            if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
        if (piv != col)
            for (int k = 0; k <= N; k++) std::swap(M[col][k], M[piv][k]);
        double d = M[col][col];
        if (std::abs(d) < 1e-8) { w[col] = 0.0f; continue; }
        for (int r = 0; r < N; r++) {
            if (r == col) continue;
            double f = M[r][col] / d;
            for (int k = col; k <= N; k++) M[r][k] -= f * M[col][k];
        }
    }
    for (int i = 0; i < N; i++)
        w[i] = (std::abs(M[i][i]) > 1e-8) ? (float)(M[i][N] / M[i][i]) : 0.0f;
}

// ---------------------------------------------------------------------------
// Mode 4 — LS+bias(W, N, NW, 1): 4 weights
// ---------------------------------------------------------------------------

static void ls4_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                        float w[4], std::vector<uint16_t>& syms) {
    double XtX[4][4] = {}, Xty[4] = {};
    for (int y = 1; y < bh; y++)
        for (int x = 1; x < bw; x++) {
            double f[4] = { (double)blk[y*bw+(x-1)], (double)blk[(y-1)*bw+x],
                            (double)blk[(y-1)*bw+(x-1)], 1.0 };
            double t = blk[y * bw + x];
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) XtX[i][j] += f[i] * f[j];
                Xty[i] += f[i] * t;
            }
        }
    ls_solve<4>(XtX, Xty, w);
    syms.resize(bw * bh);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int idx = y * bw + x;
            float p;
            if (y == 0 && x == 0) p = w[3];
            else if (y == 0)      p = w[0]*blk[x-1] + w[3];
            else if (x == 0)      p = w[1]*blk[(y-1)*bw+x] + w[3];
            else p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x]
                   + w[2]*blk[(y-1)*bw+(x-1)] + w[3];
            syms[idx] = zigzag((int16_t)(blk[idx] - (uint16_t)(int)std::max(0.0f,std::min(65535.0f,p+0.5f))));
        }
}

// ---------------------------------------------------------------------------
// Mode 6 — LS(W, N, NW, NE, NN): 5 spatial weights
// ---------------------------------------------------------------------------

static void ls5_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                        float w[5], std::vector<uint16_t>& syms) {
    double XtX[5][5] = {}, Xty[5] = {};
    for (int y = 2; y < bh; y++)
        for (int x = 1; x < bw; x++) {
            double ne = (x < bw-1) ? (double)blk[(y-1)*bw+(x+1)] : (double)blk[(y-1)*bw+x];
            double f[5] = { (double)blk[y*bw+(x-1)], (double)blk[(y-1)*bw+x],
                            (double)blk[(y-1)*bw+(x-1)], ne,
                            (double)blk[(y-2)*bw+x] };
            double t = blk[y * bw + x];
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) XtX[i][j] += f[i] * f[j];
                Xty[i] += f[i] * t;
            }
        }
    ls_solve<5>(XtX, Xty, w);
    syms.resize(bw * bh);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int idx = y * bw + x;
            float p;
            if (y == 0 && x == 0) p = 0;
            else if (y == 0)      p = (float)blk[x-1];
            else if (x == 0)      p = (float)blk[(y-1)*bw+x];
            else if (y == 1) {
                float ne = (x < bw-1) ? (float)blk[(y-1)*bw+(x+1)] : (float)blk[(y-1)*bw+x];
                p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x]
                  + w[2]*blk[(y-1)*bw+(x-1)] + w[3]*ne;
            } else {
                float ne = (x < bw-1) ? (float)blk[(y-1)*bw+(x+1)] : (float)blk[(y-1)*bw+x];
                p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x]
                  + w[2]*blk[(y-1)*bw+(x-1)] + w[3]*ne
                  + w[4]*blk[(y-2)*bw+x];
            }
            syms[idx] = zigzag((int16_t)(blk[idx] - (uint16_t)(int)std::max(0.0f,std::min(65535.0f,p+0.5f))));
        }
}

// ---------------------------------------------------------------------------
// Mode 9 — LS+bias(W, N, NW, WW, NN, 1): 6 weights (2nd-order)
// Better on smooth ramps; 24-byte overhead.
// ---------------------------------------------------------------------------

static void ls6_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                        float w[6], std::vector<uint16_t>& syms) {
    double XtX[6][6] = {}, Xty[6] = {};
    for (int y = 2; y < bh; y++)
        for (int x = 2; x < bw; x++) {
            double f[6] = { (double)blk[y*bw+(x-1)],
                            (double)blk[(y-1)*bw+x],
                            (double)blk[(y-1)*bw+(x-1)],
                            (double)blk[y*bw+(x-2)],    // WW
                            (double)blk[(y-2)*bw+x],    // NN
                            1.0 };
            double t = blk[y * bw + x];
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) XtX[i][j] += f[i] * f[j];
                Xty[i] += f[i] * t;
            }
        }
    ls_solve<6>(XtX, Xty, w);
    syms.resize(bw * bh);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int idx = y * bw + x;
            float p;
            if (y < 2 || x < 2) {
                // Fallback to GAP for border pixels
                p = (float)gap_pred(blk, x, y, bw);
            } else {
                p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x]
                  + w[2]*blk[(y-1)*bw+(x-1)] + w[3]*blk[y*bw+(x-2)]
                  + w[4]*blk[(y-2)*bw+x] + w[5];
            }
            syms[idx] = zigzag((int16_t)(blk[idx] - (uint16_t)(int)std::max(0.0f,std::min(65535.0f,p+0.5f))));
        }
}

// ---------------------------------------------------------------------------
// make_syms: residuals for simple modes 0, 1, 3, 7 (GAP), 8 (MED)
// ---------------------------------------------------------------------------

static void make_syms(const std::vector<uint16_t>& blk, int mode,
                      int bw, int bh, std::vector<uint16_t>& syms,
                      uint16_t gmean) {
    syms.resize(bw * bh);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int      idx = y * bw + x;
            uint16_t pix = blk[idx];
            uint16_t pred;
            if      (mode == 0) pred = pix;          // raw
            else if (mode == 3) pred = gmean;         // global mean
            else if (mode == 7) pred = gap_pred(blk, x, y, bw);
            else if (mode == 8) pred = med_pred(blk, x, y, bw);
            else                pred = avg_pred(blk, x, y, bw);
            syms[idx] = (mode == 0) ? pix : zigzag((int16_t)(pix - (int)pred));
        }
}

// ---------------------------------------------------------------------------
// ctx_cost: entropy estimate using hi==0 / hi>0 context split for lo byte.
// Closer to what encode_stream actually produces with context coding.
// ---------------------------------------------------------------------------

static double ctx_cost(const std::vector<uint16_t>& syms) {
    uint64_t hfreq[256] = {};
    uint64_t l0freq[256] = {}, l1freq[256] = {};
    uint64_t n0 = 0, n1 = 0;
    for (uint16_t s : syms) {
        uint8_t hi = s >> 8, lo = s & 0xFF;
        hfreq[hi]++;
        if (hi == 0) { l0freq[lo]++; n0++; }
        else         { l1freq[lo]++; n1++; }
    }
    double N = (double)syms.size();
    double H = 0.0;
    // hi entropy
    for (int i = 0; i < 256; i++)
        if (hfreq[i]) { double p = hfreq[i] / N; H -= p * std::log2(p); }
    // lo entropy conditioned on hi==0
    if (n0 > 0) for (int i = 0; i < 256; i++)
        if (l0freq[i]) { double p = l0freq[i]/(double)n0; H -= (n0/N)*p*std::log2(p); }
    // lo entropy conditioned on hi>0
    if (n1 > 0) for (int i = 0; i < 256; i++)
        if (l1freq[i]) { double p = l1freq[i]/(double)n1; H -= (n1/N)*p*std::log2(p); }
    return H;
}
