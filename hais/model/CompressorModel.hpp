#pragma once
// Compressor-only model: LS solvers, residual generation, mode selection cost.

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Predictors.hpp"

// ---------------------------------------------------------------------------
// Generic N×N Gauss-Jordan solver for least-squares weight fitting.
// Solves XtX * w = Xty; writes result into w[N].
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
// Mode 4 — LS+bias(W, N, NW, 1): fit 4 weights (bias term included).
// ---------------------------------------------------------------------------

static void ls4_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                        float w[4], std::vector<uint16_t>& syms) {
    double XtX[4][4] = {}, Xty[4] = {};
    for (int y = 1; y < bh; y++)
        for (int x = 1; x < bw; x++) {
            double f[4] = { (double)blk[y*bw+(x-1)],
                            (double)blk[(y-1)*bw+x],
                            (double)blk[(y-1)*bw+(x-1)],
                            1.0 };
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
            uint16_t pred;
            if (y == 0 && x == 0)
                pred = (uint16_t)std::max(0.0f, std::min(65535.0f, w[3] + 0.5f));
            else if (y == 0)
                pred = (uint16_t)std::max(0.0f, std::min(65535.0f,
                           w[0]*blk[x-1] + w[3] + 0.5f));
            else if (x == 0)
                pred = (uint16_t)std::max(0.0f, std::min(65535.0f,
                           w[1]*blk[(y-1)*bw+x] + w[3] + 0.5f));
            else {
                float p = w[0]*blk[y*bw+(x-1)]
                        + w[1]*blk[(y-1)*bw+x]
                        + w[2]*blk[(y-1)*bw+(x-1)]
                        + w[3];
                pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
            }
            syms[idx] = zigzag((int16_t)(blk[idx] - pred));
        }
}

// ---------------------------------------------------------------------------
// Mode 6 — LS(W, N, NW, NE, NN): 5 spatial weights, +20 bytes overhead.
// ---------------------------------------------------------------------------

static void ls5_compute(const std::vector<uint16_t>& blk, int bw, int bh,
                        float w[5], std::vector<uint16_t>& syms) {
    double XtX[5][5] = {}, Xty[5] = {};
    for (int y = 2; y < bh; y++)
        for (int x = 1; x < bw; x++) {
            double ne = (x < bw-1) ? (double)blk[(y-1)*bw+(x+1)] : (double)blk[(y-1)*bw+x];
            double f[5] = { (double)blk[y*bw+(x-1)],
                            (double)blk[(y-1)*bw+x],
                            (double)blk[(y-1)*bw+(x-1)],
                            ne,
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
            uint16_t pred;
            if (y == 0 && x == 0)  pred = 0;
            else if (y == 0)       pred = blk[x - 1];
            else if (x == 0)       pred = blk[(y-1)*bw + x];
            else if (y == 1) {
                float ne = (x < bw-1) ? (float)blk[(y-1)*bw+(x+1)] : (float)blk[(y-1)*bw+x];
                float p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x]
                        + w[2]*blk[(y-1)*bw+(x-1)] + w[3]*ne;
                pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
            } else {
                float ne = (x < bw-1) ? (float)blk[(y-1)*bw+(x+1)] : (float)blk[(y-1)*bw+x];
                float p = w[0]*blk[y*bw+(x-1)] + w[1]*blk[(y-1)*bw+x]
                        + w[2]*blk[(y-1)*bw+(x-1)] + w[3]*ne
                        + w[4]*blk[(y-2)*bw+x];
                pred = (uint16_t)(int)std::max(0.0f, std::min(65535.0f, p + 0.5f));
            }
            syms[idx] = zigzag((int16_t)(blk[idx] - pred));
        }
}

// ---------------------------------------------------------------------------
// make_syms: zigzag residuals for modes 0 (raw), 1 (avg), 3 (global_mean).
// ---------------------------------------------------------------------------

static void make_syms(const std::vector<uint16_t>& blk, int mode,
                      int bw, int bh, std::vector<uint16_t>& syms,
                      uint16_t gmean) {
    syms.resize(bw * bh);
    for (int y = 0; y < bh; y++)
        for (int x = 0; x < bw; x++) {
            int      idx = y * bw + x;
            uint16_t pix = blk[idx];
            if (mode == 0)
                syms[idx] = pix;
            else if (mode == 3)
                syms[idx] = zigzag((int16_t)(pix - (int)gmean));
            else
                syms[idx] = zigzag((int16_t)(pix - avg_pred(blk, x, y, bw)));
        }
}

// ---------------------------------------------------------------------------
// byte_cost: H(hi8) + H(lo8) — entropy estimate for mode selection.
// ---------------------------------------------------------------------------

static double byte_cost(const std::vector<uint16_t>& syms) {
    uint64_t hfreq[256] = {}, lfreq[256] = {};
    for (uint16_t s : syms) { hfreq[s >> 8]++; lfreq[s & 0xFF]++; }
    double H = 0.0, N = (double)syms.size();
    for (int i = 0; i < 256; i++) {
        if (hfreq[i]) { double p = hfreq[i] / N; H -= p * std::log2(p); }
        if (lfreq[i]) { double p = lfreq[i] / N; H -= p * std::log2(p); }
    }
    return H;
}

