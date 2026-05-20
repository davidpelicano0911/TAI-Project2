# HAIS2 Fine-Tuning Notes

## Overview

HAIS2 is a copy of HAIS with additional predictor candidates and a better cost estimator.
**Result: 3.3962 b/B vs HAIS baseline 3.3999 b/B** (−16 KB, −0.0037 b/B).

---

## Changes vs HAIS

### 1. New Predictors (model/CompressorModel.hpp + Predictors.hpp)

| Mode | Name | Neighbors | Overhead | Wins on |
|------|------|-----------|----------|---------|
| 7 | GAP | W, N, NW (gradient-adaptive blend) | 0 bytes | Never |
| 8 | MED | W, N, NW (JPEG-LS median) | 0 bytes | Never |
| 9 | **LS6** | W, N, NW, WW, NN + bias (6 weights) | 24 bytes | C, D, F, G |

LS6 is a 2nd-order least-squares predictor. It captures smooth spatial ramps and
gradients better than the previous LS4/LS5 candidates, especially for files with
near-uniform backgrounds (C, D, G) and high-gradient images (F).

### 2. Better Mode Selection Cost (`ctx_cost` replaces `byte_cost`)

`ctx_cost` conditions the lo-byte entropy estimate on hi==0 vs hi>0:

```
H = H(hi) + (n0/N)*H(lo | hi==0) + (n1/N)*H(lo | hi!=0)
```

This more accurately predicts the actual compressed size when the hi==0 fraction
is large (smooth regions), giving mode selection a more realistic signal.

### 3. Candidate Set Reduced to 6

Removed GAP (mode 7) and MED (mode 8) from the competition — analysis confirmed
they are consistently worse than `avg` across all 8 test images. Keeping them only
adds compute without changing the winner.

Final candidates: **raw, avg, mean, LS4, LS5, LS6**
(LS4 and LS5 are also effectively dominated but kept for generality.)

---

## Analysis Findings

### Per-file mode distribution (1500×1500, 3×3 blocks, BS=500)

| File | Winner | Bits/px | Notes |
|------|--------|---------|-------|
| A | avg (9/9) | 4.39 b/B | Star field; avg=8.77 bits/px, near-optimal |
| B | avg+mean | 3.24 b/B | Mixed; raw/mean competitive |
| C | mean (8/9) + LS6 (1/9) | 2.43 b/B | Smooth background |
| D | mean (6/9) + LS6 (3/9) | 2.73 b/B | Smooth gradient |
| E | raw (7/9) + mean (2/9) | 2.57 b/B | Sparse stars |
| F | LS6 (9/9) | 6.11 b/B | High-contrast; strong spatial gradients |
| G | mean (6/9) + avg+LS6 | 2.41 b/B | Smooth with features |
| H | mean (9/9) | 3.29 b/B | Uniform background |

### Why GAP/MED never win

For astronomical images, `avg = (W+N+NW)/3` already outperforms gradient-adaptive
predictors. GAP/MED are designed for natural images with sharp edges; astronomical
data is either smooth (where `mean` wins) or random-looking star fields (where no
local predictor helps vs `avg`).

### Why smaller blocks don't help

Tested BS=250 (6×6=36 blocks): result 3.4247 b/B — **worse**.

The FSE lo-stream table is ~1 KB per block regardless of block size. Switching from
9 blocks (BS=500) to 36 blocks (BS=250) adds ~27 KB/image × 8 files = 216 KB of
extra table headers, which exceeds the ~59 KB entropy gain from better local LS
adaptation.

**Optimal block size for 1500×1500: BS=500** (3×3 = 9 blocks).

### Context (ctx) mode analysis

The `ctx` flag splits the lo stream into `lo_zero` (hi==0) and `lo_nonzero` (hi>0)
sub-streams. Never triggered in practice because:

- Smooth files (C, D): 99.9% hi==0 → `lo_nonzero` sub-stream is empty; 2nd table
  header costs more than entropy gain.
- Complex files (A): 80% hi==0 → tested; ctx payload 394 bytes larger than plain.
- High-entropy files (F): ~0% hi==0 → not eligible.

### Residual entropy gap

The hi/lo byte split introduces a small coding overhead vs true 16-bit symbol coding:

| File | H(hi)+H(lo) vs H(16-bit) gap | Bytes "left on table" |
|------|------------------------------|-----------------------|
| A | 0.028 bits/px | ~79 KB |
| F | 0.019 bits/px | ~54 KB |
| B | 0.007 bits/px | ~20 KB |

Total theoretical gain from 16-bit coding: ~22 KB across all 8 files. Not worth the
65536-entry FSE table complexity.

### Remaining residual correlation

After `avg` prediction, lag-1 autocorrelation of residuals for file A = **0.79**
(raw = 0.95). This is high and indicates spatial structure that a context-adaptive
arithmetic coder (like PRISM v3) could exploit — but PRISM's single-predictor
approach still loses to HAIS's multi-predictor block selection overall.

---

## What Was Tried But Reverted

| Idea | Result |
|------|--------|
| 6-neighbor GAP (W,WW,N,NN,NW,NE) | Worse — MED fallback hurts smooth images |
| Robust GAP (hi≥8 outlier rejection) | Neutral — stars rarely trigger the threshold |
| Smaller blocks (BS=250, BS=150) | Worse — FSE header overhead > entropy gain |
| LS7 (adds NE to LS6) | Tiny gain (~0.001 b/B) on some files, not worth complexity |
| ctx split for smooth files | Never saves bytes — table overhead dominates |

---

## Compression Limits

For the current predictor set, we are already very close to the theoretical optimum:

- File A: achieved 8.77 bits/px vs theoretical residual entropy 8.70 bits/px (0.07 gap)
- File B: achieved ~6.43 bits/px vs theoretical 6.43 bits/px (essentially optimal)
- File F: achieved 12.16 bits/px vs theoretical 12.17 bits/px (essentially optimal)

Further gains require fundamentally better predictors (e.g., CALIC-style context bias
correction, or a learned/adaptive predictor) rather than tuning the entropy coder.
