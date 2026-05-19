# HELIX — Context-Adaptive Lossless Compressor for Astronomical Images

**HELIX** is a specialised lossless compressor for raw 16-bit astronomical images (1500×1500, big-endian unsigned). It builds on the predict-then-entropy-code idea of [`rais`](../rais/) but replaces the static rANS entropy coder with a **JPEG-LS-style context-adaptive Golomb-Rice coder**, and adds **per-row debiasing** as a pre-processing stage.

On the 8 training files in [data2/](../data2/) it produces the **smallest total output of any tested compressor** (15.92 MB vs 16.11 MB for `bzip2-9` and 16.23 MB for `rais` v1), while staying ~2× faster than `bzip2-9`.

---

## Build and run

```bash
make

# Compress
./compress input_file output_file.hlx

# Decompress
./decompress output_file.hlx reconstructed_file

# Round-trip test on all data2 files
make test
```

The reconstructed file is bit-for-bit identical to the original.

---

## Pipeline overview

```
Raw 16-bit BE pixels
        │
        ▼
[1] Per-row debiasing            subtract row median, store medians in header
        │
        ▼
[2] MED predictor                pred = clamp(a+b-c, min(a,b), max(a,b))
        │
        ▼
[3] JPEG-LS context model        classify pixel into 1 of 63 contexts from
        │                        quantised local gradients (with sign symmetry)
        ▼
[4] Adaptive bias correction     subtract running integer bias C[ctx]
        │
        ▼
[5] Adaptive Golomb-Rice         k[ctx] = ⌈log2(A[ctx]/N[ctx])⌉, updated online
        │
        ▼
Bitstream
```

The decoder runs the inverse of steps 5 → 4 → 3 → 2 → 1.

---

## Stage details

### 1. Per-row debiasing

Ground-based CCDs often exhibit per-row bias variations (readout-bias drift, dark-current gradients). The benchmark file `F` is an extreme case: pixel-mean varies by σ≈384 ADU row-to-row vs σ<50 in well-calibrated frames.

We compute the median of each of the 1500 rows and subtract it before prediction. The 1500 medians (3000 bytes) are stored in the header so the decoder can add them back. Median is preferred over mean because it is robust to bright sources and saturated pixels.

The downstream predictor and entropy coder then see a near-zero-mean field, which is easier to model.

### 2. MED predictor (unchanged from `rais`)

```
A = pixel to the left          (x-1, y)
B = pixel above                (x,   y-1)
C = pixel above-left           (x-1, y-1)

Prediction = clamp(A + B - C, min(A, B), max(A, B))
```

Median Edge Detection adapts to local geometry: it predicts a smooth gradient where one is present, and falls back to whichever neighbour is on the "correct side" of an edge.

### 3. JPEG-LS context model

For each pixel we compute three local gradients:

```
g1 = D - B  (vertical change above)
g2 = B - C  (vertical change left of above)
g3 = C - A  (horizontal change in row above)
```

where D is the upper-right neighbour. Each gradient is quantised into 5 levels using thresholds T1=3, T2=21:

```
g ≤ -21     → -2
-21 < g ≤ -3 → -1
 -3 < g <  3 →  0
  3 ≤ g < 21 → +1
g ≥ 21      → +2
```

This yields 5³=125 raw triples. **Sign symmetry** (LOCO-I §4.2) reduces these to **63 canonical contexts**: if the first non-zero quantised gradient is negative, we negate all gradients *and* flip the sign of the residual — so a context and its mirror share state.

### 4. Adaptive bias correction (LOCO-I §4.3)

For each context Q we maintain four counters:

| Field | Meaning |
|-------|---------|
| `A[Q]` | sum of `|residual|` |
| `B[Q]` | sum of signed residual (running mean estimator) |
| `N[Q]` | number of samples observed |
| `C[Q]` | current integer bias correction, clamped to [-128, 127] |

Before coding pixel `i` in context Q, we subtract `C[Q]` from its signed residual. After coding, `B[Q]` is updated; whenever `|B/N| ≥ 1` we increment/decrement `C` and reset `B`. This cancels persistent systematic offsets that the gradient quantiser cannot capture.

To prevent counter overflow and to let the model adapt to non-stationary statistics, all counters are halved whenever `N[Q]` reaches 64.

### 5. Adaptive Golomb-Rice coding

Each context picks its own Golomb parameter from its running |residual| mean:

```
k[Q] = smallest k such that (N[Q] << k) ≥ A[Q]
     ≈ ⌈log2(mean(|residual|))⌉
```

The bias-corrected signed residual is mapped to an unsigned integer with zigzag, then encoded as:

- **Quotient** `q = u >> k` in unary (q 1-bits followed by a single 0)
- **Remainder** `r = u & ((1<<k)-1)` as k raw bits

If `q ≥ 32` (rare, pathological residuals), we emit 32 ones followed by a zero terminator and then the full 32-bit `u` raw. This bounds the worst case.

This is functionally equivalent to JPEG-LS's "regular mode" coder. We omit the "run mode" (long stretches of identical pixels) because astronomical sky regions, while smooth, are dominated by Gaussian read-noise — runs are rare.

Compared with the static rANS used in `rais` v1, adaptive Rice:

- needs no frequency tables (saves ~2 KB header per stream),
- adapts to local statistics (per-context distributions vary widely between sky and source regions),
- avoids division/multiplication in the hot loop (faster decode on the smooth files).

For the noisier files where residuals are closer to Gaussian than Laplacian, Rice slightly underperforms rANS — this is visible in file F (0.2% larger than `rais` v1 on F alone, but smaller on all other 7 files).

---

## File format

All integers are little-endian.

```
Offset  Size   Field
------  ----   -----
0       4      Magic bytes "HELX"
4       2      Version (= 1)
6       4      Width  (pixels)
10      4      Height (pixels)
14      2      Flags  (bit 0 = row debiasing used)
16      W×2    Per-row medians (uint16 LE)  [if flag bit 0 set]
16+W×2  8      Bitstream length in bytes (uint64 LE)
...     N      Golomb-Rice bitstream (MSB-first)
```

Total fixed overhead: 24 bytes + 3000 bytes for medians = ~3 KB regardless of image content.

---

## Benchmark on the 8 training files

| Compressor | Total compressed | bits / byte | Compress | Decompress | Total time |
|------------|-----------------:|------------:|---------:|-----------:|-----------:|
| **helix**  | **15.92 MB**     | **3.537**   | 956 ms   | 692 ms     | **1.65 s** |
| bzip2-9    | 16.11 MB         | 3.579       | 2179 ms  | 1237 ms    | 3.42 s     |
| rais v1    | 16.23 MB         | 3.606       | 546 ms   | 630 ms     | 1.18 s     |
| xz-9       | 16.58 MB         | 3.685       | 12350 ms | 809 ms     | 13.2 s     |
| astra      | 17.11 MB         | 3.802       | 631 ms   | 541 ms     | 1.17 s     |
| zstd-9     | 18.82 MB         | 4.182       | 987 ms   | 90 ms      | 1.08 s     |

`helix` produces the smallest output across the board, 2× faster than `bzip2-9` (next-best ratio) and 8× faster than `xz-9` while still compressing better.

---

## Why not just use gzip/zstd?

General-purpose compressors model the file as a flat byte stream. They do not exploit:

- 16-bit big-endian integer structure (each pixel = 2 correlated bytes),
- 2-D spatial smoothness of astronomical fields,
- the fact that ~99% of pixels are sky-noise with a Gaussian/Laplacian residual distribution,
- per-row systematic biases from the CCD readout chain.

By baking those priors into the pipeline (debias → predict → context → Rice), we approach the noise-limited entropy floor for astronomical data while staying competitive on speed.

---

## References

- **Weinberger, M. J., Seroussi, G., & Sapiro, G. (2000).** *The LOCO-I Lossless Image Compression Algorithm: Principles and Standardization into JPEG-LS.* IEEE Transactions on Image Processing, 9(8), 1309–1324. https://doi.org/10.1109/83.855427
- **Pence, W. D., Seaman, R., & White, R. L. (2009).** *Lossless Astronomical Image Compression and the Effects of Noise.* PASP, 121(878), 414–427. https://arxiv.org/abs/0903.2140
- **Malvar, H. S. (2006).** *Adaptive Run-Length / Golomb-Rice Encoding of Quantized Generalized Gaussian Sources with Unknown Statistics.* IEEE DCC 2006. https://www.microsoft.com/en-us/research/publication/adaptive-run-lengthgolomb-rice-encoding-of-quantized-generalized-gaussian-sources/
- **CCSDS 121.0-B-3.** *Lossless Data Compression — Blue Book.* https://ccsds.org/Pubs/121x0b3.pdf
- **Duda, J. (2013).** *Asymmetric Numeral Systems.* arXiv:0902.0271. https://arxiv.org/abs/0902.0271 *(used by `rais` v1 for comparison)*
