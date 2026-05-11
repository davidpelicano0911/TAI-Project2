# RAIS — Raw Astronomical Image compressor

**RAIS** (Raw Astronomical Image compressor) is a lossless compressor designed specifically for raw 16-bit astronomical images (1500×1500 pixels, big-endian). It combines two well-known techniques: a **spatial predictor** and a **rANS entropy coder**.

---

## How to build and run

```bash
make

# Compress
./compress input_file output_file.rais

# Decompress
./decompress output_file.rais reconstructed_file
```

The reconstructed file is guaranteed to be **bit-for-bit identical** to the original.

---

## How it works — step by step

The pipeline has three stages:

```
Raw pixels  →  [1. Predict]  →  Residuals  →  [2. Zigzag]  →  Unsigned values  →  [3. rANS]  →  Compressed bytes
```

---

### Stage 1 — MED Predictor

**The idea:** neighbouring pixels in astronomical images are very similar. Instead of storing each pixel value directly, we store only the *error* (residual) between the actual pixel and a prediction made from its neighbours. These errors are tiny numbers concentrated around zero — much easier to compress than the original values.

**The predictor used is MED (Median Edge Detection)**, the same one used in JPEG-LS. For each pixel at position (x, y):

```
A = pixel to the left      (x-1, y)
B = pixel above            (x,   y-1)
C = pixel above-left       (x-1, y-1)
```

We compute a linear prediction `A + B - C`, then clamp it to `[min(A,B), max(A,B)]`. The clamping makes it behave like an edge detector: near a vertical edge it falls back to B (above), near a horizontal edge it falls back to A (left).

**Example** (pixel values in a smooth sky region):

```
C=1000  B=1005
A=1003  ?=1008

Prediction = A + B - C = 1003 + 1005 - 1000 = 1008
Residual   = actual - prediction = 1008 - 1008 = 0   ← tiny!
```

Without prediction, we would store `1008` (needs ~10 bits). With prediction, we store `0` (needs ~1 bit after entropy coding).

**Why MED specifically?**
A plain left-neighbour predictor works well on smooth gradients but breaks down at edges. MED adapts: when A, B, and C suggest an edge is present, it automatically switches to the better direction.

> **Reference:** Weinberger, M. J., Seroussi, G., & Sapiro, G. (2000). *The LOCO-I Lossless Image Compression Algorithm: Principles and Standardization into JPEG-LS.* IEEE Transactions on Image Processing, 9(8), 1309–1324.
> [https://doi.org/10.1109/83.855427](https://doi.org/10.1109/83.855427)

---

### Stage 2 — Zigzag mapping

The residuals are **signed** integers (e.g., −500, 0, +3, −1). Entropy coders work on **unsigned** symbols. We use a zigzag mapping to fold signed integers into non-negative ones, keeping small absolute values mapped to small unsigned values:

```
0  →  0
-1 →  1
+1 →  2
-2 →  3
+2 →  4
-3 →  5
...
```

Formula:
- If `r >= 0`: mapped value = `2 * r`
- If `r < 0`:  mapped value = `2 * |r| - 1`

**Example:**

```
Residual = -1  →  zigzag(−1) = 1   (very common near zero)
Residual = +2  →  zigzag(+2) = 4
Residual = 500 →  zigzag(500) = 1000
```

After zigzag, the distribution is heavily skewed towards 0. Most pixels in a smooth astronomical image produce residual 0 or ±1.

> **Reference:** This technique is widely used in codecs. It appears in, e.g., Protocol Buffers (Google) variable-length integer encoding and in FLIF's maniac tree coder.

---

### Stage 3 — rANS entropy coding

**ANS (Asymmetric Numeral Systems)** is a modern entropy coder that achieves compression ratios close to the theoretical Shannon entropy limit, with speed comparable to Huffman coding and accuracy closer to arithmetic coding.

**The intuition:** imagine a number `state` that encodes the entire history of symbols seen so far. To encode a symbol with frequency `f` out of total `M`:

- Renormalise `state` by pushing bytes out until it fits in a safe range.
- Update: `state = (state / f) * M + cumulative_freq(symbol) + (state % f)`

To decode, we do the reverse: read `state % M` to look up which symbol it falls in, recover the symbol, then undo the update.

**Example** (simplified, M = 8, two symbols: `A` with freq 6, `B` with freq 2):

```
Symbol table:  slots 0-5 → A,  slots 6-7 → B

Encode 'A' from state=10:
  slot = 10 % 8 = 2  (but we're encoding, not decoding)
  new state = (10 / 6) * 8 + 0 + (10 % 6) = 1*8 + 0 + 4 = 12

Decode from state=12:
  slot = 12 % 8 = 4  → symbol A  (slot 4 is in range 0-5)
  new state = 6 * (12 / 8) + 4 - 0 = 6*1 + 4 = 10  ← recovered!
```

The key property: **common symbols shrink the state slightly; rare symbols grow it a lot**. Over millions of symbols, the state size converges to exactly the Shannon entropy.

We use **rANS (range ANS)** with scale M = 65536 (16-bit). Frequencies are computed from the data in a first pass and stored in the file header (1 KB overhead).

The residuals' 16-bit values are **split into high byte and low byte**, each coded with a separate rANS stream. This is because the high byte (which captures the magnitude of the residual) has a very different distribution from the low byte (which is more uniform for large residuals).

> **Reference (ANS theory):** Duda, J. (2009). *Asymmetric numeral systems.* arXiv:0902.0271.
> [https://arxiv.org/abs/0902.0271](https://arxiv.org/abs/0902.0271)

> **Reference (rANS in practice):** Giesen, F. (2014). *Interleaved entropy coders.* (ryg's blog, widely cited in codec implementations.)
> [https://arxiv.org/abs/1402.3392](https://arxiv.org/abs/1402.3392)

---

## File format

```
Offset  Size   Field
------  ----   -----
0       4      Magic bytes: "RAIS"
4       2      Version (
= 1), little-endian
6       4      Width  in pixels, little-endian
10      4      Height in pixels, little-endian
14      4      Number of streams (= 2), little-endian

  Stream 0 (high bytes of zigzag-mapped residuals):
    18      4      nsyms (= 256)
    22      1024   Frequency table: 256 × uint32_le (sum = 65536)
    1046    8      Compressed byte count N, uint64_le
    1054    N      rANS bitstream

  Stream 1 (low bytes):
    ...     4      nsyms (= 256)
    ...     1024   Frequency table
    ...     8      Compressed byte count M
    ...     M      rANS bitstream
```

Total header overhead: ~2 KB regardless of image size.

---

## Why not just use gzip/zstd on the raw file?

General-purpose compressors treat the file as a flat byte stream. They do not know that:
- Each value is a 16-bit big-endian integer
- Neighbouring values (spatially) are highly correlated
- The data is a 2-D image with smooth gradients and sparse bright stars

By exploiting this structure with a 2-D predictor first, we transform the data into a stream that is far easier to entropy-code. Empirically, this brings us close to the **bzip2 compression ratio at 4× the speed**.

---

## Results summary (8 training files, our machine)

| Compressor  | Total compressed | bits/byte | Compress (s) | Decompress (s) |
|-------------|-----------------|-----------|-------------|----------------|
| gzip-1      | 20.2 MB         | 4.48      | 0.54        | 0.25           |
| zstd-1      | 20.2 MB         | 4.49      | 0.14        | 0.06           |
| zstd-9      | 18.8 MB         | 4.18      | 0.83        | 0.09           |
| bzip2-9     | 16.1 MB         | 3.58      | 2.25        | 1.58           |
| **rais**    | **16.2 MB**     | **3.61**  | **0.55**    | **0.52**       |
