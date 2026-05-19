# PRISM - Predictive Range-coded Image Stream

Lossless compressor for 1500x1500 raw 16-bit astronomical images. PRISM uses
a causal predictor, modular residual coding, and an inline adaptive range coder.
It does not depend on project-local coder headers.

## Algorithm

### 1. Gradient-Adaptive Predictor

PRISM predicts each pixel from three causal neighbours: W (left), N (above),
and NW (above-left).

```text
dh = |W - NW|
dv = |N - NW|

if dh > 2*dv   predict = N
elif dv > 2*dh predict = W
else           weighted blend of W and N, clamped to [min(W,N), max(W,N)]
```

The residual is computed modulo 16 bits and mapped with zigzag so small signed
residuals stay close to zero.

### 2. Hi/Lo Byte Split

Each 16-bit mapped residual is split into:

```text
hi = residual >> 8
lo = residual & 255
```

The high byte is coded with one adaptive model. The low byte is coded with one
of 8 adaptive models selected from the already reconstructed neighbouring
residual magnitudes.

### 3. Adaptive Range Coding

The entropy coder is implemented directly in `compress.cpp` and
`decompress.cpp`. Each model starts with frequency 1 for all 256 byte symbols,
updates after every symbol, and halves all frequencies when the total reaches
`1 << 16`.

The cumulative-frequency table is maintained with a Fenwick tree, so adaptive
updates and decoder symbol lookup stay logarithmic without changing the coded
probabilities. No histograms or rANS tables are stored in the file.

The implementation reads input in bulk, keeps residual/image neighbourhood
state in row buffers, and writes decompressed output a row at a time. These are
speed and memory-locality optimisations only; they do not change the stream.

## Build & Use

```bash
make
./compress   <input_raw>    <output.prism>
./decompress <input.prism>  <output_raw>
make test
```

## Format

```text
[4]  magic "PRM3"
[4]  width
[4]  height
[4]  n_lo_ctx = 8
[N]  interleaved adaptive range-coded stream:
       for each pixel in raster order:
         encode/decode hi byte with the high-byte model
         encode/decode lo byte with lo-byte model[ctx]
```

## Current Check

`make test` verifies all files in `../data2` by compressing, decompressing, and
comparing the decompressed output with the original raw input.
