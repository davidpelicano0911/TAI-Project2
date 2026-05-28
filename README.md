# TAI Project #2 — Lossless Astronomical Image Compressor

**Course:** Algorithmic Information Theory (2025/26)  

## Authors

**Group:** 4


| Name | Email |
|------|-------|
| Afonso Ferreira | afonso.ferreira@ua.pt |
| David Pelicano | davidpoetapelicano@ua.pt |
| Tomás Bras | tomasbras@ua.pt |

## Project Overview

The challenge was to design and implement a **lossless compressor/decompressor** for raw astronomical images: 1500×1500 pixel, 16-bit unsigned, big-endian raster files. The goal was not only correctness but to outperform general-purpose tools (gzip, bzip2, xz, zstd) by exploiting the spatial structure of astronomical data.

Three distinct codecs were developed, each targeting a different trade-off between compression ratio and speed.

## Repository Structure

```
.
├── all-rounder/          # best overall codec (ratio + speed balance)
├── compression-focused/  # maximum compression ratio
├── speed-focused/        # fastest compression
├── submission/           # Submitted binaries (V3–V5, compress/decompress)
├── benchmark/            # Benchmarking harness and results
├── data2/                # Benchmark dataset (8 images A–H, 1500×1500 u16be)
├── docs/                 # Project specification and reference paper
├── Poster.pdf            # Project poster
└── Report.pdf            # Final report
```

### Codecs

#### `all-rounder/` — HAIS (Hybrid Astronomical Image Compressor)

The submitted final codec. Divides the image into blocks (500×500 for 1500×1500 images) compressed independently in parallel using `std::thread`. For each block, the compressor evaluates 5 prediction modes (raw, avg, global mean, least-squares with bias, least-squares with NN context) and picks the one with the lowest estimated entropy cost. Residuals are entropy-coded with a custom **FSE/tANS** implementation (no external libraries). Achieves an average of **3.40 bits/byte** on the benchmark dataset.

#### `compression-focused/` — APX2

Targets maximum compression at the cost of speed. Uses adaptive range coding (LZMA-style carry propagation) with block-parallel encoding (format magic `APX2`). Also includes a standalone range-coder variant (`range_compress`/`range_decompress`).

#### `speed-focused/` — PRISM-Block

Targets fast compression. Uses the GAP predictor, zigzag residuals, hi/lo byte split, 8 context-adaptive range-coded lo-byte streams, and block-parallel (`std::thread`) encoding. Format magic `PRB1`. Achieves compression in ~22 ms average on the benchmark dataset.


## Building

Each codec is built independently:

```bash
cd all-rounder && make        # produces compress, decompress
cd compression-focused && make
cd speed-focused && make
cd benchmark && make
```

Requires a C++17 compiler (g++ or clang++) on Linux.

## Usage

```bash
# Compress (dimensions default to 1500×1500 if omitted)
./compress <input.raw> <output.hais>
./compress <n_rows> <n_cols> <input.raw> <output.hais>

# Decompress
./decompress <input.hais> <output.raw>
```

The encoder defaults to 1500×1500 when no dimensions are provided. The decoder reads dimensions from the file header.

## Benchmark Results (data2, 8 images, 1500×1500 u16be)

Average bits/byte across all 8 images:

| Compressor | Avg bpp | Avg compress (ms) | Avg decompress (ms) |
|---|---|---|---|
| **all-rounder (HAIS)** | **3.40** | 58.7 | 44.2 |
| compression-focused (APX2) | 3.36 | 6155.9 | 133.2 |
| speed-focused (PRISM-Block) | 3.50 | 23.1 | 33.5 |
| gzip-9 | 4.59 | 723.5 | 29.0 |
| bzip2-9 | 3.58 | 260.3 | 141.0 |
| xz-9 | 3.68 | 1344.1 | 70.6 |
| zstd-19 | 4.01 | 1491.3 | 11.1 |

Metric: `compressed_bytes × 8 / original_bytes` (bits per byte of original file).

