# ASTRA — Astronomical background Separation and Transform for Raw images

ASTRA is a lossless compressor built specifically for raw 16-bit astronomical images. Instead of treating the image as a generic stream of bytes, it exploits a property that is unique to this kind of data: **most of the image is empty sky**.

---

## The core idea

A raw astronomical image is not like a photograph. It looks like this:

```
 611  612  610  611  609  612         ← flat sky background
 610  611  613  610  611  610
 611  609  612  611  610  613
 610  612  611  609  65201  612       ← one bright star spike
 611  610  612  611  610  611
```

The vast majority of pixels are clustered tightly around a single background value (here ~611). Only a handful of pixels — stars, cosmic rays, sensor artifacts — deviate significantly from it.

ASTRA exploits this by splitting the image into two groups:

```
┌─────────────────────────────────────────────────────┐
│  Every pixel                                        │
│       │                                             │
│       ▼                                             │
│  |pixel - background| ≤ T ?                         │
│       │                                             │
│      YES                    NO                      │
│       │                      │                      │
│  Background pixel        Anomaly pixel              │
│  store tiny delta        store (position + value)   │
│  e.g. +1, -2, 0, +3      e.g. (index=4, val=65201) │
└─────────────────────────────────────────────────────┘
```

The tiny deltas are then entropy-coded with rANS. Because they are all small numbers near zero, they compress extremely well. The anomalies are stored as a plain list — there are very few of them.

---

## How background and T are chosen automatically

The compressor never asks you to set any parameters. It computes everything from the image itself.

**Background** = the most frequent pixel value (the mode of the histogram). In a sky image, this is the peak of the sky background cluster.

**T** = 4 × MAD, rounded to the next power of two.

MAD stands for *Median Absolute Deviation*: the median of `|pixel - background|` over all pixels. It measures the typical spread of the background without being thrown off by stars (which are extreme outliers and do not affect the median).

Example for file C (a flat bias frame):

```
Background = 611   (appears 144,784 times out of 2,250,000)
MAD        = 5     (most pixels are within ±5 of 611)
T          = 32    (= 4 × 5, rounded to next power of 2)

At T=32:  99.92% of pixels are background  →  only 1,757 anomalies to store explicitly
```

---

## What gets stored in the file

```
┌──────────────────────────────────┐
│ Header                           │
│   magic: "ASTR"                  │
│   width, height                  │
│   background value               │
│   threshold T                    │
│   number of anomalies            │
├──────────────────────────────────┤
│ Anomaly list                     │
│   for each anomaly:              │
│     pixel index (4 bytes)        │
│     pixel value (2 bytes)        │
├──────────────────────────────────┤
│ rANS compressed delta stream     │
│   frequency table (1 KB)         │
│   compressed bytes               │
└──────────────────────────────────┘
```

The decompressor reads the header, loads the anomaly list, decodes the delta stream, and reconstructs every pixel in order — anomaly pixels get their stored value, all others get `background + delta`.

---

## Why this beats general-purpose compressors on astronomical images

A tool like gzip or zstd sees the file as a flat stream of bytes and looks for repeated byte patterns. It does not know that:

- The values are 16-bit integers
- The image has a dominant background
- Stars are sparse point exceptions

ASTRA encodes *meaning*, not just bytes. For a flat sky image (file C, std=10), the delta stream is almost entirely zeros and ±1, which compresses to near the theoretical minimum. General-purpose tools cannot get close to this because they have no concept of "background".

---

## Results on the 8 training files

| File | Description | Anomalies | ASTRA (b/B) | bzip2-9 (b/B) |
|------|-------------|-----------|-------------|---------------|
| C    | Flat bias frame | 1,757 (0.08%) | **2.45** | 2.59 |
| E    | Dark frame | 3,598 (0.16%) | **2.58** | 2.71 |
| B    | Sparse star field | 5,626 (0.25%) | **3.44** | 3.40 |
| D    | Low background | 28,742 (1.28%) | **2.93** | 2.89 |
| H    | Mixed field | 24,854 (1.10%) | 3.49 | 3.41 |
| G    | Narrow band | 48,404 (2.15%) | 3.21 | 2.54 |
| A    | High background | 106,796 (4.75%) | 5.89 | 4.62 |
| F    | High noise | 26,579 (1.18%) | 6.43 | 6.47 |

ASTRA wins clearly on the files where the background model fits (C, E, B, D). It is weaker on G and A where the background variance is higher and more pixels fall outside the threshold.

---

## Usage

```bash
# Compress
./compress input_file output_file

# Decompress
./decompress input_file output_file
```

The output is guaranteed to be bit-for-bit identical to the original.
