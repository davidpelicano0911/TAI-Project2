#!/usr/bin/env python3
"""
View or save a raw 16-bit big-endian astronomical image.

Usage:
    python view.py <file>          # show interactively
    python view.py <file> out.png  # save to PNG instead
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

WIDTH, HEIGHT = 1500, 1500

path = sys.argv[1]
out  = sys.argv[2] if len(sys.argv) > 2 else None

# Read raw big-endian uint16
data = np.fromfile(path, dtype='>u2').reshape((HEIGHT, WIDTH))

# Stretch: use 1st–99th percentile so stars don't wash out the background
lo, hi = np.percentile(data, 1), np.percentile(data, 99)
stretched = np.clip((data.astype(np.float32) - lo) / (hi - lo), 0, 1)

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(stretched, cmap='gray', origin='upper')
ax.set_title(path)
ax.axis('off')
plt.tight_layout()

if out:
    plt.savefig(out, dpi=150)
    print(f"Saved: {out}")
else:
    plt.show()
