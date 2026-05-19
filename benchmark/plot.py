#!/usr/bin/env python3
"""
Plot compression benchmark results.

Usage:
    python plot.py [results.csv] [results_plot.png]

Produces a scatter plot of bits/byte vs total time (compress + decompress),
one point per compressor (averaged over all files).
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "results.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results_plot.png"

df = pd.read_csv(CSV)
df = df[df["lossless"] == True].copy()

df["t_total_s"] = (df["compress_ms"] + df["decompress_ms"]) / 1000.0

agg = (
    df.groupby("compressor")
    .agg(
        bits_per_byte=("bpp", "mean"),
        t_total_s=("t_total_s", "sum"),
    )
    .reset_index()
)

FAMILY_COLORS = {
    "gzip":  "#2563eb",
    "bzip2": "#ea580c",
    "xz":    "#7c3aed",
    "zstd":  "#16a34a",
    "rais":  "#e11d48",
    "astra": "#f59e0b",
    "hais":  "#0d9488",
    "helix": "#be185d",
    "prism": "#dc2626",
    "src-fast": "#4ade80",
    "src-ratio": "#22c55e",
    "src-balanced": "#15803d",
}

def family(name):
    for k in FAMILY_COLORS:
        if name.startswith(k):
            return k
    return "other"

agg["family"] = agg["compressor"].apply(family)
agg["color"]  = agg["family"].map(FAMILY_COLORS).fillna("#6b7280")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(13, 7))

ax.grid(True, color="#e5e7eb", linewidth=0.6, which="both", zorder=0)
ax.set_facecolor("white")
fig.patch.set_facecolor("white")
ax.set_axisbelow(True)

# Scatter + labels
# Place each label at a small fixed offset in axis-fraction space, then
# iteratively push overlapping labels apart — no connector lines drawn.
for _, row in agg.iterrows():
    ax.scatter(
        row["t_total_s"], row["bits_per_byte"],
        color=row["color"], s=130, zorder=3,
        edgecolors="#374151", linewidths=0.5,
    )

texts = []
for _, row in agg.iterrows():
    t = ax.text(
        row["t_total_s"], row["bits_per_byte"],
        "   " + row["compressor"],
        fontsize=8, color="#111827",
        va="center", ha="left",
    )
    texts.append(t)

# Simple repulsion loop — nudge overlapping labels apart (no arrows)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

def get_box(t):
    return t.get_window_extent(renderer=renderer)

MAX_ITER = 60
for _ in range(MAX_ITER):
    moved = False
    for i, ti in enumerate(texts):
        bi = get_box(ti)
        for j, tj in enumerate(texts):
            if i == j:
                continue
            bj = get_box(tj)
            # overlap in both axes?
            ox = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            oy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if ox > 0 and oy > 0:
                # push vertically
                dy = (oy + 1) / 2
                xi, yi = ti.get_position()
                xj, yj = tj.get_position()
                # convert dy pixels to data coords
                inv = ax.transData.inverted()
                p0 = inv.transform((0, 0))
                p1 = inv.transform((0, dy))
                delta = abs(p1[1] - p0[1])
                ti.set_position((xi, yi + delta))
                tj.set_position((xj, yj - delta))
                moved = True
    if not moved:
        break

# Legend
handles = [
    plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=col,
               markeredgecolor="#374151", markeredgewidth=0.5,
               markersize=9, label=fam, linewidth=0)
    for fam, col in FAMILY_COLORS.items()
]
ax.legend(handles=handles, fontsize=9, loc="upper right",
          framealpha=0.9, edgecolor="#d1d5db")

ax.set_xscale("log")
ax.set_xlabel("t_total (seconds, log scale)  [compress + decompress, summed over all files]", fontsize=10)
ax.set_ylabel("bits / byte of original  (total over all files)", fontsize=10)
ax.set_title("Compression Benchmark — bits/byte vs total time", fontsize=12, fontweight="bold")

plt.tight_layout()
plt.savefig(OUT, dpi=150, facecolor="white")
print(f"Saved: {OUT}")
plt.show()
