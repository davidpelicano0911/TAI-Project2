#!/usr/bin/env python3
"""
Plot compression benchmark results — dark theme, Pareto frontier curve.

Usage:
    python plot.py [results.csv] [results_plot.png]
"""

import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

CSV = sys.argv[1] if len(sys.argv) > 1 else "results.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results_plot.png"

groups = {}
with open(CSV, newline="") as f:
    for row in csv.DictReader(f):
        if row["lossless"].lower() != "true":
            continue
        name = row["compressor"]
        groups.setdefault(name, {"compressor": name, "bpp": [], "t_total_s": 0.0})
        groups[name]["bpp"].append(float(row["bpp"]))
        groups[name]["t_total_s"] += (float(row["compress_ms"]) + float(row["decompress_ms"])) / 1000.0

agg = []
for item in groups.values():
    agg.append({
        "compressor": item["compressor"],
        "bits_per_byte": float(np.mean(item["bpp"])),
        "t_total_s": item["t_total_s"],
    })

FAMILY_COLORS = {
    "gzip":  "#2563eb",
    "bzip2": "#ea580c",
    "xz":    "#7c3aed",
    "zstd":  "#16a34a",
    "rais":  "#dc2626",
    "astra": "#d97706",
    "hais":  "#0d9488",
    "helix": "#be185d",
    "prism-block": "#991b1b",
    "prism": "#dc2626",
}

def family(name):
    for k in FAMILY_COLORS:
        if name.startswith(k):
            return k
    return "other"

for row in agg:
    row["family"] = family(row["compressor"])
    row["color"] = FAMILY_COLORS.get(row["family"], "#9ca3af")

# ---------------------------------------------------------------------------
# Pareto frontier (lower-left envelope: best bpp for each time budget)
# ---------------------------------------------------------------------------
pts = np.array([[row["t_total_s"], row["bits_per_byte"]] for row in agg])
sorted_idx = np.argsort(pts[:, 0])
sorted_pts = pts[sorted_idx]
pareto = [sorted_pts[0]]
for p in sorted_pts[1:]:
    if p[1] < pareto[-1][1]:
        pareto.append(p)
pareto = np.array(pareto)

# ---------------------------------------------------------------------------
# Light theme
# ---------------------------------------------------------------------------
BG    = "white"
GRID  = "#e5e7eb"
EDGE  = "#9ca3af"
LABEL = "#374151"
TICK  = "#6b7280"

plt.rcParams.update({
    "figure.facecolor": BG,
    "axes.facecolor":   BG,
    "text.color":       "#111827",
    "axes.labelcolor":  LABEL,
    "xtick.color":      TICK,
    "ytick.color":      TICK,
    "axes.edgecolor":   EDGE,
    "grid.color":       GRID,
    "grid.linewidth":   0.6,
})

fig, ax = plt.subplots(figsize=(13, 7))
ax.grid(True, which="both", zorder=0)
ax.set_axisbelow(True)

# ---------------------------------------------------------------------------
# Pareto curve
# ---------------------------------------------------------------------------
if len(pareto) >= 2:
    ax.plot(pareto[:, 0], pareto[:, 1], color="#6b7280",
            linewidth=1.6, alpha=0.8, zorder=2)

# ---------------------------------------------------------------------------
# Scatter points
# ---------------------------------------------------------------------------
for row in agg:
    ax.scatter(
        row["t_total_s"], row["bits_per_byte"],
        color=row["color"], s=170, zorder=4,
        edgecolors="#374151", linewidths=0.5,
    )
    ax.annotate(
        row["compressor"],
        (row["t_total_s"], row["bits_per_byte"]),
        xytext=(7, 0),
        textcoords="offset points",
        va="center",
        ha="left",
        fontsize=8,
        color="#111827",
        zorder=5,
        clip_on=False,
    )

# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------
handles = [
    plt.Line2D([0], [0], marker="o", color="w",
               markerfacecolor=col, markeredgecolor="#374151",
               markeredgewidth=0.5, markersize=9,
               label=fam, linewidth=0)
    for fam, col in FAMILY_COLORS.items()
]
ax.legend(handles=handles, fontsize=9, loc="upper right",
          framealpha=0.9, edgecolor=EDGE, facecolor="white",
          labelcolor="#111827")

# ---------------------------------------------------------------------------
# Axes
# ---------------------------------------------------------------------------
ax.set_xscale("log")
ax.set_xlabel("t_total (seconds)", fontsize=10)
ax.set_ylabel("bits/byte", fontsize=10)

ax.xaxis.set_major_formatter(ticker.FuncFormatter(
    lambda v, _: f"{v:g}"
))

ax.set_title("Evaluation Plot", fontsize=14, fontweight="bold",
             color="#111827", pad=18)
ax.text(0.5, 1.028, "Showing results for Overall",
        transform=ax.transAxes, ha="center", fontsize=9, color=LABEL)

plt.tight_layout()
plt.savefig(OUT, dpi=150, bbox_inches="tight", facecolor=BG)
print(f"Saved: {OUT}")
plt.show()
