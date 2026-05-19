#!/usr/bin/env python3
"""
Plot compression benchmark results — dark theme, Pareto frontier curve.

Usage:
    python plot.py [results.csv]
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import PchipInterpolator

CSV = sys.argv[1] if len(sys.argv) > 1 else "results.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results_plot.png"

df = pd.read_csv(CSV)
df = df[df["lossless"] == True].copy()
df["t_total_s"] = (df["compress_ms"] + df["decompress_ms"]) / 1000.0

agg = (
    df.groupby("compressor")
    .agg(bits_per_byte=("bpp", "mean"), t_total_s=("t_total_s", "sum"))
    .reset_index()
)

FAMILY_COLORS = {
    "gzip":  "#2563eb",
    "bzip2": "#ea580c",
    "xz":    "#7c3aed",
    "zstd":  "#16a34a",
    "rais":  "#dc2626",
    "astra": "#d97706",
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
agg["color"]  = agg["family"].map(FAMILY_COLORS).fillna("#9ca3af")

# ---------------------------------------------------------------------------
# Pareto frontier (lower-left envelope: best bpp for each time budget)
# ---------------------------------------------------------------------------
pts = agg[["t_total_s", "bits_per_byte"]].values
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
# Pareto curve — smooth monotone spline in log-log space (PchipInterpolator)
# ---------------------------------------------------------------------------
if len(pareto) >= 2:
    x_all = agg["t_total_s"].values
    x_min = min(x_all) * 0.6
    x_max = max(x_all) * 1.4

    log_x_p   = np.log10(pareto[:, 0])
    log_bpp_p = np.log10(pareto[:, 1])

    # Pad boundary points so the curve extends across the full data range
    # Left pad: extrapolate slope; right pad: gentle slope continuing down
    slope_l = (log_bpp_p[1]  - log_bpp_p[0])  / (log_x_p[1]  - log_x_p[0])  if len(pareto) > 1 else -0.1
    slope_r = -0.04  # gentle right tail — no compressor beats the last pareto point significantly

    lx_l = np.log10(x_min);  lbpp_l = log_bpp_p[0]  + slope_l * (lx_l - log_x_p[0])
    lx_r = np.log10(x_max);  lbpp_r = log_bpp_p[-1] + slope_r * (lx_r - log_x_p[-1])

    log_x_ext   = np.concatenate([[lx_l],   log_x_p,   [lx_r]])
    log_bpp_ext = np.concatenate([[lbpp_l], log_bpp_p, [lbpp_r]])

    interp   = PchipInterpolator(log_x_ext, log_bpp_ext)
    x_smooth = np.logspace(np.log10(x_min), np.log10(x_max), 400)
    y_smooth = 10.0 ** interp(np.log10(x_smooth))
    ax.plot(x_smooth, y_smooth, color="#6b7280", linewidth=1.6, alpha=0.8, zorder=2)

# ---------------------------------------------------------------------------
# Scatter points
# ---------------------------------------------------------------------------
for _, row in agg.iterrows():
    ax.scatter(
        row["t_total_s"], row["bits_per_byte"],
        color=row["color"], s=170, zorder=4,
        edgecolors="#374151", linewidths=0.5,
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
out = CSV.replace(".csv", "_plot.png")
plt.savefig(out, dpi=150, bbox_inches="tight", facecolor=BG)
print(f"Saved: {out}")
plt.show()
