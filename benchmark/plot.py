#!/usr/bin/env python3
"""
Plot compression benchmark results — light theme, Pareto frontier curve.
Usage: python plot.py [results.csv] [results_plot.png]
"""

import sys, csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
try:
    from adjustText import adjust_text
    HAS_ADJUST = True
except ImportError:
    HAS_ADJUST = False

CSV = sys.argv[1] if len(sys.argv) > 1 else "results.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results_plot.png"

# ---------------------------------------------------------------------------
# Aggregate: mean bpp, total time per compressor
# ---------------------------------------------------------------------------
groups = {}
with open(CSV, newline="") as f:
    for row in csv.DictReader(f):
        if row["lossless"].lower() != "true":
            continue
        name = row["compressor"]
        g = groups.setdefault(name, {"bpp": [], "t_s": 0.0})
        g["bpp"].append(float(row["bpp"]))
        g["t_s"] += (float(row["compress_ms"]) + float(row["decompress_ms"])) / 1000.0

agg = [{"name": n, "bpp": float(np.mean(v["bpp"])), "t_s": v["t_s"]}
       for n, v in groups.items()]

COLORS = {
    "gzip":        "#2563eb",
    "bzip2":       "#ea580c",
    "xz":          "#7c3aed",
    "zstd":        "#16a34a",
    "rais":        "#dc2626",
    "astra":       "#d97706",
    "hais":        "#0d9488",
    "helix":       "#be185d",
    "prism-block": "#6366f1",
    "prism":       "#9ca3af",
}

def family(n):
    for k in COLORS:
        if n.startswith(k):
            return k
    return "other"

for r in agg:
    r["family"] = family(r["name"])
    r["color"]  = COLORS.get(r["family"], "#6b7280")

# ---------------------------------------------------------------------------
# Pareto frontier
# ---------------------------------------------------------------------------
pts = np.array([[r["t_s"], r["bpp"]] for r in agg])
idx = np.argsort(pts[:, 0])
pareto = [pts[idx[0]]]
for p in pts[idx[1:]]:
    if p[1] < pareto[-1][1]:
        pareto.append(p)
pareto = np.array(pareto)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
BG, GRID, EDGE, LABEL = "white", "#e5e7eb", "#9ca3af", "#374151"

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": BG,
    "text.color": "#111827", "axes.labelcolor": LABEL,
    "xtick.color": "#6b7280", "ytick.color": "#6b7280",
    "axes.edgecolor": EDGE, "grid.color": GRID, "grid.linewidth": 0.5,
})

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.grid(True, which="both", zorder=0, alpha=0.7)
ax.set_axisbelow(True)

if len(pareto) >= 2:
    ax.plot(pareto[:, 0], pareto[:, 1], color="#9ca3af",
            linewidth=1.4, alpha=0.7, zorder=2, linestyle="--")

texts = []
for r in agg:
    ax.scatter(r["t_s"], r["bpp"], color=r["color"], s=90, zorder=4,
               edgecolors="#374151", linewidths=0.4)
    label = f"{r['name']}\n{r['bpp']:.2f} b/B"
    texts.append(ax.text(r["t_s"], r["bpp"], label,
                         fontsize=6.5, color="#111827", zorder=5,
                         va="bottom", ha="left"))

if HAS_ADJUST:
    try:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#9ca3af", lw=0.5),
                    expand_points=(1.4, 1.6), force_text=(0.4, 0.6))
    except Exception:
        pass

# Legend
handles = [plt.Line2D([0],[0], marker="o", color="w",
                      markerfacecolor=c, markeredgecolor="#374151",
                      markeredgewidth=0.4, markersize=7,
                      label=f, linewidth=0)
           for f, c in COLORS.items()]
ax.legend(handles=handles, fontsize=7.5, loc="upper right",
          framealpha=0.9, edgecolor=EDGE, facecolor="white")

ax.set_xscale("log")
ax.set_xlabel("t total (seconds)", fontsize=9)
ax.set_ylabel("bits / byte", fontsize=9)

# X axis: major ticks at every 0.1, 0.2, 0.5, 1, 2, 5, 10, 20 + minor grid
ax.xaxis.set_major_locator(ticker.LogLocator(base=10, subs=[1,2,3,4,5,6,7,8,9], numticks=20))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(
    lambda v, _: f"{v:g}" if v in (0.1,0.2,0.3,0.5,1,2,3,5,10,20,30) else ""))
ax.xaxis.set_minor_locator(ticker.NullLocator())

# Y axis: fine-grained ticks every 0.1 b/B, minor every 0.05
all_bpp = [r["bpp"] for r in agg]
ymin = max(0, min(all_bpp) - 0.2)
ymax = max(all_bpp) + 0.2
ax.set_ylim(ymin, ymax)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
ax.grid(True, which="minor", color=GRID, linewidth=0.3, alpha=0.5)

ax.set_title("Compression Benchmark — Overall", fontsize=11,
             fontweight="bold", pad=10)

plt.tight_layout()
plt.savefig(OUT, dpi=150, bbox_inches="tight", facecolor=BG)
print(f"Saved: {OUT}")
