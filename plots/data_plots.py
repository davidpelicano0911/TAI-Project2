import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

DATA = '/home/alof/Desktop/TAI/TAI-Project2/data2/'
H, W = 1500, 1500

def load(name):
    raw = np.fromfile(DATA + name, dtype='>u2')
    return raw.reshape(H, W).astype(np.float32)

# Use image B (mixed sky) — smooth gradients, strong spatial correlation
img = load('B')

GREEN_BG  = '#f0faf4'
TEXT_DARK = '#222222'

# ── Fig 1: image crop (top-left 400×400) ─────────────────────────────────────
crop = img[:400, :400]

fig1, ax1 = plt.subplots(figsize=(4, 4))
fig1.patch.set_facecolor('white')
ax1.imshow(crop, cmap='gray', interpolation='nearest',
           vmin=np.percentile(crop, 1), vmax=np.percentile(crop, 99))
ax1.set_title('Image B', fontsize=11, color=TEXT_DARK, pad=8)
ax1.axis('off')
plt.tight_layout()
plt.savefig('/home/alof/Desktop/TAI/TAI-Project2/plots/data_crop.png',
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("crop saved")

# ── Fig 2: residual histogram using GAP predictor ────────────────────────────
def gap_pred(W, N, NW):
    dh = abs(int(W) - int(NW))
    dv = abs(int(N) - int(NW))
    if dv > 2 * dh:
        return int(W)
    elif dh > 2 * dv:
        return int(N)
    else:
        p = (int(W) * (dv + 1) + int(N) * (dh + 1)) / (dv + dh + 2)
        lo, hi = min(int(W), int(N)), max(int(W), int(N))
        return int(max(lo, min(hi, p)))

p  = img[1:, 1:].astype(np.int32)
W  = img[1:, :-1].astype(np.int32)
N  = img[:-1, 1:].astype(np.int32)
NW = img[:-1, :-1].astype(np.int32)

dh = np.abs(W - NW)
dv = np.abs(N - NW)

blend = (W * (dv + 1) + N * (dh + 1)).astype(np.float64) / (dv + dh + 2)
lo = np.minimum(W, N)
hi = np.maximum(W, N)
blend = np.clip(blend, lo, hi).astype(np.int32)

pred = np.where(dv > 2 * dh, W, np.where(dh > 2 * dv, N, blend))
diff = (p - pred).ravel()
pct = np.mean(np.abs(diff) <= 32) * 100
clip = 200

fig2, ax2 = plt.subplots(figsize=(5, 4))
fig2.patch.set_facecolor('white')
ax2.set_facecolor(GREEN_BG)

bins = np.arange(-clip - 0.5, clip + 1.5, 1)
ax2.hist(diff.clip(-clip, clip), bins=bins, color='#27ae60', edgecolor='none', zorder=3)

ax2.axvline(0, color='#e74c3c', linewidth=1.2, linestyle='--', zorder=5)

# compute and annotate 90th percentile
p90 = int(np.percentile(np.abs(diff), 90))
ax2.axvspan(-p90, p90, color='#27ae60', alpha=0.15, zorder=2)
ax2.axvline(-p90, color='#27ae60', linewidth=1, linestyle=':', zorder=4)
ax2.axvline( p90, color='#27ae60', linewidth=1, linestyle=':', zorder=4)
ax2.text(p90 + 4, 15000, f'90% of differences\nwithin ±{p90}', fontsize=8, color='#1a8a48', va='bottom', linespacing=1.4)

ax2.set_xlabel(r'$p - \hat{p}_{GAP}$', fontsize=10, color=TEXT_DARK)
ax2.set_ylabel('Count (log scale)', fontsize=10, color=TEXT_DARK)
ax2.set_title('GAP prediction residuals — image B', fontsize=10, color=TEXT_DARK, pad=8)
ax2.set_xlim(-clip, clip)
ax2.set_yscale('log')
ax2.tick_params(colors=TEXT_DARK)
for spine in ax2.spines.values():
    spine.set_edgecolor('#cccccc')
ax2.grid(True, axis='y', color='#d4edd9', linewidth=0.5, linestyle='--', zorder=0)

plt.tight_layout()
plt.savefig('/home/alof/Desktop/TAI/TAI-Project2/plots/data_scatter.png',
            dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("residual histogram saved")
