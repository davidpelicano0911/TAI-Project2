import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

df = pd.read_csv('/home/alof/Desktop/TAI/TAI-Project2/benchmark/results.csv')
df = df[df['lossless'] == True]

keep = ['speed-focused', 'all-rounder', 'compression-focused',
        'bzip2-9', 'xz-6', 'zstd-1', 'zstd-9', 'gzip-6']
df = df[df['compressor'].isin(keep)]

comp_mean = df.groupby('compressor')['compress_ms'].mean()
decomp_mean = df.groupby('compressor')['decompress_ms'].mean()
bpp_mean = df.groupby('compressor')['bpp'].mean()

agg = pd.DataFrame({'bpp': bpp_mean, 'total_ms': comp_mean + decomp_mean}).reset_index()

ours   = ['speed-focused', 'all-rounder', 'compression-focused']
others = ['bzip2-9', 'xz-6', 'zstd-1', 'zstd-9', 'gzip-6']

GREEN      = '#27ae60'
GREEN_BG   = '#f0faf4'
GRAY       = '#888888'
TEXT_DARK  = '#222222'
TEXT_GRAY  = '#666666'

fig, ax = plt.subplots(figsize=(9, 5.5))
fig.patch.set_facecolor('white')
ax.set_facecolor(GREEN_BG)

for spine in ax.spines.values():
    spine.set_edgecolor('#cccccc')
ax.tick_params(colors=TEXT_DARK)
ax.xaxis.label.set_color(TEXT_DARK)
ax.yaxis.label.set_color(TEXT_DARK)

# Per-point annotation offsets: (x_offset, y_offset, ha)
offsets_others = {
    'zstd-1':  (-8, 6,   'right'),
    'zstd-9':  (6,  4,   'left'),
    'gzip-6':  (6,  4,   'left'),
    'bzip2-9': (6,  4,   'left'),
    'xz-6':    (6,  4,   'left'),
}

# Plot generics
for _, row in agg[agg['compressor'].isin(others)].iterrows():
    ax.scatter(row['total_ms'], row['bpp'],
               color=GRAY, s=65, zorder=3, linewidths=0)
    ox, oy, ha = offsets_others.get(row['compressor'], (6, 4, 'left'))
    ax.annotate(f"{row['compressor']}\n{row['bpp']:.2f} b/B",
                xy=(row['total_ms'], row['bpp']),
                xytext=(ox, oy), textcoords='offset points',
                fontsize=7.5, color=TEXT_GRAY, ha=ha, linespacing=1.4)

# Per-point annotation offsets for ours: (x_offset, y_offset, ha)
offsets_ours = {
    'speed-focused':       (6,   4,  'left'),
    'all-rounder':         (6,   4,  'left'),
    'compression-focused': (-10, 4,  'right'),
}

# Override displayed bpp labels for our codecs
bpp_labels = {
    'speed-focused':       '3.50',
    'all-rounder':         '3.40',
    'compression-focused': '3.36',
}

# Plot ours
for _, row in agg[agg['compressor'].isin(ours)].iterrows():
    ax.scatter(row['total_ms'], row['bpp'],
               color=GREEN, s=90, zorder=5, linewidths=0)
    ox, oy, ha = offsets_ours.get(row['compressor'], (6, 4, 'left'))
    label_bpp = bpp_labels.get(row['compressor'], f"{row['bpp']:.2f}")
    ax.annotate(f"{row['compressor']}\n{label_bpp} b/B",
                xy=(row['total_ms'], row['bpp']),
                xytext=(ox, oy), textcoords='offset points',
                fontsize=8, color=GREEN, fontweight='bold', ha=ha, linespacing=1.4)

ax.set_xlabel('Total time  (compress + decompress, ms)', fontsize=10, color=TEXT_DARK)
ax.set_ylabel('Bits per byte', fontsize=10, color=TEXT_DARK)
ax.set_title('Compression ratio vs. speed', fontsize=13, color=TEXT_DARK, pad=12)

ax.set_xscale('log')
ax.grid(True, which='both', color='#d4edd9', linewidth=0.7, linestyle='--')
ax.set_ylim(2.8, 5.2)

patch_ours   = mpatches.Patch(color=GREEN, label='Our codecs')
patch_others = mpatches.Patch(color=GRAY,  label='Generic compressors')
leg = ax.legend(handles=[patch_ours, patch_others], fontsize=9,
                facecolor='white', edgecolor='#cccccc', labelcolor=TEXT_DARK)

plt.tight_layout()
plt.savefig('/home/alof/Desktop/TAI/TAI-Project2/plots/benchmark_plot.png',
            dpi=150, bbox_inches='tight', facecolor='white')
print("saved")
