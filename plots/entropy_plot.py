import matplotlib.pyplot as plt
import numpy as np

images = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
h0     = [10.27, 6.81, 4.87, 5.43, 5.11, 12.52, 5.57, 6.61]
types  = ['Stellar field', 'Mixed sky', 'Dark frame', 'Bias frame',
          'Flat field', 'Science\n(sources)', 'Bias frame', 'Calibration']

type_colors = {
    'Bias frame':      '#27ae60',
    'Flat field':      '#2ecc71',
    'Dark frame':      '#a8d5b5',
    'Calibration':     '#a8d5b5',
    'Mixed sky':       '#f39c12',
    'Stellar field':   '#e67e22',
    'Science\n(sources)': '#e74c3c',
}

colors = [type_colors[t] for t in types]

GREEN_BG  = '#f0faf4'
TEXT_DARK = '#222222'
TEXT_GRAY = '#666666'

fig, ax = plt.subplots(figsize=(8, 4.2))
fig.patch.set_facecolor('white')
ax.set_facecolor(GREEN_BG)

bars = ax.bar(images, h0, color=colors, edgecolor='white', linewidth=0.8, zorder=3)

for bar, val in zip(bars, h0):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
            f'{val:.2f}', ha='center', va='bottom', fontsize=8, color=TEXT_DARK)

for i, (bar, t) in enumerate(zip(bars, types)):
    ax.text(bar.get_x() + bar.get_width() / 2, 0.2,
            t, ha='center', va='bottom', fontsize=6.5, color='white',
            fontweight='bold', linespacing=1.3)

ax.set_xlabel('Image', fontsize=10, color=TEXT_DARK)
ax.set_ylabel('Zero-order entropy $H_0$ (bpp)', fontsize=10, color=TEXT_DARK)
ax.set_title('Benchmark dataset entropy per image', fontsize=12, color=TEXT_DARK, pad=10)
ax.set_ylim(0, 14)
ax.tick_params(colors=TEXT_DARK)
for spine in ax.spines.values():
    spine.set_edgecolor('#cccccc')
ax.grid(True, axis='y', color='#d4edd9', linewidth=0.7, linestyle='--', zorder=0)

# Legend
from matplotlib.patches import Patch
legend_items = [
    Patch(color='#27ae60',   label='Bias frame'),
    Patch(color='#2ecc71',   label='Flat field'),
    Patch(color='#a8d5b5',   label='Dark / Calibration'),
    Patch(color='#f39c12',   label='Mixed sky'),
    Patch(color='#e67e22',   label='Stellar field'),
    Patch(color='#e74c3c',   label='Science (sources)'),
]
ax.legend(handles=legend_items, fontsize=7.5, facecolor='white',
          edgecolor='#cccccc', labelcolor=TEXT_DARK, loc='upper left')

plt.tight_layout()
plt.savefig('/home/alof/Desktop/TAI/TAI-Project2/plots/entropy_plot.png',
            dpi=150, bbox_inches='tight', facecolor='white')
print("saved")
