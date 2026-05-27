import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse

fig, ax = plt.subplots(figsize=(10, 4.5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.axis('off')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

GREEN_DARK = '#1a5c36'
GREEN_MID  = '#27ae60'
GREEN_BG   = '#f0faf4'
GREEN_ACC  = '#a8dfc0'
TEXT_DARK  = '#1a1a1a'
TEXT_MID   = '#555555'

BOX_W  = 0.18
BOX_H  = 0.13
OVAL_W = 0.13
OVAL_H = 0.09

# Layout:
#   Row 1 (y=0.82): Start(oval)  -->  Raw Image  -->  Spatial Prediction
#   Row 2 (y=0.42): Compressed Output(oval)  <--  Entropy Coding  <--  Hi/Lo Split  <--  Residual

nodes = {
    'raw':        ('box',     'Raw Image',           '1500×1500, 16-bit',          0.18, 0.68),
    'pred':       ('box_key', 'Spatial Prediction',  'GAP / avg / gmean / LS',     0.42, 0.68),
    'resid':      ('box',     'Residual',            r'$r = p - \hat{p}$, zigzag', 0.66, 0.68),
    'split':      ('box',     'Hi / Lo Split',       'hi = r>>8,  lo = r&0xFF',    0.66, 0.38),
    'entropy':    ('box',     'Entropy Coding',      'FSE / range coder',          0.42, 0.38),
    'output':     ('oval',    'Compressed\nOutput',  '',                           0.18, 0.38),
}

arrows = [
    ('raw',     'pred',    'right'),
    ('pred',    'resid',   'right'),
    ('resid',   'split',   'down'),
    ('split',   'entropy', 'left'),
    ('entropy', 'output',  'left'),
]

def draw_node(ax, kind, label, subtitle, cx, cy):
    if kind == 'oval':
        el = Ellipse((cx, cy), OVAL_W, OVAL_H,
                     facecolor=GREEN_ACC, edgecolor=GREEN_MID, linewidth=1.5, zorder=3)
        ax.add_patch(el)
        ax.text(cx, cy, label, ha='center', va='center',
                fontsize=8.5, fontweight='bold', color=GREEN_DARK, zorder=4, linespacing=1.3)
    elif kind == 'box_key':
        x0, y0 = cx - BOX_W/2, cy - BOX_H/2
        shadow = FancyBboxPatch((x0+0.004, y0-0.004), BOX_W, BOX_H,
                                boxstyle='round,pad=0.01,rounding_size=0.015',
                                linewidth=0, facecolor='#7ecba1', zorder=2)
        ax.add_patch(shadow)
        box = FancyBboxPatch((x0, y0), BOX_W, BOX_H,
                             boxstyle='round,pad=0.01,rounding_size=0.015',
                             linewidth=0, facecolor=GREEN_MID, zorder=3)
        ax.add_patch(box)
        ax.text(cx, cy + 0.026, label, ha='center', va='center',
                fontsize=9, fontweight='bold', color='white', zorder=4)
        ax.text(cx, cy - 0.026, subtitle, ha='center', va='center',
                fontsize=7, color='#d4f0e0', style='italic', zorder=4)
    else:
        x0, y0 = cx - BOX_W/2, cy - BOX_H/2
        box = FancyBboxPatch((x0, y0), BOX_W, BOX_H,
                             boxstyle='round,pad=0.01,rounding_size=0.015',
                             linewidth=1.3, edgecolor=GREEN_MID, facecolor=GREEN_BG,
                             zorder=3)
        ax.add_patch(box)
        ax.text(cx, cy + 0.026, label, ha='center', va='center',
                fontsize=9, fontweight='bold', color=GREEN_DARK, zorder=4)
        if subtitle:
            ax.text(cx, cy - 0.026, subtitle, ha='center', va='center',
                    fontsize=7, color=TEXT_MID, style='italic', zorder=4)

for key, (kind, label, subtitle, cx, cy) in nodes.items():
    draw_node(ax, kind, label, subtitle, cx, cy)

def node_edge(kind, cx, cy, side):
    hw = (OVAL_W if kind == 'oval' else BOX_W) / 2
    hh = (OVAL_H if kind == 'oval' else BOX_H) / 2
    if side == 'right': return (cx + hw, cy)
    if side == 'left':  return (cx - hw, cy)
    if side == 'top':   return (cx, cy + hh)
    if side == 'bottom':return (cx, cy - hh)

for (fk, tk, direction) in arrows:
    fkind, _, _, fcx, fcy = nodes[fk]
    tkind, _, _, tcx, tcy = nodes[tk]
    if direction == 'right':
        src = node_edge(fkind, fcx, fcy, 'right')
        dst = node_edge(tkind, tcx, tcy, 'left')
    elif direction == 'left':
        src = node_edge(fkind, fcx, fcy, 'left')
        dst = node_edge(tkind, tcx, tcy, 'right')
    elif direction == 'down':
        src = node_edge(fkind, fcx, fcy, 'bottom')
        dst = node_edge(tkind, tcx, tcy, 'top')

    ax.annotate('',
                xy=(dst[0], dst[1]),
                xytext=(src[0], src[1]),
                arrowprops=dict(arrowstyle='->', color=GREEN_MID,
                                lw=1.5, mutation_scale=13),
                zorder=2)

plt.tight_layout(pad=0.2)
plt.savefig('/home/alof/Desktop/TAI/TAI-Project2/plots/pipeline_diagram.png',
            dpi=150, bbox_inches='tight', facecolor='white')
print("saved")
