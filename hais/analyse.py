#!/usr/bin/env python3
"""
Análise visual do efeito do preditor MED por ficheiro e por bloco.
Gera 3 gráficos:
  1. Entropia original vs resíduos MED por ficheiro
  2. Mapa de calor: ganho do MED bloco a bloco (ficheiro C vs ficheiro A)
  3. Distribuição de resíduos (histograma) para um bloco de fundo vs estrela
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

WIDTH, HEIGHT, BLOCK = 1500, 1500, 150
BLOCKS_X = BLOCKS_Y = WIDTH // BLOCK
DATA_DIR = Path("../data2")

# ---------------------------------------------------------------------------
# Funções auxiliares
# ---------------------------------------------------------------------------

def read_raw(path):
    return np.fromfile(path, dtype='>u2').reshape(HEIGHT, WIDTH).astype(np.int32)

def med_predict_block(block):
    """Aplica preditor MED a um bloco 2D, devolve resíduos."""
    pred = np.zeros_like(block)
    # primeira linha: preditor = vizinho esquerdo
    pred[0, 1:] = block[0, :-1]
    # primeira coluna: preditor = vizinho acima
    pred[1:, 0] = block[:-1, 0]
    # resto: MED
    A = block[1:, :-1]  # esquerda
    B = block[:-1, 1:]  # cima
    C = block[:-1, :-1] # cima-esquerda
    med = A + B - C
    lo  = np.minimum(A, B)
    hi  = np.maximum(A, B)
    pred[1:, 1:] = np.clip(med, lo, hi)
    return (block - pred).astype(np.int32)

def empirical_entropy(values):
    """Entropia empírica em bits/símbolo."""
    vals = values.flatten()
    _, counts = np.unique(vals, return_counts=True)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log2(p)))

def block_entropy_gain(image):
    """Calcula ganho de entropia do MED para cada bloco (positivo = MED ajuda)."""
    gain = np.zeros((BLOCKS_Y, BLOCKS_X))
    for by in range(BLOCKS_Y):
        for bx in range(BLOCKS_X):
            blk = image[by*BLOCK:(by+1)*BLOCK, bx*BLOCK:(bx+1)*BLOCK]
            H_raw = empirical_entropy(blk)
            H_med = empirical_entropy(med_predict_block(blk))
            gain[by, bx] = H_raw - H_med  # positivo = MED ganha
    return gain

# ---------------------------------------------------------------------------
# Gráfico 1 — Entropia original vs MED para os 8 ficheiros
# ---------------------------------------------------------------------------

files = list('ABCDEFGH')
H_raw_all, H_med_all = [], []

print("A calcular entropias...")
for f in files:
    path = DATA_DIR / f
    if not path.exists():
        print(f"  {f}: não encontrado, a saltar")
        H_raw_all.append(0); H_med_all.append(0)
        continue
    img = read_raw(path)
    residuals = med_predict_block(img)
    H_raw_all.append(empirical_entropy(img))
    H_med_all.append(empirical_entropy(residuals))
    print(f"  {f}: raw={H_raw_all[-1]:.3f}  med={H_med_all[-1]:.3f}  ganho={H_raw_all[-1]-H_med_all[-1]:+.3f}")

fig, axes = plt.subplots(1, 1, figsize=(10, 5))
x = np.arange(len(files))
w = 0.35
bars1 = axes.bar(x - w/2, H_raw_all, w, label='Pixels originais', color='steelblue')
bars2 = axes.bar(x + w/2, H_med_all, w, label='Resíduos MED',     color='tomato')
axes.set_xlabel('Ficheiro')
axes.set_ylabel('Entropia (bits/símbolo)')
axes.set_title('Entropia original vs resíduos MED — 8 imagens astronómicas')
axes.set_xticks(x); axes.set_xticklabels(files)
axes.legend()
axes.grid(axis='y', alpha=0.3)
# anotar ganho
for i in range(len(files)):
    g = H_raw_all[i] - H_med_all[i]
    axes.text(i, max(H_raw_all[i], H_med_all[i]) + 0.05,
              f'{g:+.2f}', ha='center', fontsize=8,
              color='green' if g > 0 else 'red')
plt.tight_layout()
plt.savefig('plot_entropy.png', dpi=150)
print("Guardado: plot_entropy.png")
plt.close()

# ---------------------------------------------------------------------------
# Gráfico 2 — Mapa de calor do ganho MED por bloco (ficheiro C e ficheiro A)
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, fname, title in zip(axes, ['C', 'A'],
                             ['Ficheiro C (flat bias — ruído puro)',
                              'Ficheiro A (alto brilho — estrutura)']):
    path = DATA_DIR / fname
    if not path.exists():
        ax.set_title(f'{title}\n(ficheiro não encontrado)')
        continue
    img = read_raw(path)
    gain = block_entropy_gain(img)
    vmax = max(abs(gain.min()), abs(gain.max()))
    im = ax.imshow(gain, cmap='RdYlGn', vmin=-vmax, vmax=vmax, origin='upper')
    plt.colorbar(im, ax=ax, label='Ganho MED (bits/símbolo)\nverde=MED ajuda  vermelho=MED piora')
    ax.set_title(title)
    ax.set_xlabel('Bloco X'); ax.set_ylabel('Bloco Y')
    # mostrar valores
    for by in range(BLOCKS_Y):
        for bx in range(BLOCKS_X):
            ax.text(bx, by, f'{gain[by,bx]:.1f}', ha='center', va='center',
                    fontsize=6, color='black')

plt.suptitle('Ganho do preditor MED por bloco 150×150', fontsize=13)
plt.tight_layout()
plt.savefig('plot_heatmap.png', dpi=150)
print("Guardado: plot_heatmap.png")
plt.close()

# ---------------------------------------------------------------------------
# Gráfico 3 — Histograma de resíduos: bloco de fundo vs bloco com estrela
# ---------------------------------------------------------------------------

path_b = DATA_DIR / 'B'
if path_b.exists():
    img_b = read_raw(path_b)
    # bloco de fundo (canto superior esquerdo)
    blk_bg  = img_b[0:BLOCK, 0:BLOCK]
    res_bg  = med_predict_block(blk_bg)
    # bloco com mais variância (provavelmente tem estrela)
    variances = np.array([[img_b[by*BLOCK:(by+1)*BLOCK, bx*BLOCK:(bx+1)*BLOCK].var()
                           for bx in range(BLOCKS_X)] for by in range(BLOCKS_Y)])
    by_s, bx_s = np.unravel_index(variances.argmax(), variances.shape)
    blk_star = img_b[by_s*BLOCK:(by_s+1)*BLOCK, bx_s*BLOCK:(bx_s+1)*BLOCK]
    res_star = med_predict_block(blk_star)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, blk, res, label in zip(
            axes,
            [blk_bg, blk_star],
            [res_bg, res_star],
            ['Bloco de fundo (sem estrela)', f'Bloco com estrela (bloco {by_s},{bx_s})']):
        ax.hist(blk.flatten(),  bins=80, alpha=0.6, label='Original',     color='steelblue')
        ax.hist(res.flatten(),  bins=80, alpha=0.6, label='Resíduos MED', color='tomato')
        H_r = empirical_entropy(blk)
        H_m = empirical_entropy(res)
        ax.set_title(f'{label}\nH_orig={H_r:.2f}  H_med={H_m:.2f}  ganho={H_r-H_m:+.2f}')
        ax.set_xlabel('Valor'); ax.set_ylabel('Frequência')
        ax.legend(); ax.set_yscale('log')

    plt.suptitle('Ficheiro B — distribuição original vs resíduos MED', fontsize=13)
    plt.tight_layout()
    plt.savefig('plot_histograms.png', dpi=150)
    print("Guardado: plot_histograms.png")
    plt.close()

print("\nPronto. Ficheiros gerados:")
print("  plot_entropy.png   — entropia por ficheiro")
print("  plot_heatmap.png   — mapa de calor por bloco")
print("  plot_histograms.png — distribuição de resíduos")
