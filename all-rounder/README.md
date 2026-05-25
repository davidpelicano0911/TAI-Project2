# HAIS — Hybrid Astronomical Image Compressor

Compressor/descompressor lossless para imagens astronómicas raw de 16 bits em big-endian.

## Compilar

```
cd hais/
make
```

Gera os binários `compress` (56 KB) e `decompress` (27 KB).

## Utilização

```
./compress   <n_rows> <n_cols> <input.raw> <output.hais>
./decompress <input.hais> <output.raw>
```

As dimensões são obrigatórias no encoder. O decoder lê as dimensões diretamente do cabeçalho do ficheiro `.hais`.

## Estratégia

A imagem é dividida em blocos. O tamanho de bloco é escolhido automaticamente: o maior valor ≤ 512 que divida exatamente `W` e `H` com pelo menos 4 blocos; se não existir, usa 256 (com blocos parciais). Para as imagens 1500×1500 do benchmark, o bloco é 500×500.

Cada bloco é comprimido de forma independente e em paralelo (`std::thread`). A descompressão é também paralela: o ficheiro é lido integralmente para memória, os offsets dos blocos são determinados numa passagem sequencial, e os blocos são descomprimidos em paralelo.

### Seleção de modo por bloco

Para cada bloco, o compressor avalia 5 modos e escolhe o que minimiza o custo estimado:

| Modo | Preditor | Overhead por bloco | Quando é escolhido |
|------|----------|--------------------|--------------------|
| 0 — raw      | nenhum                        | 0 B   | blocos com ruído puro, sem correlação espacial |
| 1 — avg      | média de W, N, NW             | 0 B   | blocos planos com ruído residual de leitura |
| 3 — mean     | média global da imagem        | 0 B   | blocos com fundo uniforme não-nulo (pedestal de bias) |
| 4 — LS+bias  | LS linear (W, N, NW, 1)       | 16 B  | blocos com gradiente suave e offset DC significativo |
| 6 — LS+NN    | LS linear (W, N, NW, NE, NN)  | 20 B  | blocos com correlação vertical forte (duas linhas de contexto) |

O critério de seleção é `byte_cost = H(hi8) + H(lo8)`, que estima o custo dos dois streams FSE. O overhead fixo em bytes dos modos LS é somado antes da comparação.

Adicionalmente, o compressor testa sempre uma variante **Context FSE**: o stream `lo8` é dividido em dois sub-streams separados consoante `hi8 == 0` ou `hi8 ≠ 0`. Se esta variante resultar num payload menor, é preferida e sinalizada com o bit 7 do byte de modo.

### Preditores

Os vizinhos causais são `W` (esquerda), `N` (cima), `NW` (diagonal), `NE` (cima-direita) e `NN` (dois acima) — todos já reconstruídos no momento da previsão.

**avg**: `(W + N + NW + 1) / 3`. Amortece o ruído de leitura em regiões planas.

**mean**: usa a média global de todos os pixels da imagem como preditor constante. A média é calculada uma vez e guardada no cabeçalho (2 bytes). Eficaz em imagens com pedestal de bias uniforme.

**LS+bias** (modo 4): resolve `min_w Σ(pixel − w·[W, N, NW, 1])²` por Gauss-Jordan usando os pixels interiores do bloco. Os 4 pesos float32 (16 bytes) são guardados no bitstream e lidos pelo decoder.

**LS+NN** (modo 6): igual mas com vetor `[W, N, NW, NE, NN]` — 5 pesos float32 = 20 bytes. Captura correlação vertical de segunda ordem (gradiente).

O resíduo é codificado com mapeamento zigzag: `v ≥ 0 → 2v`, `v < 0 → −2v − 1`.

### Codificação entrópica: FSE/tANS

Cada símbolo de 16 bits é separado em `hi8 = sym >> 8` e `lo8 = sym & 0xFF`. Cada stream é comprimido independentemente com **FSE/tANS** (Finite State Entropy / tabled ANS), implementado de raiz em C++ sem bibliotecas externas.

A tabela FSE usa `SCALE = 2^16 = 65536` estados. As frequências dos 256 bytes são normalizadas para soma `SCALE`. A fórmula de distribuição dos símbolos na tabela (`step = (SCALE>>1) + (SCALE>>3) + 3`) é proveniente do trabalho de Yann Collet sobre FSE [1].

O encoder processa os símbolos em ordem inversa e emite bits de transição; o decoder reconstrói os símbolos em ordem normal a partir do estado inicial (primeiros 2 bytes do stream).

### Formato das tabelas de frequência

| Flag | Formato | Quando |
|------|---------|--------|
| `0x01` | Símbolo único: 1 byte de símbolo, sem bitstream | apenas 1 símbolo distinto |
| `0x02..0xFE` | Sparse: `nnz × (sym:u8 + freq:u32le)` | 2 a 254 símbolos distintos |
| `0x00` | Dense: `256 × u32le` | 255 ou 256 símbolos distintos |

## Formato do ficheiro

```
Header (18 bytes):
  magic      : 4 B  — "HAIS"
  width      : 4 B  — u32 little-endian
  height     : 4 B  — u32 little-endian
  block_size : 4 B  — u32 little-endian
  gmean      : 2 B  — u16 little-endian, média global da imagem

Por bloco (row-major, left-to-right top-to-bottom):
  mode_byte  : 1 B  — bits [6:0] = modo (0/1/3/4/6), bit 7 = context FSE flag
  [pesos]    : 16 B se modo==4, 20 B se modo==6, 0 B caso contrário
  stream hi8 : flag(1B) + tabela + size(8B) + bits
  stream lo8 ou [lo_zero + lo_nonzero] se context FSE
```

Todos os inteiros multibyte são little-endian.

## Resultados (benchmark `data2`, 8 imagens 1500×1500 u16be)

| Imagem | HAIS | bzip2-9 | xz-9 | zstd-19 |
|--------|------|---------|------|---------|
| A | 4.39 | 4.62 | 4.93 | 5.46 |
| B | 3.24 | 3.40 | 3.55 | 3.77 |
| C | 2.43 | 2.59 | 2.61 | 2.74 |
| D | 2.73 | 2.89 | 2.96 | 3.16 |
| E | 2.56 | 2.71 | 2.76 | 2.91 |
| F | 6.14 | 6.47 | 6.56 | 7.06 |
| G | 2.41 | 2.54 | 2.56 | 2.71 |
| H | 3.29 | 3.41 | 3.56 | 4.20 |
| **Média** | **3.40** | **3.58** | **3.68** | **4.00** |

Métrica: `compressed_bytes × 8 / original_bytes` (bits por byte do ficheiro original).

## Referências

[1] Y. Collet, "Finite State Entropy — a new entropy coder," 2013.
    https://github.com/Cyan4973/FiniteStateEntropy
