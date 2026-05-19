# HAIS — Hybrid Astronomical Image Compressor

Compressor/descompressor lossless para imagens astronómicas raw de 16 bits em big-endian.

## Compilar

Executar `make` na diretoria `hais/`. Gera os binários `compress` e `decompress`.

## Utilização

```
./compress   <input.raw> <output.hais> [width height]
./decompress <input.hais> <output.raw>
```

Se `width` e `height` não forem indicados, o compressor assume `1500 × 1500` (tamanho do benchmark `data2`). Para imagens com outras dimensões, os argumentos são obrigatórios.

O compressor valida que `width × height × 2` corresponde ao tamanho real do ficheiro.

## Estratégia

A imagem é dividida em blocos de **500 × 500** pixels. Cada bloco é comprimido de forma independente e em paralelo (`std::thread`). A descompressão é também paralela — o ficheiro é lido integralmente para memória, os offsets dos blocos são determinados numa passagem sequencial, e os blocos são descomprimidos em paralelo.

### Seleção de modo por bloco

Para cada bloco, o compressor testa 5 modos e escolhe o que minimiza o `byte_cost`:

| Modo | Preditor | Overhead por bloco | Quando é escolhido |
|------|----------|--------------------|--------------------|
| 0 — raw       | nenhum                   | 0 B    | blocos com ruído puro, sem correlação espacial |
| 1 — avg       | média de 3 vizinhos      | 0 B    | blocos planos com ruído residual de leitura |
| 2 — LS        | regressão linear (W,N,NW)| 12 B   | blocos com gradientes suaves e estrutura linear |
| 3 — mean      | média global da imagem   | 0 B    | imagens com fundo uniforme não-nulo |
| 4 — LS+bias   | LS com bias (W,N,NW,1)   | 16 B   | blocos com offset DC significativo |

O critério de seleção é `byte_cost = H(hi8) + H(lo8)`, que estima diretamente o custo dos dois streams que o encoder vai pagar. O overhead em bytes dos modos LS e LS+bias é contabilizado antes da comparação.

### Preditores

Os vizinhos causais são `A` (esquerda), `B` (cima) e `C` (diagonal cima-esquerda) — todos já reconstruídos no momento da previsão.

**avg**: calcula `(A + B + C + 1) / 3`. Amortece o ruído de leitura em regiões planas, onde preditores extrapolativos amplificam pequenas flutuações.

**LS** (Least Squares per block): resolve localmente `min_w Σ(pixel − w·[A,B,C])²` usando todos os pixels interiores do bloco. Os pesos ótimos `w` são guardados no bitstream (3 × float32 = 12 bytes). O decoder lê os pesos e aplica a mesma predição linear, reproduzindo os resíduos exatos.

**LS+bias**: igual ao LS mas com vetor de features `[A,B,C,1]`, adicionando um bias independente do contexto. Guarda 4 × float32 = 16 bytes. Útil quando o bloco tem um offset DC que não é capturado pelos vizinhos.

**mean**: usa a média global de todos os pixels da imagem como preditor constante. Elimina a maior parte da variância em imagens com fundo uniforme de intensidade não-nula (frequente em imagens astronómicas com pedestal de bias). A média global é calculada uma vez e guardada no cabeçalho (2 bytes).

O resíduo é codificado com mapeamento zigzag signed→unsigned: `v≥0 → 2v`, `v<0 → -2v-1`.

### Codificação entrópica: FSE/tANS

Cada símbolo de 16 bits é separado em dois bytes (`hi8 = sym >> 8`, `lo8 = sym & 0xFF`) e cada stream é comprimido independentemente com **FSE/tANS** (Finite State Entropy / tabled ANS), implementado de raiz sem bibliotecas externas.

A tabela FSE tem `SCALE = 2^12 = 4096` estados. As frequências dos 256 bytes são normalizadas para soma `SCALE`. O encoder processa os símbolos em ordem inversa e guarda os bits de transição; o payload começa com o estado final de 16 bits. O decoder reconstrói os símbolos em ordem normal a partir desse estado.

### Formato das tabelas de frequência

| Flag | Formato | Quando |
|------|---------|--------|
| `0x01..0xFE` | Sparse: `nnz × (sym:u8 + freq:u32le)` | 1 a 254 símbolos distintos |
| `0x00` | Dense: `256 × u32le` | 255 símbolos distintos |
| `0xFF` | Uniforme implícito | 256 símbolos distintos |

O caso uniforme evita guardar 1024 bytes de tabela quando a distribuição é praticamente plana — comum no byte baixo de imagens com ruído de leitura.

## Formato do ficheiro

```
Header (18 bytes):
  magic      : 4 B  — "HAIS"
  width      : 4 B  — little-endian u32
  height     : 4 B  — little-endian u32
  block_size : 4 B  — little-endian u32
  gmean      : 2 B  — little-endian u16, média global da imagem

Por bloco (row-major):
  mode       : 1 B
  [pesos]    : 12 B se mode==2, 16 B se mode==4, 0 B caso contrário
  stream hi8 : flag(1B) + table_data + size(8B) + bits
  stream lo8 : flag(1B) + table_data + size(8B) + bits
```

Todos os inteiros multibyte são little-endian.

## Paralelismo

A compressão é paralela: um contador atómico distribui índices de blocos pelas threads disponíveis (`hardware_concurrency()`). A descompressão é também paralela: o ficheiro é lido para memória, os offsets de cada bloco são determinados sequencialmente, e os blocos são descomprimidos em paralelo com o mesmo mecanismo. A escrita final é sequencial para manter a ordem dos blocos.

## Resultados (benchmark `data2`, 8 imagens 1500×1500 u16be)

| Imagem | HAIS (bpp) | Balanced (bpp) | zstd-19 (bpp) |
|--------|-----------|----------------|---------------|
| A      | 4.424     | 4.51           | —             |
| B      | 3.266     | 3.33           | —             |
| C      | 2.439     | 2.46           | —             |
| D      | 2.747     | 2.72           | —             |
| E      | 2.572     | 2.60           | —             |
| F      | 6.152     | 6.21           | —             |
| G      | 2.419     | 2.42           | —             |
| H      | 3.314     | 3.27           | —             |
| **Média** | **3.417** | **3.430**   | —             |

Métrica: bits por byte do ficheiro original não-comprimido (`compressed_bytes × 8 / original_bytes`).

## Diferenças face a abordagens existentes

**Face ao JPEG-LS / LOCO-I:** o JPEG-LS usa o preditor MED com codificador Golomb-Rice adaptativo. O HAIS substitui o MED pelo avg (mais robusto a pixels isolados de alta intensidade sobre fundo escuro) e o Golomb-Rice pelo FSE/tANS, mais eficiente quando a distribuição dos resíduos não é Laplaciana.

**Face ao Quintas-Torra 2026:** o paper usa preditores lineares ponderados treinados globalmente por dataset. O HAIS resolve o LS localmente por bloco sem fase de treino — os pesos são derivados do próprio bloco e guardados no bitstream. A seleção do preditor é feita por bloco com critério `byte_cost`, não globalmente.

**Face a compressores genéricos (zstd, xz):** não exploram a estrutura espacial 2D. O HAIS aplica preditores causais que removem redundância espacial antes da codificação entrópica.

**O que é próprio desta implementação:**
- Preditor global_mean: elimina a variância de fundo em imagens com pedestal não-nulo, sem overhead por bloco
- Seleção automática de 5 modos por bloco com critério `byte_cost = H(hi8) + H(lo8)`
- LS per-block sem treino: pesos ótimos locais derivados de cada bloco, guardados no bitstream
- LS+bias: acrescenta um termo de bias independente ao preditor linear
- Split hi8/lo8: streams independentes com tabelas FSE adaptadas à sua distribuição
- FSE/tANS implementado de raiz em C++, sem bibliotecas externas
- Blocos 500×500 (divisor de 1500): sem blocos parciais no benchmark padrão, maximizando os dados de treino FSE por bloco
