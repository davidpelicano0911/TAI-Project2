# HAIS — Hybrid Astronomical Image Compressor

Compressor/descompressor lossless para imagens astronómicas raw de 16 bits em big-endian.

## Compilar

Executar `make` na diretoria `hais/`. Gera os binários `compress` e `decompress`.

## Utilização

`./compress <input.raw> <output.hais> [width height]`  
`./decompress <input.hais> <output.raw>`

Se `width` e `height` não forem indicados, o compressor assume `1500 x 1500` (tamanho do benchmark `data2`). Para imagens com outras dimensões, os argumentos são obrigatórios.

O nome dos ficheiros do paper segue a convenção `u16be-1x[rows]x[cols]`, pelo que se deve passar `width=cols` e `height=rows`. O compressor valida que `width × height × 2` corresponde ao tamanho real do ficheiro.

## Estratégia

A imagem é dividida em blocos de `150 × 150` pixels. Cada bloco é comprimido de forma independente e em paralelo (`std::thread`).

### Seleção de modo por bloco

Para cada bloco, o compressor testa 3 modos e escolhe o que minimiza o `byte_cost`:

| Modo | Preditor | Quando é escolhido |
|------|----------|--------------------|
| 0 — raw | nenhum | blocos com ruído puro, sem correlação espacial |
| 1 — MED | Median Edge Detector (LOCO-I) | blocos com arestas ou gradientes lineares |
| 2 — avg | média de 3 vizinhos causais | blocos planos com ruído residual de leitura |

O critério de seleção é `byte_cost = H(hi8) + H(lo8)`, que mede diretamente o custo dos dois streams que o encoder vai pagar, em vez de calcular entropia sobre os 65536 símbolos possíveis.

O modo raw vence quando qualquer preditor só introduz resíduos piores que os pixels originais. O avg vence em regiões planas porque a média de três vizinhos amortece o ruído, enquanto o MED pode amplificá-lo. O MED vence onde há estrutura espacial clara com arestas.

### Preditores

Os vizinhos causais são `A` (esquerda), `B` (cima) e `C` (diagonal cima-esquerda) — todos já reconstruídos no momento da previsão.

**MED** (Median Edge Detector): preditor não-linear que devolve `min(A,B)` se `C >= max(A,B)`, `max(A,B)` se `C <= min(A,B)`, e `A+B−C` caso contrário. Detecta arestas e extrapola gradientes lineares com precisão. O resultado está sempre em `[min(A,B), max(A,B)]`, sem necessidade de clamping.

**avg**: calcula `(A + B + C + 1) / 3`. Amortece o ruído de leitura em regiões planas, onde a extrapolação do MED pode amplificar pequenas flutuações.

O modo é gravado em 1 byte antes dos streams do bloco. O descompressor usa o mesmo preditor para reconstruir o pixel exato.

### Codificação entrópica: FSE/tANS

Cada símbolo de 16 bits é separado em dois bytes (`hi8 = sym >> 8`, `lo8 = sym & 0xFF`) e cada stream é comprimido independentemente com **FSE/tANS** (Finite State Entropy / tabled ANS), implementado de raiz sem bibliotecas externas.

A tabela FSE tem `SCALE = 2^16 = 65536` estados. As frequências dos 256 bytes são normalizadas para soma `SCALE`. O encoder processa os símbolos em ordem inversa e guarda os bits de transição; o payload começa com o estado final de 16 bits. O decoder reconstrói os símbolos em ordem normal a partir desse estado.

### Formato das tabelas de frequência

| Flag | Formato | Quando |
|------|---------|--------|
| `0x01..0xFE` | Sparse: `nnz × (sym:u8 + freq:u32le)` | 1 a 254 símbolos distintos |
| `0x00` | Dense: `256 × u32le` | 255 símbolos distintos |
| `0xFF` | Uniforme implícito | 256 símbolos distintos |

O caso uniforme evita guardar 1024 bytes de tabela quando a distribuição é praticamente plana — comum no byte baixo de imagens com ruído de leitura.

## Formato do ficheiro

Cabeçalho: magic `HAIS` (4B), width (4B), height (4B), block_size (4B). Para cada bloco em ordem row-major: mode (1B), stream hi8 (table_flag + table_data + stream_size:8B + bits), stream lo8 (mesmo formato). Todos os inteiros multibyte são little-endian.

## Paralelismo

A compressão é paralela: um contador atómico distribui índices de blocos pelas threads disponíveis (`hardware_concurrency()`). A escrita final é sequencial para manter a ordem dos blocos. A descompressão é sequencial.

## Diferenças face a abordagens existentes

O HAIS combina um conjunto de decisões de design que, individualmente, existem em outros sistemas, mas cuja combinação específica não existe em nenhuma solução publicada.

**Face ao JPEG-LS / LOCO-I:** o JPEG-LS usa o preditor MED com um codificador Golomb-Rice de parâmetro adaptativo. O HAIS usa o mesmo MED mas substitui o Golomb-Rice por FSE/tANS, que se adapta melhor a distribuições não-geométricas, e adiciona o modo avg como alternativa por bloco.

**Face ao Quintas-Torra 2026:** o paper usa preditores lineares ponderados (2×2, 9×9, ou treinados por dataset) com um codificador aritmético binário contextual (CCSDS 123.0-B-2). O HAIS usa preditores não-lineares (MED) e uma média simples (avg), sem fase de treino, com FSE/tANS. A seleção do preditor é feita por bloco em vez de ser global para toda a imagem.

**Face ao zstd / compressores genéricos:** os compressores genéricos não exploram a estrutura espacial 2D da imagem. O HAIS aplica preditores causais que removem redundância espacial antes da codificação entrópica, reduzindo a entropia dos resíduos face aos pixels brutos.

**O que é próprio desta implementação:**
- Seleção automática de modo por bloco com critério `byte_cost = H(hi8) + H(lo8)`, mais preciso que entropia sobre 16 bits
- Split hi8/lo8: os dois bytes do símbolo de 16 bits são comprimidos em streams independentes, cada um com a sua tabela FSE adaptada à sua distribuição
- FSE/tANS implementado de raiz em C++, sem bibliotecas externas
- Paralelismo por bloco com contador atómico (`std::thread`)
- Três modos competitivos (raw / MED / avg) onde nenhum domina globalmente — a escolha depende das características locais de cada bloco
