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

A imagem é dividida em blocos de `150 × 150` pixels. Cada bloco é comprimido de forma independente e em paralelo (`std::thread`). A descompressão é também paralela — o ficheiro é lido integralmente para memória, os offsets dos blocos são determinados numa passagem sequencial, e os blocos são descomprimidos em paralelo.

### Seleção de modo por bloco

Para cada bloco, o compressor testa 3 modos e escolhe o que minimiza o `byte_cost`:

| Modo | Preditor | Quando é escolhido |
|------|----------|--------------------|
| 0 — raw | nenhum | blocos com ruído puro, sem correlação espacial |
| 1 — avg | média de 3 vizinhos causais | blocos planos com ruído residual de leitura |
| 2 — LS  | regressão linear por bloco | blocos com gradientes suaves e estrutura linear |

O critério de seleção é `byte_cost = H(hi8) + H(lo8)`, que mede diretamente o custo dos dois streams que o encoder vai pagar. O modo LS tem um overhead de 12 bytes (3 pesos float32) contabilizado no custo antes da comparação.

### Preditores

Os vizinhos causais são `A` (esquerda), `B` (cima) e `C` (diagonal cima-esquerda) — todos já reconstruídos no momento da previsão.

**avg**: calcula `(A + B + C + 1) / 3`. Amortece o ruído de leitura em regiões planas, onde preditores extrapolativos amplificam pequenas flutuações.

**LS** (Least Squares per block): resolve localmente `min_w Σ(pixel − w·[A, B, C])²` usando todos os pixels interiores do bloco. Os pesos ótimos `w` são guardados no bitstream (3 × float32 = 12 bytes por bloco). O decoder lê os pesos e aplica a mesma predição linear, reproduzindo os resíduos exatos. Não requer fase de treino — os pesos adaptam-se às características locais de cada bloco.

O modo é gravado em 1 byte antes dos streams do bloco.

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

Cabeçalho: magic `HAIS` (4B), width (4B), height (4B), block_size (4B). Para cada bloco em ordem row-major: mode (1B), e para o modo LS os 3 pesos float32 (12B), depois stream hi8 e stream lo8 (cada um com table_flag + table_data + stream_size:8B + bits). Todos os inteiros multibyte são little-endian.

## Paralelismo

A compressão é paralela: um contador atómico distribui índices de blocos pelas threads disponíveis (`hardware_concurrency()`). A descompressão é também paralela: o ficheiro é lido para memória, os offsets de cada bloco são determinados sequencialmente, e os blocos são descomprimidos em paralelo com o mesmo mecanismo. A escrita final é sequencial para manter a ordem dos blocos.

## Diferenças face a abordagens existentes

O HAIS combina um conjunto de decisões de design que, individualmente, existem em outros sistemas, mas cuja combinação específica não existe em nenhuma solução publicada.

**Face ao JPEG-LS / LOCO-I:** o JPEG-LS usa o preditor MED com um codificador Golomb-Rice de parâmetro adaptativo. O HAIS não usa MED — usa avg e LS — e substitui o Golomb-Rice por FSE/tANS.

**Face ao Quintas-Torra 2026:** o paper usa preditores lineares ponderados treinados globalmente por dataset (2×2, 9×9, ou otimizados). O HAIS resolve o LS localmente por bloco, sem fase de treino — os pesos são derivados do próprio bloco e guardados no bitstream. A seleção do preditor é feita por bloco com critério byte_cost, não globalmente.

**Face ao zstd / compressores genéricos:** os compressores genéricos não exploram a estrutura espacial 2D da imagem. O HAIS aplica preditores causais que removem redundância espacial antes da codificação entrópica.

**O que é próprio desta implementação:**
- Seleção automática de modo por bloco com critério `byte_cost = H(hi8) + H(lo8)`
- LS per-block sem treino: pesos ótimos locais derivados de cada bloco, guardados no bitstream
- Split hi8/lo8: streams independentes com tabelas FSE adaptadas à sua distribuição
- FSE/tANS implementado de raiz em C++, sem bibliotecas externas
- Paralelismo em compressão e descompressão com contador atómico (`std::thread`)
