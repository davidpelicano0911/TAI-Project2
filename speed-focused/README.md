# PRISM Block

Block-parallel PRISM variant for 1500x1500 raw 16-bit astronomical images.

This compressor keeps PRISM's predictor, modular zigzag residuals, high/low
byte split, and 8 low-byte contexts, but divides the image into independent
horizontal row blocks. Each block has its own adaptive range-coded payload and
can be compressed or decompressed independently.

The design follows the same practical idea used by FITS tile compression:
independent image tiles trade a little model continuity for parallelism,
bounded memory, and robustness.

## Build & Use

```bash
make
./compress   <n_rows> <n_cols> <input_raw>     <output.prismb>
./decompress <input.prismb>  <output_raw>
make test
```

## Format

```text
[4]  magic "PRB1"
[4]  width
[4]  height
[4]  block_rows
[4]  n_blocks
for each block:
  [4]  row0
  [4]  rows
  [8]  payload byte count
  [N]  adaptive range-coded payload
```

Blocks are horizontal stripes. Inside each block, prediction and residual
contexts only use pixels/residuals already decoded inside the same block.
