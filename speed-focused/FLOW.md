# PRISM-Block Flow Diagram

```mermaid
flowchart TD
    A([Raw 16-bit image file]) --> B[Read & byte-swap to uint16 array]
    B --> C[Split into horizontal strips\nfull width × 150 rows each]

    C --> D[Thread pool\none band per thread]

    subgraph ENCODE ["encode_block  —  runs independently per band"]
        direction TB
        E[Init models\n1× hi AdaptModel\n8× lo AdaptModel] --> F

        F{First row\nof band?} -- yes --> G[pred = 0\nfor pixel 0\npred = left pixel\nfor rest]
        F -- no --> H[pred = N\nfor pixel 0\npred = GAP W,N,NW\nfor rest]

        G --> I[Compute residual\nu = pixel − pred]
        H --> I

        I --> J[Zigzag encode u\nzz = u×2 if u≥0\nzz = −u×2−1 if u<0]
        J --> K[Split zz\nhi = zz >> 8\nlo = zz & 0xFF]

        K --> L[Compute lo context\nmW = abs prev_res left\nmN = abs prev_res above\nctx = LO_CTX_TAB mW+mN]

        L --> M[Range-encode hi\nusing hi AdaptModel\nFenwick prefix_less]
        L --> N[Range-encode lo\nusing lo AdaptModel ctx\nFenwick prefix_less]

        M --> O[Update hi model\nrescale if total ≥ 65536]
        N --> P[Update lo model ctx\nrescale if total ≥ 65536]

        O --> Q[Store residual\ncurr_res gx = res]
        P --> Q
        Q --> R{More pixels\nor rows?}
        R -- yes --> F
        R -- no --> S[enc.finish\nflush 5 carry bytes]
    end

    D --> ENCODE
    ENCODE --> T[Block payload bytes]

    T --> U[Write .prismb file\nMAGIC + W + H + block_rows + nblocks\nthen per block: row0 + rows + size + payload]

    U --> V([Compressed .prismb file])

    style ENCODE fill:#f0f4ff,stroke:#4466cc,stroke-width:1.5px
    style A fill:#e8f5e9,stroke:#2e7d32
    style V fill:#fff3e0,stroke:#e65100
    style D fill:#fce4ec,stroke:#c62828
```

## Decoder Flow

```mermaid
flowchart TD
    A([Compressed .prismb file]) --> B[Read header\nW × H × block_rows × nblocks]
    B --> C[Read all block payloads into memory]
    C --> D[Thread pool\none block per thread]

    subgraph DECODE ["decode_block  —  runs independently per band"]
        direction TB
        E[Init models\n1× hi + 8× lo AdaptModel] --> F
        F[Init RangeDecoder\nload 5 seed bytes] --> G

        G{First row\nof band?} -- yes --> H[pred = 0 pixel 0\npred = left pixel rest]
        G -- no --> I[pred = prev_img 0 pixel 0\npred = GAP W N NW rest]

        H --> J[Compute lo context\nmW + mN from residual history]
        I --> J

        J --> K[Decode hi byte\nfrom hi AdaptModel]
        K --> L[Decode lo byte\nfrom lo AdaptModel ctx]

        L --> M[Reconstruct zz\nzz = hi<<8 | lo]
        M --> N[Zigzag decode\nu = zz/2 if even\nu = 65536 − zz+1 /2 if odd]
        N --> O[Recover pixel\npixel = pred + u]

        O --> P[Update both models]
        P --> Q{More pixels\nor rows?}
        Q -- yes --> G
        Q -- no --> R[Write row big-endian\nto output buffer at row offset]
    end

    D --> DECODE
    DECODE --> S[Assembled output buffer\nwidth × height × 2 bytes]
    S --> T([Decoded raw file])

    style DECODE fill:#f0f4ff,stroke:#4466cc,stroke-width:1.5px
    style A fill:#fff3e0,stroke:#e65100
    style T fill:#e8f5e9,stroke:#2e7d32
    style D fill:#fce4ec,stroke:#c62828
```

## GAP Predictor Detail

```mermaid
flowchart TD
    A[Given W, N, NW] --> B[dh = abs W − NW\ndv = abs N − NW]
    B --> C{dh > 2×dv?}
    C -- yes --> D[pred = N\nstrong horizontal gradient]
    C -- no --> E{dv > 2×dh?}
    E -- yes --> F[pred = W\nstrong vertical gradient]
    E -- no --> G[wW = dv+1\nwN = dh+1\npred = wW×W + wN×N / wW+wN]
    G --> H[Clamp pred to min W,N .. max W,N]
    D --> I([Return pred])
    F --> I
    H --> I
```

## Adaptive Model & Range Coder Detail

```mermaid
flowchart LR
    subgraph MODEL ["AdaptModel  Fenwick tree, 256 symbols"]
        A[freq 256 entries\nall init to 1\ntotal = 256] --> B[Fenwick tree\nO log 256 prefix sums]
        B --> C[prefix_less sym\ncumulative count below sym]
        C --> D[update sym\nfreq sym ++\nif total ≥ 65536 → halve all freqs]
    end

    subgraph RANGE ["RangeEncoder  LZMA carry-propagation"]
        E[low = 0\nrange = 0xFFFFFFFF] --> F[encode sym\nr = range / total\nlow += cumul × r\nrange = freq × r]
        F --> G{range < 2²⁴?}
        G -- yes --> H[shift: emit carry byte\nrange <<= 8]
        H --> G
        G -- no --> I[next symbol]
        I --> F
    end

    MODEL --> RANGE
```
