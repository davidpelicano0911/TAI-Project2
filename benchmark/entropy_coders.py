#!/usr/bin/env python3
"""Benchmark RAIS entropy coder variants.

This compares the original rais/compress.cpp rANS executable with the
experimental rais/coder_compare variants that keep the same predictor and
residual representation while swapping only the entropy coder.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


WIDTH = 1500
HEIGHT = 1500
NPIX = WIDTH * HEIGHT
RAW_BYTES = NPIX * 2


@dataclass(frozen=True)
class Variant:
    name: str
    ext: str
    compress: tuple[str, ...]
    decompress: tuple[str, ...]


def run_timed(cmd: list[str]) -> tuple[int, float]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t1 = time.perf_counter()
    return proc.returncode, (t1 - t0) * 1000.0


def verify(a: Path, b: Path) -> bool:
    return subprocess.run(["cmp", "-s", str(a), str(b)]).returncode == 0


def resolve_input_path(raw: str, bench_dir: Path) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    cwd_path = path.resolve()
    if cwd_path.exists():
        return cwd_path
    return (bench_dir / path).resolve()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", nargs="?", default="../data2")
    parser.add_argument("--tmp-dir", default="/tmp/tai_entropy_coders")
    parser.add_argument("-o", "--csv", default="entropy_coder_results.csv")
    args = parser.parse_args()

    bench_dir = Path(__file__).resolve().parent
    repo = bench_dir.parent
    data_dir = resolve_input_path(args.data_dir, bench_dir)
    tmp_dir = Path(args.tmp_dir)
    csv_path = (bench_dir / args.csv).resolve()

    subprocess.check_call(["make", "-C", str(repo / "rais")], stdout=subprocess.DEVNULL)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    variants = [
        Variant(
            "compress.cpp-rans",
            ".rais",
            (str(repo / "rais" / "compress"), "{src}", "{dst}"),
            (str(repo / "rais" / "decompress"), "{src}", "{dst}"),
        ),
        Variant(
            "compare-rans",
            ".raic",
            (str(repo / "rais" / "coder_compare"), "c", "rans", "{src}", "{dst}"),
            (str(repo / "rais" / "coder_compare"), "d", "{src}", "{dst}"),
        ),
        Variant(
            "compare-huffman",
            ".raic",
            (str(repo / "rais" / "coder_compare"), "c", "huffman", "{src}", "{dst}"),
            (str(repo / "rais" / "coder_compare"), "d", "{src}", "{dst}"),
        ),
        Variant(
            "compare-arithmetic",
            ".raic",
            (str(repo / "rais" / "coder_compare"), "c", "arith", "{src}", "{dst}"),
            (str(repo / "rais" / "coder_compare"), "d", "{src}", "{dst}"),
        ),
    ]

    data_files = sorted(p for p in data_dir.iterdir() if p.is_file() and p.stat().st_size == RAW_BYTES)
    if not data_files:
        raise SystemExit(f"No {RAW_BYTES}-byte data files found in {data_dir}")

    rows: list[dict[str, object]] = []
    print("coder,file,comp_bytes,ratio,bits_per_pixel,compress_ms,decompress_ms,lossless")
    for variant in variants:
        for src in data_files:
            comp = tmp_dir / f"{src.name}{variant.ext}"
            dec = tmp_dir / f"{src.name}.dec"
            if comp.exists():
                comp.unlink()
            if dec.exists():
                dec.unlink()

            ccmd = [part.format(src=src, dst=comp) for part in variant.compress]
            dcmd = [part.format(src=comp, dst=dec) for part in variant.decompress]
            c_code, c_ms = run_timed(ccmd)
            d_code, d_ms = run_timed(dcmd)
            comp_bytes = comp.stat().st_size if comp.exists() else -1
            lossless = c_code == 0 and d_code == 0 and verify(src, dec)
            ratio = comp_bytes / RAW_BYTES if comp_bytes > 0 else -1.0
            bits_per_pixel = comp_bytes * 8.0 / NPIX if comp_bytes > 0 else -1.0

            row = {
                "coder": variant.name,
                "file": src.name,
                "orig_bytes": RAW_BYTES,
                "comp_bytes": comp_bytes,
                "ratio": ratio,
                "bits_per_pixel": bits_per_pixel,
                "compress_ms": c_ms,
                "decompress_ms": d_ms,
                "lossless": str(lossless).lower(),
            }
            rows.append(row)
            print(
                f"{variant.name},{src.name},{comp_bytes},{ratio:.6f},"
                f"{bits_per_pixel:.4f},{c_ms:.1f},{d_ms:.1f},{lossless}"
            )

            if comp.exists():
                comp.unlink()
            if dec.exists():
                dec.unlink()

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {csv_path}")
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
