#!/usr/bin/env python3
"""Single entry point that regenerates every output/ deliverable.

The heavy physics was executed once with the locked CPython 3.12 / NumPy 2.5.2 /
SciPy 1.18.0 environment through run_reference.py, producing the two audited
frozen products reference_143K.json and reference_250K.json (shipped alongside
this file and reproduced from src/ + config_*.json). This entry point rebuilds
results.csv, figure2_subset.png, diagnostics.json, method.md and conclusions.md
deterministically from those frozen results.

    python3 run_all.py [--reference-dir .] [--output-dir .]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import build_outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--src-dir", type=Path, default=None)
    args = parser.parse_args()

    ref = args.reference_dir
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    build_outputs.generate(
        ref / "reference_143K.json",
        ref / "reference_250K.json",
        args.src_dir or (ref / "src"),
        out,
    )
    print("run_all.py: wrote all deliverables to", out)


if __name__ == "__main__":
    main()
