#!/usr/bin/env python3
"""Write a deterministic byte-size and SHA-256 manifest for the module."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "MANIFEST.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    paths = sorted(p for p in ROOT.rglob("*") if p.is_file() and p != OUT)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["module_relative_path", "bytes", "sha256"])
        for path in paths:
            writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, sha256(path)])
    print(f"wrote {len(paths)} entries to {OUT}")


if __name__ == "__main__":
    main()
