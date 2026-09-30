#!/usr/bin/env python3
"""Build deterministic source and release-asset manifests for this module."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


MODULE = Path(__file__).resolve().parents[2]
SOURCE_ROOT = MODULE / "data" / "02_frozen_input"
OUT = MODULE / "data" / "00_source_manifest" / "source_manifest.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    rows = []
    for path in sorted(p for p in SOURCE_ROOT.rglob("*") if p.is_file()):
        rel = path.relative_to(MODULE).as_posix()
        category = path.relative_to(SOURCE_ROOT).parts[0]
        rows.append(
            {
                "module_relative_path": rel,
                "source_category": category,
                "source_release": "Mouse M frozen release 2026-09-23",
                "accession": "PENDING_DEPOSITION" if category != "ASC" else "NOT_APPLICABLE_LOCAL_ASC",
                "redistribution_status": "derived_or_mapping_only",
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} entries to {OUT}")


if __name__ == "__main__":
    main()
