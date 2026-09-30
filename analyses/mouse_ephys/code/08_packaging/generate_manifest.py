#!/usr/bin/env python3
"""Generate the module file index and SHA-256 list."""

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "audit"
SKIP = {"audit/FILE_INDEX.csv", "audit/SHA256SUMS.txt"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


files = [
    path for path in sorted(ROOT.rglob("*"))
    if (
        path.is_file()
        and path.relative_to(ROOT).as_posix() not in SKIP
        and "__pycache__" not in path.parts
        and path.suffix.lower() not in {".pyc", ".pyo"}
    )
]
rows = [
    {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }
    for path in files
]
AUDIT.mkdir(parents=True, exist_ok=True)
with (AUDIT / "FILE_INDEX.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
(AUDIT / "SHA256SUMS.txt").write_text(
    "".join(f"{row['sha256']}  {row['path']}\n" for row in rows),
    encoding="utf-8",
    newline="\n",
)
print(f"Indexed {len(rows)} files ({sum(row['bytes'] for row in rows)} bytes).")
