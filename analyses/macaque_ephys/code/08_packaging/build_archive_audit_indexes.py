#!/usr/bin/env python3
"""Build deterministic file, checksum, size and portability audit indexes."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

MODULE = Path(__file__).resolve().parents[2]
AUDIT = MODULE / "audit"
GENERATED = {
    "audit/FILE_INDEX.csv",
    "audit/SHA256SUMS.txt",
    "audit/PRECOMMIT_AUDIT.json",
    "audit/ABSOLUTE_PATH_SCAN.txt",
}
TEXT_SUFFIXES = {".py", ".r", ".R", ".md", ".txt", ".json", ".csv", ".tsv", ".html", ".js", ".css"}
ABSOLUTE = re.compile(r"(?:[A-Za-z]:[\\/](?:Users|Documents and Settings)[\\/]|C:[\\/]Users[\\/]|OneDrive[\\/])", re.I)
SECRET_PATTERNS = {
    "github_token": re.compile(r"(?:ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})"),
    "openai_key": re.compile(r"sk-[A-Za-z0-9_-]{20,}"),
    "declared_secret": re.compile(r"(?:api[_-]?key|password|access[_-]?token)\s*[:=]\s*['\"][^'\"]{8,}", re.I),
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def main() -> None:
    AUDIT.mkdir(exist_ok=True)
    files = sorted(
        p for p in MODULE.rglob("*")
        if p.is_file()
        and "__pycache__" not in p.parts
        and p.suffix.lower() != ".pyc"
        and p.relative_to(MODULE).as_posix() not in GENERATED
    )
    rows = []
    absolute_hits = []
    secret_hits = []
    categories = Counter()
    for path in files:
        rel = path.relative_to(MODULE).as_posix()
        size = path.stat().st_size
        rows.append((rel, size, digest(path)))
        categories[rel.split("/", 1)[0]] += 1
        if path.suffix in TEXT_SUFFIXES:
            try:
                text = path.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            for number, line in enumerate(text.splitlines(), 1):
                if ABSOLUTE.search(line):
                    absolute_hits.append(f"{rel}:{number}: {line.strip()[:300]}")
                for name, pattern in SECRET_PATTERNS.items():
                    if pattern.search(line):
                        secret_hits.append({"file": rel, "line": number, "pattern": name})

    with (AUDIT / "FILE_INDEX.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["relative_path", "bytes", "sha256"])
        writer.writerows(rows)
    (AUDIT / "SHA256SUMS.txt").write_text(
        "".join(f"{sha}  {rel}\n" for rel, _, sha in rows), encoding="utf-8"
    )
    (AUDIT / "ABSOLUTE_PATH_SCAN.txt").write_text(
        "Historical absolute-path references (not secrets)\n"
        "=================================================\n"
        + ("\n".join(absolute_hits) if absolute_hits else "None") + "\n",
        encoding="utf-8",
    )
    report = {
        "indexed_files_excluding_generated_indexes": len(rows),
        "indexed_bytes": sum(size for _, size, _ in rows),
        "counts_by_top_level_directory": dict(sorted(categories.items())),
        "files_over_25_MB": [
            {"path": rel, "bytes": size, "sha256": sha}
            for rel, size, sha in rows if size > 25 * 1024 * 1024
        ],
        "files_over_50_MB": [
            {"path": rel, "bytes": size, "sha256": sha}
            for rel, size, sha in rows if size > 50 * 1024 * 1024
        ],
        "absolute_path_occurrences": len(absolute_hits),
        "absolute_path_files": sorted({hit.split(":", 1)[0] for hit in absolute_hits}),
        "suspected_secret_hits": secret_hits,
        "note": "Generated index files are intentionally excluded from their own checksum index.",
    }
    (AUDIT / "PRECOMMIT_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
