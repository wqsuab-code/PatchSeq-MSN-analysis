#!/usr/bin/env python3
"""Build strict D1/D2/ambiguous composition from frozen consensus labels."""

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data/04_frozen_classification/e_strict_t_identity.csv"
OUTPUT = ROOT / "results/04_T_E_RRR/strict_T_identity_by_Etype.csv"

rows = list(csv.DictReader(INPUT.open(encoding="utf-8-sig", newline="")))
if len(rows) != 450:
    raise RuntimeError(f"Expected 450 GC-HC consensus cells, found {len(rows)}")

counts = Counter()
for row in rows:
    stable = row["major_class_agreement"].strip().lower() == "true"
    label = row["T_identity"] if stable and row["T_identity"] in {"D1", "D2"} else "Ambiguous"
    counts[(row["E_class"], label)] += 1

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["E_class", "T_identity", "N"], lineterminator="\n")
    writer.writeheader()
    for e_class in ["E1", "E2", "E3", "E4", "E5"]:
        for label in ["D1", "D2", "Ambiguous"]:
            writer.writerow({"E_class": e_class, "T_identity": label, "N": counts[(e_class, label)]})

print(f"Wrote {OUTPUT}")
