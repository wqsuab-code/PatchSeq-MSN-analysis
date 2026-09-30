#!/usr/bin/env python3
"""Validate the archived GC15 merge, HC-GC concordance and frozen E labels."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

MODULE = Path(__file__).resolve().parents[2]
PARAMS = json.loads((MODULE / "config/analysis_parameters.json").read_text(encoding="utf-8"))


def main() -> None:
    root = MODULE / "data/04_frozen_classification"
    cells = pd.read_csv(root / "frozen_cell_assignments_and_tsne_n390.csv")
    mapping = pd.read_csv(root / "00_res3_GC15_to_C4_mapping.csv")
    mapped = cells["GC_raw"].map(mapping.set_index("raw_GC")["merged_class"])
    if not mapped.equals(cells["GC_merged"]):
        raise SystemExit("GC15-to-GC4 mapping does not reproduce archived merged labels")
    recomputed = cells["HC_class"].eq(cells["GC_merged"])
    if not recomputed.equals(cells["Consensus"]):
        raise SystemExit("Consensus indicator is inconsistent with HC and merged GC labels")

    consensus = cells.loc[cells.Consensus].copy()
    consensus["E_class"] = consensus.HC_class.map(PARAMS["class_label_map"])
    class_counts = consensus.E_class.value_counts().sort_index().to_dict()
    t_counts = consensus.T_class.value_counts().to_dict()
    report = {
        "complete_case_cells": int(len(cells)),
        "consensus_cells": int(len(consensus)),
        "consensus_fraction": float(len(consensus) / len(cells)),
        "class_counts": class_counts,
        "T_counts": t_counts,
    }
    print(json.dumps(report, indent=2))
    if class_counts != PARAMS["class_counts"]:
        raise SystemExit("Frozen E-class counts do not match parameters")
    if t_counts != PARAMS["T_counts_consensus"]:
        raise SystemExit("Consensus transcriptomic-class counts do not match parameters")


if __name__ == "__main__":
    main()
