#!/usr/bin/env python
"""Freeze the reversible n=187 Morph discovery cohort after six-cell ASC review."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "outputs/morph_qc/morph_taxonomy_round6_blind_de_novo/01_blinded_193cell_morphology_input.csv"
CONFIG = ROOT / "config/morph_current_reconstruction_spe_local_pyramid_quarantine.json"
OUT = ROOT / "outputs/morph_qc/morph187_after_incomplete_reconstruction_quarantine"
ID = "MSN_unique_ID"
NEW_SIX = [
    "D251161_Batch8",
    "A20256247_Batch7",
    "A20257158_Batch7",
    "A2025797_Batch7",
    "A20256245_Batch7",
    "O9192_Batch2",
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dat = pd.read_csv(INPUT)
    if len(dat) != 193 or dat[ID].duplicated().any():
        raise RuntimeError("Expected the current 193-cell unique Morph cohort")
    missing = sorted(set(NEW_SIX) - set(dat[ID]))
    if missing:
        raise RuntimeError(f"New quarantine cells absent from n=193 input: {missing}")

    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    configured = set(config["suspected_reconstruction_failure"]["possible_under_reconstruction"])
    if not set(NEW_SIX).issubset(configured):
        raise RuntimeError("The six cells are not all frozen in the quarantine configuration")

    audit = dat[[ID]].copy()
    audit["Previous_primary_Morph_status"] = "included_n193"
    audit["Incomplete_reconstruction_quarantine"] = audit[ID].isin(NEW_SIX)
    audit["Current_status"] = audit["Incomplete_reconstruction_quarantine"].map(
        {True: "temporary_quarantine_incomplete_reconstruction", False: "included_n187"}
    )
    audit["Exclusion_scope"] = audit["Incomplete_reconstruction_quarantine"].map(
        {True: "Morph discovery only; reversible", False: "not excluded"}
    )
    audit.to_csv(OUT / "00_n193_to_n187_transition_audit.csv", index=False)

    retained = dat.loc[~dat[ID].isin(NEW_SIX)].copy()
    if len(retained) != 187:
        raise RuntimeError(f"Expected 187 retained cells, found {len(retained)}")
    retained.to_csv(OUT / "01_morph187_discovery_input.csv", index=False)
    pd.DataFrame({
        ID: NEW_SIX,
        "Adjudication": "incomplete reconstruction",
        "Status": "temporary reversible quarantine",
        "Permanent_global_exclusion": False,
    }).to_csv(OUT / "02_six_newly_quarantined_cells.csv", index=False)

    summary = {
        "previous_primary_cells": 193,
        "newly_quarantined_n": 6,
        "current_primary_discovery_cells": 187,
        "newly_quarantined_cells": NEW_SIX,
        "reason": "user-adjudicated incomplete reconstruction after PC1-PC3 three-SD morphology review",
        "scope": "Morphology discovery/classification only",
        "reversible": True,
        "permanent_global_exclusion_changed": False,
    }
    (OUT / "00_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
