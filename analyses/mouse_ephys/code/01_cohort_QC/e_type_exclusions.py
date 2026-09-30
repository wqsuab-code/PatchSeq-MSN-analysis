"""Shared permanent cell exclusions for all E-type analyses."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


def exclusion_file(root: Path) -> Path:
    return Path(root) / "config" / "e_type_global_exclusions.txt"


def load_excluded_cell_ids(root: Path) -> set[str]:
    path = exclusion_file(root)
    if not path.exists():
        raise FileNotFoundError(f"Global E-type exclusion file is missing: {path}")
    ids = {
        line.strip()
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    if not ids:
        raise RuntimeError(f"Global E-type exclusion file is empty: {path}")
    return ids


def filter_excluded_cells(
    frame: pd.DataFrame,
    root: Path,
    id_column: str = "MSN_unique_ID",
    *,
    require_unique_ids: bool = True,
) -> tuple[pd.DataFrame, list[str]]:
    if id_column not in frame.columns:
        raise KeyError(f"ID column {id_column!r} is absent from the input table.")
    ids = frame[id_column].astype(str)
    if require_unique_ids and ids.duplicated().any():
        duplicated = sorted(ids.loc[ids.duplicated(keep=False)].unique())
        raise ValueError(f"Duplicated cell IDs: {duplicated[:10]}")
    excluded = load_excluded_cell_ids(root)
    removed = sorted(set(ids).intersection(excluded))
    keep = ~ids.isin(excluded)
    return frame.loc[keep].copy(), removed


def assert_cells_not_excluded(cell_ids: Iterable[str], root: Path) -> None:
    excluded = load_excluded_cell_ids(root)
    overlap = sorted(set(map(str, cell_ids)).intersection(excluded))
    if overlap:
        raise ValueError(
            "Globally excluded cells cannot enter an E-type analysis: "
            + ", ".join(overlap)
        )
