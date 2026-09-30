from __future__ import annotations

import base64
import csv
import io
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MORPH = ROOT / "interactive/mouse-morph-feature-explorer/dist/data.csv"
AUDIT = ROOT / "outputs/20260820_REFonly_D1D2_marker_classifier_scan/D1D2_marker_classifier_Query_predictions.csv"
OUT = ROOT / "interactive/mouse-morph-feature-explorer/dist/strict-identity-data.js"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


morph = read_rows(MORPH)
audit = {row["cell_id"]: row for row in read_rows(AUDIT)}

columns = [
    "MSN_unique_ID",
    "strict_D1D2_identity",
    "strict_D1D2_pass",
    "audit_present",
    "agree_RPCA",
    "agree_Pearson",
    "marker_major",
    "pearson_major",
    "bootstrap_top1_fraction",
    "bootstrap_margin",
]
buffer = io.StringIO(newline="")
writer = csv.DictWriter(buffer, fieldnames=columns, lineterminator="\n")
writer.writeheader()

for cell in morph:
    cell_id = cell["MSN_unique_ID"]
    row = audit.get(cell_id)
    stable = row is not None and row.get("Final_CellType_stability") in {"Stable_D1", "Stable_D2"}
    rpca = row is not None and row.get("agree_RPCA", "").upper() == "TRUE"
    pearson = row is not None and row.get("agree_Pearson", "").upper() == "TRUE"
    strict = stable and rpca and pearson
    writer.writerow(
        {
            "MSN_unique_ID": cell_id,
            "strict_D1D2_identity": row.get("Final_consensus_CellType", "") if strict else "",
            "strict_D1D2_pass": str(strict),
            "audit_present": str(row is not None),
            "agree_RPCA": row.get("agree_RPCA", "") if row else "",
            "agree_Pearson": row.get("agree_Pearson", "") if row else "",
            "marker_major": row.get("Marker_major", "") if row else "",
            "pearson_major": row.get("Pearson_major", "") if row else "",
            "bootstrap_top1_fraction": row.get("bootstrap_top1_fraction", "") if row else "",
            "bootstrap_margin": row.get("bootstrap_margin", "") if row else "",
        }
    )

payload = base64.b64encode(buffer.getvalue().encode("utf-8")).decode("ascii")
OUT.write_text(f"window.MOUSE_M_STRICT_IDENTITY_DATA=atob('{payload}');\n", encoding="utf-8")
print(f"Wrote {OUT} with {len(morph)} morphology cells; {sum(1 for r in morph if audit.get(r['MSN_unique_ID']))} audit matches.")
