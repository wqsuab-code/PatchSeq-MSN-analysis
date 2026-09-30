# Portability audit

The portable integrity test runs from a clean clone using only the Python
standard library. A complete raw-data rerun is not yet one-command portable.

## Remaining absolute paths

Absolute Windows paths remain in the following archived provenance or historical
files:

- `data/04_frozen_classification/GC_parameters_and_run_log.txt`
- `results/05_display_source/E1-E5_core6_radar_parameters_and_QA.txt`
- `results/05_display_source/run_log.txt`
- scripts under `legacy/historical/`

These paths are retained as execution provenance or in explicitly historical
scripts. Current heatmap, radar and strict-441 RRR scripts resolve module-relative
paths. `config/paths.example.json` defines the external-data locations.

## External dependencies

- The strict-441 RRR requires `MOUSE_E_RNA_RDS` to point to the external RNA RDS
  listed in the source manifest.
- Representative-trace regeneration requires the ten external ABF files.
- Historical interactive editors expect a local static web server; they do not
  alter frozen analytical results.
- The exact recorded R and Python package versions are archived under
  `environment/`.

## Verification level

The final 24-result concordance audit is preserved from the recorded full rerun.
This repository commit additionally performs a frozen consistency audit of
hashes, row counts, feature order, IDs, labels, statistics, ML and RRR outputs.
It does not repeat the 999-permutation or 200-split computations during CI.
