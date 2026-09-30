# Run order

The archive separates frozen consistency checks from computationally expensive
reruns. Commands below are run from `analyses/mouse_morphology/`.

1. **Fast archive validation**

   ```text
   python tests/test_archive_integrity.py
   ```

   This validates manifests, cohort transitions, feature order, IDs, class
   counts, PCA products, HC-GC consensus, ML/RRR cohorts, missingness and module
   boundaries without recomputing the 500-repeat suites.

2. **Cohort and reconstruction QC**

   Review scripts in `code/01_cohort_QC/` in numerical workflow order:
   morphology availability/QC, year-aware ASC mapping, E/T audit, 193-cell
   outlier screen and the six-cell incomplete-reconstruction quarantine.

3. **Feature preprocessing and redundancy audit**

   Use `code/02_preprocessing/`. The authoritative transform is recorded in
   `config/frozen_analysis.json`. Historical robust-scaling/candidate feature
   configurations under `legacy/` are not primary inputs.

4. **Classification and consensus**

   Use `code/03_classification/` to reproduce the NPC/HC/GC scan, final
   HC-derived taxonomy, label alignment and bootstrap stability. Frozen inputs
   and outputs are under `data/02_frozen_input/`, `data/03_transformed/`,
   `data/04_frozen_classification/` and `results/clustering_scan/`.

5. **Statistics and figures**

   Use `code/04_statistics_figures/`. Figure source tables and frozen exports
   are retained under `results/final_taxonomy/` and `figures/`.

6. **Machine-learning validation**

   Use `code/05_ML_validation/`. All transforms and PCA fits must be estimated
   inside each training split. Recording-day groups must never cross train/test
   partitions. The saved 500-repeat results are validated by row counts, split
   logs, summaries and hashes; they are not recomputed by the fast test.

7. **Strict T-M RRR**

   Use `code/06_cross_modal/`. Only the strict n=168 RPCA-Pearson-concordant
   D1/D2 plus HC-GC-consensus cohort is valid. Do not substitute the broad
   descriptive n=180 D1/D2 table.

8. **Interactive editors and packaging**

   Editors are self-contained under `interactive/`. Packaging and metadata
   normalization utilities are under `code/08_packaging/`.

## Portability note

Some archived provenance scripts and source tables retain original Windows
paths. They document the historical execution but are not invoked by the clean-
clone integrity test. Remaining paths are enumerated in
`audit/absolute_path_scan.txt`; primary frozen data and tests are repository-
relative.
