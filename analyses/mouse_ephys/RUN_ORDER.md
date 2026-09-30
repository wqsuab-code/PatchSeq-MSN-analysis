# Frozen mouse electrophysiology pipeline order

1. Verify the source manifest and place controlled source files outside the Git
   checkout as described in `data/00_source_manifest/README.md`.
2. Apply the documented stage-1 exclusions to the 556-cell main cohort to obtain
   549 cells.
3. Apply the E-QC and complete-case criteria to obtain 494 cells with 25
   candidate electrophysiological features.
4. Apply the project-wide frozen exclusion list to define the 493-cell active
   taxonomy cohort.
5. Audit skewness and redundancy; construct the frozen 18-feature
   shift-sum-log-Z matrix and the three-component PCA representation.
6. Fit HC and GC independently, align their five labels and retain the 450-cell
   GC-HC consensus cohort.
7. Generate raw-feature statistics, heatmap, six-feature radar panels and fixed
   t-SNE displays from the frozen labels and coordinates.
8. Run nested stratified linear-SVM validation and the independent
   recording-date-grouped machine-learning robustness analysis.
9. Apply the frozen strict D1/D2 identity rule to the 450 consensus cells,
   retaining 441 strict cells and nine ambiguous cells.
10. Fit the descriptive rank-3 T-E RRR in the 441 strict cells.
11. Run `python analyses/mouse_ephys/tests/test_archive_integrity.py` from the
    repository root, then regenerate `audit/FILE_INDEX.csv` and
    `audit/SHA256SUMS.txt`.

The archived scripts preserve their original analysis logic. Some scripts still
contain workstation-era paths and are not claimed as a one-command clean-clone
rerun. The integrity test is fully portable and checks all frozen conclusions
without the controlled raw files.
