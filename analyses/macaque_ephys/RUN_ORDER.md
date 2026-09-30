# Audited run order

1. Verify the public data release using
   `data/00_source_manifest/SOURCE_ARCHIVE_MANIFEST.json`. Raw third-party data
   are not included in Git.
2. Build and audit the Ca/Pu/NAc MSN stage-1 cohort with
   `code/01_cohort_QC/run_dstr_dstrvstr_rheobase_qc.py` and the accompanying QC
   scripts.
3. Audit redundancy, missingness and skewness; verify the ordered final19
   feature list in `config/analysis_parameters.json`.
4. Run `code/02_preprocessing/reproduce_frozen_preprocessing_pca_hc.py` to
   reconstruct the 390x19 transformed matrix, full PCA and Ward K=4 partition.
5. Verify the archived GC15 labels, surjective merge and HC-GC consensus with
   `code/03_classification/validate_frozen_consensus.py`.
6. Generate E-class statistics and static figures from the frozen 368-cell
   consensus cohort using `code/04_statistics_figures`.
7. Run or audit the 500-repeat donor-grouped validation in
   `code/05_ML_validation`. The primary transform has
   `yeo_johnson=false`; older Yeo-Johnson runs are sensitivity analyses.
8. Refit or verify the 346-cell D1/D2 T-E RRR using
   `code/06_cross_modal/plot_ca_pu_nac_rrr_fourpanel.py`.
9. Build interactive pages from archived result tables using
   `code/07_interactive`.
10. Run `tests/test_archive_integrity.py`, regenerate `audit/FILE_INDEX.csv`
    and `audit/SHA256SUMS.txt`, and perform the portability and sensitive-data
    scans before packaging.

The final19 Seurat GC step cannot currently be regenerated from a single
authoritative historical script. Do not substitute a generic scan script or an
older NAc-only analysis. The frozen assignments and seed scans are the
authoritative record until the exact historical R source is recovered or a
new, independently validated reproduction script is approved.
