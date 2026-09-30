# Frozen pipeline order

1. Verify and extract `data/00_source/Macaque-PatchSeq-BG.zip`.
2. Run the cohort and feature-QC scripts in `code/01_cohort_QC`.
3. Run `code/02_preprocessing_PCA/run_m18_adaptive_transform_pca.R` in all-126 mode.
4. Run `code/03_M_classification/scan_HC_GC_complete_cases.R` with the 126-cell cohort.
5. Freeze NPC5/HC-K4/GC-resolution-2.3 with `freeze_npcs5_hck4_res23.R`.
6. Reproduce HC and GC labels and build the 486-cell audit.
7. Generate class statistics and figures from the frozen 117-cell consensus cohort.
8. Run the 500-repeat donor-grouped validation.
9. Run M-T RRR from the source H5AD, frozen 117 IDs and frozen 18-feature Z matrix.
10. Build interactive pages and regenerate manifests/checksums.

The exact historical scripts are retained for audit. A fully portable one-command
runner should only be declared after the absolute paths listed in
`audit/PORTABILITY_AUDIT.md` have been replaced by configuration-driven paths.
