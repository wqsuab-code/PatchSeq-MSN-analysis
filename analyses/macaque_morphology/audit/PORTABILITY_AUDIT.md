# Portability audit

The files below preserve exact historical code but require path/module cleanup before claiming a clean-machine one-command reproduction.

| File | Finding |
|---|---|
| `code/01_cohort_QC/audit_dedup19_missingness.py` | hard-coded Windows path |
| `code/01_cohort_QC/qc_macaque.py` | hard-coded Windows path |
| `code/01_cohort_QC/qc_macaque_morphology_msn_capunac.py` | hard-coded Windows path |
| `code/01_cohort_QC/run_m18_pca_outliers_PROVENANCE_ONLY.R` | hard-coded Windows path |
| `code/02_preprocessing_PCA/run_m18_adaptive_transform_pca.R` | hard-coded Windows path |
| `code/03_M_classification/freeze_npcs5_hck4_res23.R` | hard-coded Windows path |
| `code/03_M_classification/reproduce_macaque_M4_frozen_gc.R` | hard-coded Windows path |
| `code/03_M_classification/scan_HC_GC_complete_cases.R` | hard-coded Windows path |
| `code/04_statistics_figures/compare_all18_M_classes.py` | hard-coded Windows path |
| `code/04_statistics_figures/make_panels_A_to_I.R` | hard-coded Windows path |
| `code/04_statistics_figures/optimize_frozen_tsne.R` | hard-coded Windows path |
| `code/04_statistics_figures/plot_consensus_gc_tclass_tsne.R` | hard-coded Windows path |
| `code/04_statistics_figures/plot_consensusM4_Tclass_stacked.R` | hard-coded Windows path |
| `code/04_statistics_figures/plot_tsne_display_contraction_85pct.R` | hard-coded Windows path |
| `code/05_ML_validation/validate_macaque_M4_full_ml500.py` | hard-coded Windows path |
| `code/06_M_T_RRR/rrr_computational_helpers.py` | hard-coded Windows path; temporary extraction path |
| `code/06_M_T_RRR/rrr_plotting_helpers.py` | RRR helper imported by historical module name |
| `code/06_M_T_RRR/run_frozen_M_T_RRR.py` | RRR helper imported by historical module name |
| `code/07_interactive/build_macaque_m_rrr_editor.py` | hard-coded Windows path |
| `code/07_interactive/build_macaque_m_tsne_comparison_editor.py` | hard-coded Windows path |
| `code/08_packaging/package_frozen_main_bundle.py` | hard-coded Windows path |
| `code/08_packaging/package_ML_writing_bundle.py` | hard-coded Windows path |
| `code/08_packaging/package_tsne_comparison_statistics.py` | hard-coded Windows path |

Recommended release correction: replace absolute paths with repository-relative paths or `config/analysis_parameters.json`, extract the H5AD from the bundled ZIP into a configurable work directory, and move shared RRR functions into neutral modules that do not imply use of E-analysis data.
