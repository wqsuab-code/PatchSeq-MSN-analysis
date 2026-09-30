# Macaque MSN cohort and morphology quality control

Analyses were restricted to macaque medium spiny neurons annotated as STR D1 MSN, STR D2 MSN or STR Hybrid MSN and sampled from the caudate nucleus, putamen or nucleus accumbens. Cell metadata and morphology measurements were read directly from `Macaque-PatchSeq-BG.zip`, whose SHA-256 checksum was verified before analysis. Mouse morphology results and macaque electrophysiology features were not used.

The source metadata contained 486 eligible MSN profiles. Morphology classification required complete observations across 18 nonredundant dendritic and somatic measurements. Axonal variables, axonal Sholl measurements, the circular variable `axon_exit_theta_coronal`, and the derived `3_Sholl_PC1` score were excluded. A total of 126 cells from 42 donors met the complete-case criterion. No imputation was performed and no cell was removed solely because of an extreme PCA score.
