# Panel source map

The attached images define the submission layout. The entries below identify the underlying scientific asset and whether an exact reproducible assembly exists.

| Figure | Panel | Content | Status | Frozen source |
| --- | --- | --- | --- | --- |
| Fig. 1 | a | Reference MSN t-SNE, 16 subtypes | FROZEN | `figures/main/Fig1ab_reference_and_query_source.*`; coordinates in `data/msn/MSN_PC16_per15_eta100_Pearson_Top750_bootstrap_plot_data.csv`; colors in `data/msn/MSN_fixed_subtype_color_mapping.csv` |
| Fig. 1 | b | 588 stable MSN query cells, modal Pearson label | FROZEN | Same t-SNE source and per-cell bootstrap table; labels are independent of t-SNE |
| Fig. 1 | c | Pearson versus RPCA subtype chord | FROZEN | `figures/main/Fig1c_Pearson_RPCA_chord.*`; links in `data/msn/MSN_Pearson_Top750_RPCA_subtype_chord_links.csv` |
| Fig. 1 | d | Reference/query split expression distributions and subtype proportions | FROZEN | `figures/main/Fig1d_RefQuery_expression_and_distribution.*`; counts in `data/msn/MSN_RefQuery_paired_subtype_counts_and_percentages.csv`; maxima in `data/msn/MSN_RefQuery_normalized_stacked_violin_gene_max.csv` |
| ED Fig. 1 | a | Patch-seq electrophysiology, transcriptomics and morphology workflow | LAYOUT-FROZEN | Present only in the supplied composite. Constituent microscopy, trace, schematic and reconstruction files are not yet registered |
| ED Fig. 1 | b | Detected genes in reference neurons and 643 Patch-seq libraries | SOURCE-FOUND | Reference metadata: `YZ_all_final_fixed_RNA_clean.rds`; query counts: combined gene-body matrix plus 643-cell index. Verified medians: reference 3,010; query 9,394 |
| ED Fig. 1 | c | Mitochondrial and ribosomal read percentages in 643 libraries | SOURCE-FOUND | Combined gene-body count matrix and 643-cell index; exact panel assembly script still required |
| ED Fig. 1 | d | Pearson neural versus non-neural centroid QC | FROZEN | `data/qc/Top15_risk_x040_x050_QC_table.csv`; 630 retained, 9 low-risk, 2 high-risk |
| ED Fig. 1 | e | Five-PC neural versus non-neural centroid QC | FROZEN | `data/qc/Top15CleanMarker_MaxPC_NPC5_QC_table.csv`; 623 retained, 3 low-risk, 15 high-risk |
| ED Fig. 1 | f | Sample attrition | SOURCE-FOUND | Cohort constants in `config/frozen_parameters.yml`; exact diagram assembly remains to be scripted |
| ED Fig. 1 | g | Detected genes versus reads, Good/Risk | FROZEN | `figures/extended_data_1/ED1g_reads_vs_genes.*`; plot data in `data/qc/Query_QC_reads_genes_plot_data.csv` |
| ED Fig. 1 | h | Neural Pearson correlation versus detected genes | FROZEN | `figures/extended_data_1/ED1h_neural_corr_vs_genes.*`; plot data in `data/qc/Query_QC_genes_vs_neural_Pearson_plot_data.csv` |
| ED Fig. 1 | i | Intron versus exon counts | FROZEN | `figures/extended_data_1/ED1i_exon_vs_intron.*`; plot data in `data/qc/Query_QC_exon_vs_intron_plot_data.csv` |
| ED Fig. 1 | j | Neural and non-neural marker expression | FROZEN | `figures/extended_data_1/ED1j_marker_expression.*`; plot data in `data/qc/marker_expression_Ref_Good_Risk_plot_data.csv` |
| ED Fig. 2 | a | RPCA anchor counts across the global PC scan | FROZEN | `figures/extended_data_2/ED2abc_global_RPCA_stability_correct_nPC.*`; `data/global/stage02B_RPCA_scan_summary.csv`; actual nPC values are 10, 12, ..., 30 |
| ED Fig. 2 | b | Alluvial CellType assignments across PC settings | FROZEN | Same corrected figure; `data/global/stage02B_RPCA_all_predictions.csv`; final stability excludes nPC10 |
| ED Fig. 2 | c | Stable and fluctuating broad/root assignments | FROZEN | Same corrected figure; summary tables and `data/global/ED2c_stable_variable_counts.csv` |
| ED Fig. 2 | d | Global neuronal UMAP with stable/disputed query cells | FROZEN | `figures/extended_data_2/ED2d_global_neuron_UMAP.*`; plot data in `data/global/GlobalCellType_B_Query_RPCA_Root_consensus_disputed_plot_data.csv` |
| ED Fig. 2 | e | 500-bootstrap Pearson subtype support heatmap | FROZEN-DATA | Vote fractions and per-cell tiers in `data/msn/`; supplied composite adds the Good/Risk side annotation |

## Important distinction

The t-SNE and UMAP coordinates are visualizations. They do not determine the frozen labels. Fig. 1b and Fig. 1d use modal Pearson subtype labels from the 500 gene-bootstrap analysis. Fig. 1c compares those labels with the independently frozen RPCA transfer labels.
