# Raw and prepared source registry

Large source objects are not duplicated inside this review bundle. Their exact local paths, sizes and SHA-256 hashes are frozen here so that the copied compact tables can be traced back to the authoritative inputs.

| Role | Local source | Bytes | SHA-256 |
| --- | --- | ---: | --- |
| Reference Seurat object containing the 21,359 D1/D2/IN neurons | `C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input/05_Codex_branch_specific_mapping_input/00_original_inputs/YZ_all_final_fixed_RNA_clean.rds` | 450,509,398 | `A0FD5DFA183FE5EA966F6A2E7100392F3D995E9981AD98E61AE1F6B205C678A4` |
| Query Seurat object after two low-read removals (n = 641) | `C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input/05_Codex_branch_specific_mapping_input/00_original_inputs/query_geneBody_symbol_raw_QC_no39cells_noO9273_Batch2_noO9268_Batch2_n641.rds` | 25,919,107 | `24509C44AE2EAAA7B17014A38DF073BB4D2200CAFBE4925D2A6A3B9F63B91314` |
| Exon-only gene-by-cell count matrix | `C:/Users/53461/Downloads/Patchseq_Step4_for_local/Patchseq_Step4_for_local/Step4_Symbol_conversion/matrices/exon_only_counts_matrix_symbol_sum.tsv` | 82,739,382 | `D5797B7B8AA00BC33B4F9801E70DDA8F5D234B9FDA418A6F865DB6AFFAE3EEA5` |
| Combined gene-body exon-plus-intron count matrix | `C:/Users/53461/Downloads/Patchseq_Step4_for_local/Patchseq_Step4_for_local/Step4_Symbol_conversion/matrices/gene_body_exon_intron_counts_matrix_symbol_sum.tsv` | 83,489,342 | `C6E354CBD923DEDA9389A9AD02A4E15793F1EF2BF9562771FAB3BD96E0865BC5` |
| Query PCA object entering global RPCA | `C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input/05_Codex_branch_specific_mapping_input/01_global_neuron_objects/stage02A_query_PCA_object.rds` | 131,886,010 | `4CD70F9DE4A16239140071061143631FEC282D172828E395855104B9366BF355` |
| Prepared MSN Top750/721-feature reference-query pair | `C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs/05_branch_specific_subtype_mapping/60_MSN_PC16_tSNE_complete_package/inputs_prepared/MSN_top750_dispersion_z_pca_pair.rds` | 60,587,460 | `D0F0ADCC3DE6E58A0A45D2FFD664C9A411A3F84274A66D72F5B0F3182D6C3AF3` |

## External metadata fields still pending

The following items cannot be inferred safely from these expression objects and must be confirmed from the primary sequencing workflow records: adapter-trimming software and parameters, genome assembly, aligner version and parameters, gene annotation release, and the exact read-counting software/version.

## Historical analysis provenance

The shared conversation `https://chatgpt.com/share/6abc9d75-4fb0-83ea-b358-c04be90ffb24` documents earlier server-side source paths and a broad-class Seurat/SingleR bootstrap workflow. See `SHARED_CHAT_PROVENANCE.md`. Those paths and settings are historical provenance, not replacements for the checksum-verified inputs above.
