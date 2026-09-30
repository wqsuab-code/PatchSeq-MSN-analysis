# Conflict and concordance audit

## Overall decision

The latest frozen Mouse M analysis is internally concordant. The 187-cell raw table, final assignment table and interactive feature table contain the same identifiers and identical values for all ten frozen features. Recomputed PCA and Ward.D2 K=4 clustering reproduce the stored taxonomy exactly. GC at resolution 0.50 agrees with HC for 181 of 187 cells (96.79%). The strict T-M RRR score table contains exactly the 168 cells defined by HC-GC consensus plus stable RPCA and Pearson D1/D2 agreement.

## Verified agreements

- Cohort: 187 rows, 187 unique IDs; class counts {'M1': 67, 'M2': 23, 'M3': 62, 'M4': 35}.
- Raw versus final frozen features: maximum absolute difference 0.
- Final versus interactive frozen features: maximum absolute difference 0.
- PCA reproduction: maximum absolute difference after allowing component sign flips 4e-14.
- HC reproduction: 187/187 labels, ARI=1.000.
- HC-GC agreement: 181/187 (96.79%).
- Strict T-M RRR: n=168, D1=74, D2=94; score IDs match exactly.

## Resolved semantic conflict

The file `11_D1D2_counts_by_final_M_class.csv` and the `D1_D2` column in the 187-cell table use a broader descriptive identity and yield 180 HC-GC-consensus D1/D2 cells. Strict T-M RRR requires stable D1/D2 plus agreement of both RPCA and Pearson classifiers and yields 168 cells. The difference is 12 cells. These files are not numerically interchangeable. The package labels the broad table as descriptive and the strict n=168 table as authoritative for RRR.

## Version conflict resolved

The 19 September source package contained the earlier `Mouse_M_Ver09182026.tif`. This release supersedes it with `Mouse_M_20260920_Ver3.svg` and the final 20 September feature-comparison layout. The three downloaded copies of the final M-feature layout are byte-identical; only one canonical copy is retained.

## Raw reconstruction limitation

The external ASC root was not available during packaging, so zero ASC files were copied. The frozen mapping audit contains 183 selected filenames for the 187-cell cohort, and the archived NeuroM audit documents 183 year-aware matches and 182 successful geometry recomputations. `RAW_ASC_acquisition_status_for_final187.csv` is the exact checklist for adding the files later without changing IDs or analytical outputs.
