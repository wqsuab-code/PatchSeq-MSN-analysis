# Archive scope, conflict disposition and open gaps

## Freeze decision

The 2026-09-23 release is authoritative and supersedes the 2026-09-19 source
package. Frozen assignments, summary, strict RRR tables and ML web data are
byte-identical to their final workspace counterparts. The final composite is
`figures/final/Mouse_M_20260920_Ver3.svg`.

## Conflicts found and resolved without changing conclusions

1. **Broad versus strict D1/D2 scope.** The broad descriptive table contains 180
   HC-GC-consensus D1/D2-labelled cells. Strict T-M RRR contains 168 cells after
   requiring RPCA-Pearson agreement. The tables remain separate and tests prevent
   substituting n=180 for n=168.
2. **Historical feature configurations.** `morph_primary_features_v1.json` and
   candidate core-feature files describe earlier robust-scaling or candidate
   analyses. They are isolated under `legacy/config/`; the primary ten-feature
   order and shift-sum-10,000/log1p/Z-score transform are frozen in
   `config/frozen_analysis.json`.
3. **Figure version.** Earlier `Mouse_M_Ver09182026.tif` is superseded by the
   20 September Ver3 SVG. The earlier image is not presented as final.

None of these dispositions changes the frozen 187-cell taxonomy or the strict
168-cell RRR conclusions.

## Open provenance and portability gaps

- Original ASC files were on an unmounted external drive and are not included.
- A permanent mouse data accession has not yet been assigned.
- The standalone training generator for the strict RPCA-Pearson identity
  classifier has not been recovered. Frozen predictions and the downstream
  168-cell RRR selection are retained and auditable.
- Exact R/Seurat package versions and an `renv.lock` were absent from the frozen
  materials. Python versions are recorded.
- Some provenance scripts/tables retain original absolute Windows paths; these
  are enumerated in `absolute_path_scan.txt` and are not required by the fast
  clean-clone integrity test.
- True mouse IDs are unavailable. Recording day (105 groups) is an explicit
  proxy because one mouse was used per day; no animal identities were invented.
