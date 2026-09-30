# Frozen-version conflict review

No unresolved conflict was found that changes the frozen cell labels, feature
matrix, classification counts, statistics or T-E RRR conclusions. The
23 September frozen bundle and the 24 September recorded rerun agree on all 24
prespecified checks.

| Topic | Candidate versions | Resolution |
| --- | --- | --- |
| Stage-1 cohort | 549 and historical 551 | Use 549; the JSON and final Methods agree. Retain only the 551 filename and checksum as historical provenance. |
| D1/D2 composition | 450 direct major-class assignments versus strict 441 plus 9 ambiguous | Use the strict stability rule for current figures and cross-modal analysis. |
| Transformation | shift-sum-log-Z versus Yeo-Johnson | Use shift-sum-log-Z for the primary taxonomy; retain Yeo-Johnson as sensitivity only. |
| ML performance | 92.0% nested cell-stratified accuracy versus 91.33% aggregated recording-date-grouped accuracy | Report each with its own validation design and estimand. |
| Radar features | six-feature main panel versus optional nine-feature editor | Use six features in the main figure; mark nine-feature controls as optional display functionality. |
| RRR class version | merged GC K = 5/HC K = 5 versus GC-unmerged K = 11/HC K = 13 | Use merged K = 5 as primary; keep the unmerged model as exploratory. |
| Representative traces | original candidate screen included an expert-label criterion | Retain only as descriptive display provenance; do not use as evidence defining GC-HC classes. |

Legacy `HC_GC` filenames are retained for traceability, while narrative text uses
GC-HC to match the stated analysis order.
