# Archive review outcome

## Frozen authority

The 23 September 2026 frozen bundle is the authoritative output snapshot. The
24 September 2026 recorded R/Python rerun reproduced all 24 prespecified results
and therefore validates, rather than supersedes, that snapshot. The 29 September
review package supplied the public directory structure and source-data subset.

## Conflict decision

No conflict was found that would change the five frozen labels, cohort sizes,
statistics, machine-learning conclusions or strict-441 RRR. Historical and
sensitivity alternatives are separated in `legacy/` and documented in
`CONFLICT_RESOLUTION.md`.

## Known limitations

- Permanent raw-data accession identifiers and access dates are pending.
- A unique animal identifier is not available in the archived E-analysis tables;
  recording date is the grouped-validation unit and sequencing batch is metadata.
- Two current plotting scripts still require path refactoring for clean-machine
  reruns.
- The representative traces are descriptive provenance and were not selected by
  a fully label-blind GC-HC-only procedure.
- License, authorship metadata and repository DOI remain pending.
