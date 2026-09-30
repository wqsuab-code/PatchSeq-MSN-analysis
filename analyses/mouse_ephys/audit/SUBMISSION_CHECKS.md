# Submission checks

Checks performed on 30 September 2026 before the first module commit.

- Integrity test: PASS.
- Python syntax compilation: PASS.
- Module file count: 204 files, including manifest and checksum list.
- Module size: approximately 15.4 MiB; exact byte counts are reported from the
  final Git tree in the submission report.
- Code files: 21 R, Python or JavaScript files under `code/`.
- Data and result files: 84 files under `data/` and `results/`.
- Primary figure files: 13 under `figures/`; two representative-trace files are
  retained under `legacy/display_only/`.
- Interactive files: 43.
- Files larger than 25 MB in Git module: none.
- Files larger than 50 MB in Git module: none.
- External release assets: RNA RDS, 143,345,471 bytes; frozen review ZIP,
  8,932,345 bytes. Both match `data/00_source_manifest/release_assets.csv`.
- Absolute-path scan: five provenance or historical files identified; see
  `PORTABILITY_AUDIT.md`.
- Credential and email scan: no credential, token, password, browser cache or
  personal email detected.
- Module-mixing scan: no macaque or morphology predictor is present in the
  frozen 18-feature matrix. A historical script name records a visual/method
  reference only; see `MODULE_SEPARATION_AUDIT.md`.
- Git scope: only `analyses/mouse_ephys/` and the repository `README.md` changed.
- `analyses/macaque_morphology/`: unchanged.

The 999-permutation nested-SVM test and 200-split grouped-ML analysis were not
recomputed during this lightweight submission check. Their saved raw tables,
parameters, row counts, summaries, hashes and the 24-check recorded rerun were
validated instead.
