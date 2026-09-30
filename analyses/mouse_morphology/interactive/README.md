# Interactive editors

- `feature-explorer/index.html`: frozen ten-feature comparison editor.
- `feature-explorer/rrr-editor.html`: strict 168-cell T-M RRR editor.
- `feature-explorer/d1d2-distribution-editor.html`: D1/D2 morphology distribution editor.
- `feature-explorer/distribution-statistics.html`: cohort and class counts.
- `ml-studio/index.html`: ML validation figure editor.

Serve this directory with a local static web server; opening pages directly by
`file://` may prevent the browser from loading companion CSV/JavaScript data.
The editors are visualization layers over frozen data and do not refit the
taxonomy, ML models or RRR.
