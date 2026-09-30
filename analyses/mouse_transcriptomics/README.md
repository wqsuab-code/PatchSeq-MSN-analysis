# Frozen transcriptomic (T) analysis framework

Freeze date: 2026-09-30

This directory is the control layer for the mouse Patch-seq transcriptomic analysis. It links the original inputs, formal classification tables, panel-level plot data, plotting code and the three supplied submission layouts. It does not redefine labels or rerun mapping.

## Authoritative analysis flow

1. **Sequenced input**: 643 Patch-seq libraries.
2. **Mapping input**: two extremely low-read libraries were removed, leaving 641 cells.
3. **Transcriptomic QC**: two independent reference-correlation diagnostics were retained. Their intersection defines 622 Good and 19 Risk cells. Risk cells remain in the analysis and are annotated.
4. **Global neuronal RPCA**: the full scan used 10, 12, 14, ..., 30 PCs. PC10 was an exploratory boundary condition. Formal stability used the ten retained settings from 12 to 30 PCs.
5. **Stable major identities**: 272 D1, 279 D2 and 34 IN cells were invariant at the broad CellType level. At the root level, 588 cells were stable MSN and 34 were stable IN; 19 were MSN-IN disputed.
6. **MSN subtype RPCA**: 20,029 reference MSNs and 588 RPCA-stable Patch-seq MSNs were analyzed using the frozen Top750/721-feature, PC16 mapping.
7. **MSN Pearson support**: subtype centroids were calculated from independently log-normalized reference data. Five hundred bootstrap replicates resampled all 721 effective genes with replacement. The modal Pearson label is the reporting label in Fig. 1b and Fig. 1d; RPCA remains the independent anchor-transfer comparison.
8. **Submission figures**: the supplied composite images are layout targets. Individual frozen panels and their data remain the authoritative scientific assets.

## Frozen headline results

- Reference neurons: 21,359 (D1 = 12,342; D2 = 7,687; IN = 1,330).
- Patch-seq libraries: 643 sequenced; 641 mapped; 622 Good; 19 Risk.
- Stable MSN query: 588; stable IN query: 34.
- MSN reference: 20,029 cells across 16 subtypes.
- Formal MSN RPCA: Top750 reference-only dispersion ranking; 721 effective genes; PC16; `k.anchor=10`; `k.filter=200`; `k.weight=30`; 2,331 anchors.
- Pearson-RPCA agreement: 422/588 exact subtype; 572/588 after D1/D2 collapse.
- Pearson bootstrap tiers: High 339, Moderate 175, Low 74.
- QC composition within those tiers: High 336 Good + 3 Risk; Moderate 167 Good + 8 Risk; Low 70 Good + 4 Risk.

## Directory map

- `reference_layouts/`: the three supplied final composite layouts.
- `figures/`: panel-level frozen source figures copied from existing analysis outputs.
- `data/`: compact tables required to audit counts, labels and panel values.
- `code/`: canonical plotting scripts used by the copied panel assets.
- `config/frozen_parameters.yml`: machine-readable parameters and expected counts.
- `docs/PANEL_SOURCE_MAP.md`: panel-by-panel source and status table.
- `docs/FREEZE_DECISIONS.md`: definitions that must not drift between figures and text.
- `docs/KNOWN_GAPS.md`: remaining items that require reconstruction or manuscript correction.
- `docs/NEXT_STEPS.md`: controlled sequence for completing the submission figures.
- `docs/SHARED_CHAT_PROVENANCE.md`: audited context recovered from the shared `Ref IN MSN Mapping` conversation, with historical and frozen workflows kept explicitly separate.
- `manifest/SHA256SUMS.txt`: checksum manifest generated after packaging.

## Status language

- **FROZEN**: data, code and output are identified and the numerical result is internally checked.
- **LAYOUT-FROZEN**: the supplied composite image is accepted as the visual target, but its assembly script is not yet present.
- **SOURCE-FOUND**: underlying data and/or component figure are located, but the final panel needs a clean assembly script.
- **PENDING**: an input, exact method field or reconstruction step is still missing.

The current MSN subtype result is frozen. IN subtype mapping is not yet frozen and must not inherit MSN-specific HVG or PC settings.

## First rebuilt panel

`figures/extended_data_2/ED2abc_global_RPCA_stability_correct_nPC.png` and its PDF are regenerated from the authoritative global tables. They use the actual nPC axis (`10,12,...,30`) and replace the ambiguous consecutive-PC labels in the supplied draft layout.
