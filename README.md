# Patch-seq MSN analysis

Reproducible electrophysiological, morphological and cross-modal analyses of
mouse and macaque striatal medium spiny neurons.

## Available analysis modules

### Mouse transcriptomic mapping and MSN T-type assignment

The frozen mouse transcriptomic module is located in
[`analyses/mouse_transcriptomics`](analyses/mouse_transcriptomics). It contains
the audited global neuronal RPCA stability analysis, the 588-cell MSN subtype
mapping workflow, Pearson gene-bootstrap validation, figure source data,
publication panels, provenance records and integrity checks. IN subtype mapping
is explicitly retained as pending work and is not represented as frozen.

### Macaque morphology M1-M4 and transcriptome-to-morphology analysis

The first frozen module is located in [`analyses/macaque_morphology`](analyses/macaque_morphology).
It contains the complete derived-data chain from the 486-cell target cohort to
126 complete morphology reconstructions, the frozen M1-M4 classification in 117
HC-GC consensus cells, donor-grouped machine-learning validation, M-T RRR source
data, figures and interactive editors.

Additional macaque electrophysiology, mouse electrophysiology, mouse morphology
and cross-modal modules will be added as separately audited releases.

## Data policy

Analysis-ready derived data and figures are versioned in this repository. Large
run-level outputs are distributed as versioned release assets. Primary third-party
data remain linked to their authoritative repositories unless redistribution is
explicitly permitted.

## Release status

This repository is initially intended for private collaborator review. License,
author list, permanent accessions and archival DOI must be finalized before the
first public release.
