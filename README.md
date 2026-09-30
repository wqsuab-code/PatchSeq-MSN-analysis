Patch-seq MSN analysis

Reproducible transcriptomic, electrophysiological, morphological, and cross-modal analyses of mouse and macaque striatal neurons, with a primary focus on medium spiny neurons (MSNs).

Overview

This repository contains curated analysis code, derived data, figure source data, provenance records, and integrity checks supporting multimodal Patch-seq analyses of striatal neuron identity and phenotypic diversity. The project integrates transcriptomic, electrophysiological, and morphological measurements to examine whether D1- and D2-associated MSNs form discrete morpho-electric classes or instead occupy overlapping continua of phenotypes.

The repository is organized into species- and modality-specific analysis modules. Each module is designed to be reproducible from documented input tables or derived intermediate files, with frozen outputs, validation records, and figure source data retained for manuscript-associated review.

Large run-level outputs and primary raw sequencing files are not stored directly in this repository. These files are distributed through public data archives, release assets, or institutional storage as appropriate.

Analysis modules

Mouse transcriptomics

Location: analyses/mouse_transcriptomics

This module contains mouse transcriptomic mapping and MSN T-type assignment workflows. It includes audited neuronal RPCA stability analysis, MSN subtype mapping, Pearson gene-bootstrap validation, transcriptomic annotation tables, figure source data, publication panels, provenance records, and integrity checks.

Interneuron annotations may be retained in metadata where present, but interneuron subtype mapping is not treated as a frozen primary result unless explicitly stated in the module documentation.

Mouse electrophysiology

Location: analyses/mouse_ephys

This module contains mouse electrophysiological feature extraction, quality control, clustering, and transcriptome-to-electrophysiology analyses. It includes the derived chain from the stage-1 electrophysiology cohort to active-cell, consensus-cluster, and strict D1/D2 analysis sets; the electrophysiological taxonomy; nested and recording-date-grouped machine-learning validation; T-E RRR source data; publication figures; interactive editors; and integrity tests.

Mouse morphology

Location: analyses/mouse_morphology

This module contains mouse neuronal morphology processing, reconstruction-derived features, morphology quality control, clustering, and cross-modal analyses. It is intended to document the derived morphology cohorts, morphology-type assignments, feature matrices, validation analyses, figure source data, and provenance records used for mouse morphology and transcriptome-to-morphology comparisons.

Macaque electrophysiology

Location: analyses/macaque_ephys

This module contains macaque electrophysiological feature extraction, quality control, clustering, validation, and cross-species comparison analyses. It is intended to support macaque electrophysiology type assignments, morpho-electric comparisons where applicable, machine-learning validation, figure source data, and provenance records.

Macaque morphology

Location: analyses/macaque_morphology

This module contains macaque morphology and transcriptome-to-morphology analyses. It includes the derived-data chain from the target morphology cohort to complete reconstructions, morphology-type classification, donor-grouped machine-learning validation, morphology-transcriptome RRR source data, figures, interactive editors, and integrity checks.

Data policy

Analysis-ready derived data, figure source data, provenance records, and integrity checks are versioned in this repository. Large intermediate files, run-level outputs, raw sequencing files, and other storage-intensive artifacts are distributed separately through appropriate data archives or release assets.

Primary third-party datasets are linked to their authoritative repositories unless redistribution is explicitly permitted. Raw and processed sequencing files associated with this study should be accessed through the relevant public archive once permanent accessions are available.

Reproducibility

Each analysis module is expected to include:

documented input files or input manifests

scripts or notebooks required to regenerate derived outputs

frozen result tables used for figures and manuscript claims

provenance records describing software versions and key parameters

integrity checks for important derived files

Module-specific README files may provide additional instructions for rerunning analyses, validating outputs, or reproducing publication figures.

Release status

This repository is prepared for collaborator, reviewer, and manuscript-associated evaluation. License terms, author list, permanent data accessions, and archival DOI should be finalized before the first public release.
