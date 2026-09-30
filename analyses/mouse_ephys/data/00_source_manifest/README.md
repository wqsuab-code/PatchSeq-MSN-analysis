# Controlled source files

Original electrophysiology workbooks, ABF recordings and the transcriptomic
count object are not committed to ordinary Git history. `source_files.csv`
records their immutable filenames, byte sizes and SHA-256 hashes. Permanent
repository accessions and redistribution terms remain to be supplied before a
public release.

The analysis-ready `stage1_filtered.json` is committed because it is the exact
549-cell derived input used by the archived preprocessing script. The 143 MB RNA
RDS and the frozen review ZIP are staged outside Git in the local
`PatchSeq-MSN-analysis-mouse-ephys-release-assets` directory.

The ten ABF files listed in the manifest support the displayed representative
traces only. The trace set is retained as display provenance because its original
candidate screen included an expert-label criterion; it did not define the
GC-HC taxonomy.
