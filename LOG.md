# tWhEEL Development Log

This file summarizes repository history. It is not the usage manual; see
`README.md` for the behavior and configuration of the current source tree.

Version dates below are taken from Git tags. `src/twheel_metadata.py` is the
single source for the runtime release version, and release tags are expected to
match its `__version__` value.

## Unreleased (current `main` after v1.3.1)

- Fixed CSV field escaping in standard and pathway browser exports and expanded
  regression coverage for commas, quotes, and line breaks.
- Corrected Docker documentation to describe plain `docker build`, distinguish
  compatibility pins from current upstream releases, and record the Python 3.9
  base-image migration requirement.
- Centralized release and step metadata, updated `--version` to 1.3.1, and
  expanded `--help` and `--list-steps` to cover both workflows and all eight
  implemented steps.
- Replaced historical local configurations with five portable minimal profiles
  covering de novo GO screening and existing-protein GO/PW analysis with and
  without a genome off-target reference.
- Added prominent GO/PW report warnings that distinguish preliminary
  transcriptome-only screening from genome-backed off-target searches while
  documenting the remaining transcript-derived splice-junction limitation.

## v1.3.1 (2026-07-17)

- Added immutable, manifest-led normalized gene/guide generations with atomic
  current-generation selection and strict source/artifact validation.
- Added content-tagged multi-gene alignment images, bounded MAFFT workers,
  validated per-group consensus caching, and shared workspace locks.
- Updated GO helper contracts to consume normalized gene rows and require an
  explicit normalized-schema capability label where applicable.
- Added versioned browser-payload generations, deferred guide chunks, atomic
  HTML publication, and script-safe serialization for standard and pathway
  reports.
- Strengthened pathway identifier reconciliation, Docker workspace path
  translation, report integrity checks, and automated coverage.

## v1.3.0 (2026-07-16)

- Added immutable, checksum-verified report-source snapshots under
  `web_data/report_sources/`.
- Reduced `twheel_data.db` to a minimal GO label/definition and report-metadata
  cache; native analysis tables remain in TSV, CSV, FASTA, BED, text, and GPML
  files.
- Added atomic cache publication and integrity checks.
- Refactored the pathway report to use the shared standard-report shell and
  static assets.
- Added report regeneration contracts and validation coverage for standard and
  pathway generators.

## v1.2.1 (2026-07-16)

- Consolidated protein input-mode resolution across the orchestrator and
  annotation modules.
- Added strict duplicate-key rejection for YAML configuration files.
- Strengthened output-directory reuse, optional-genome handling, and pathway
  analysis handoff validation.
- Added namespace-specific, validated FASTA-search caches.
- Tightened Multi-group alignment and site-specific guide contracts.
- Added exact pathway identifier matching, canonical gene-unit counting, audit
  outputs, and report validation tests.
- Expanded automated coverage across GO, guide, pathway, report, and Docker
  helper behavior.

## v1.2.0 (2025-11-18)

- Added the SGD-centered pathway workflow.
- Added WikiPathways GPML parsing and PW-KFI calculation.
- Added pathway result tables, mappings, and a static pathway report.
- Added `sgd_annotation`, `pathway_analysis`, and
  `pathway_report_generation` step dispatch.

## v1.1.0 (2025-09-10)

- Added support for supplied protein sequences as an alternative to
  TransDecoder output.
- Added protein-to-transcript mapping support for downstream transcript
  sequence extraction.
- Added organism-specific configuration examples and annotation changes.

## v1.0.0 (2025-08-12)

- Established the transcriptome-to-GO workflow:
  gene prediction, protein annotation, guide design, GO-KFI analysis, and a
  static report.
- Added Docker integration for the external bioinformatics tools and GO helper
  scripts.

## Current limitations relevant to release history

- Release preparation must update `src/twheel_metadata.py` and create a matching
  Git tag; this consistency is not yet enforced by release automation.
- tWhEEL does not yet reconcile transcriptome-derived guide candidates to exon
  structure or contiguous genomic target coordinates.
- There is no `CONTRIBUTING.md`; contribution changes should currently be
  discussed through a GitHub issue.

## License

tWhEEL is distributed under the MIT License. See `LICENSE`.
