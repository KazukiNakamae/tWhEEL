# tWhEEL Development Log

This file summarizes repository history. It is not the usage manual; see
`README.md` for the behavior and configuration of the current source tree.

Version dates below are taken from Git tags. The CLI currently reports a stale,
hard-coded `1.0.0`, so `python src/twheel.py --version` is not the release-history
source of truth.

## Unreleased (current `main` after v1.2.1)

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

- `--version` remains hard-coded to `1.0.0` and must be synchronized separately.
- `--list-steps` omits the three implemented SGD/pathway steps.
- Dated configuration files are local run records whose `LOCAL_DATA` inputs are
  not distributed in Git.
- There is no `CONTRIBUTING.md`; contribution changes should currently be
  discussed through a GitHub issue.

## License

tWhEEL is distributed under the MIT License. See `LICENSE`.
