# tWhEEL Development Log

This file summarizes repository history. It is not the usage manual; see
`README.md` for the behavior and configuration of the current source tree.

Version dates below are taken from Git tags. `src/twheel_metadata.py` is the
single source for the runtime release version, and release tags are expected to
match its `__version__` value.

## v1.3.5 (2026-09-27)

- Corrected GO-KFI chart labels to use exact category counts instead of
  reconstructing S and M from percentages rounded to one decimal place.
- Derived bar widths from exact counts so rounding neither removes nonzero
  categories nor produces a stacked total above or below 100%.
- Retained support for exact counts in older report payloads, added explicit
  errors for missing or inconsistent chart counts, and handled empty data.
- Corrected the GO chart's row label from `Contigs` to `GO terms`.
- Added generated-JavaScript regression tests covering all eleven publication
  GO statistics, small categories, empty data, and invalid or legacy payloads.
- Left GO-KFI calculation, guide filtering, and pathway analysis unchanged.
  Existing analyses need only report regeneration, with the original reports
  preserved and the analysis and renderer versions recorded separately.

## v1.3.4

- Converted supported IUPAC ambiguity codes to `N` in aligned Multi consensus
  guide-source FASTA files before CrisFlash execution.
- Added fail-closed CrisFlash validation for embedded error messages, complete
  source-record scanning, and BED output creation.
- Added a checksum-verified publication benchmark restart plan from
  `guide_design` for the twelve in-scope GO/PW analyses.

- Changed the genome-backed target-locus report notice from a prominent red
  warning to a subdued informational notice while retaining the critical red
  warning for transcriptome-only preliminary screening.
- Corrected the expanded project name to “transcript-aware Whole-genomic
  Evaluation for Editing Locus (tWhEEL)” in the README, CLI help, and generated
  GO and pathway report headers, and placed the existing logo at the top of the
  README.
- Added a Dockerized RefSeq/Ensembl `tx2protein.tsv` helper with strict FASTA
  reconciliation, audit JSON output, gzip support, selectable RefSeq XP
  inclusion, and matched protein FASTA generation for explicit exclusions.

## v1.3.2 (2026-07-17)

- Fixed CSV field escaping in standard and pathway browser exports and expanded
  regression coverage for commas, quotes, and line breaks.
- Corrected Docker documentation to describe plain `docker build`, distinguish
  compatibility pins from current upstream releases, and record the Python 3.9
  base-image migration requirement.
- Centralized release and step metadata, updated `--version` to 1.3.2, and
  expanded `--help` and `--list-steps` to cover both workflows and all eight
  implemented steps.
- Replaced historical local configurations with five portable minimal profiles
  covering de novo GO screening and existing-protein GO/PW analysis with and
  without a genome off-target reference.
- Added prominent GO/PW report warnings that distinguish preliminary
  transcriptome-only screening from genome-backed off-target searches while
  documenting `MM0 == 0` exclusion and the remaining lack of intended-locus
  and GFF/GTF exon-context reconciliation.
- Removed the `crispr.max_mismatch` YAML option, fixed CrisFlash mismatch
  searching at `-m 2`, and made stale configurations fail with an explicit
  migration error.

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
- Genome-backed runs exclude `MM0 == 0` candidates, but tWhEEL does not yet
  verify that a retained `MM0 == 1` match belongs to the intended gene locus or
  reconcile that match with GFF/GTF exon structure.
- There is no `CONTRIBUTING.md`; contribution changes should currently be
  discussed through a GitHub issue.

## License

tWhEEL is distributed under the MIT License. See `LICENSE`.
