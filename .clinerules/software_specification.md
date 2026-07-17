# tWhEEL Implementation Specification

This document describes the behavior implemented in the current source tree.
It is an implementation reference, not a roadmap. User-facing setup and usage
belong in `README.md`.

## 1. Scope

tWhEEL is a sequential, configuration-driven transcriptome workflow with two
supported analysis profiles:

```text
GO profile
gene_prediction -> annotation -> guide_design -> go_analysis
                -> report_generation

Pathway profile
gene_prediction -> sgd_annotation -> guide_design -> pathway_analysis
                -> pathway_report_generation
```

The two profiles should use separate run directories. `annotation` and
`sgd_annotation` publish the same intermediate artifact names, so they are
alternatives rather than safe parallel branches within one run.

## 2. Runtime architecture

`src/twheel.py` provides the CLI and the `tWhEELPipeline` orchestrator. It:

1. loads YAML with duplicate-key rejection;
2. derives a timestamped or reused run directory from `output.base_dir` and
   `--name`;
3. sets `output.intermediate` to `<run>/intermediate`;
4. creates the run, intermediate, and web-data directories;
5. validates full-run inputs; and
6. invokes configured steps sequentially.

Most external bioinformatics tools run through the Docker CLI via
`subprocess`. The report modules are local Python code and generate static HTML,
JSON, assets, and a small SQLite GO metadata cache.

### 2.1 Implemented step dispatch

`tWhEELPipeline.run_step()` accepts:

- `gene_prediction`
- `annotation`
- `sgd_annotation`
- `guide_design`
- `go_analysis`
- `pathway_analysis`
- `pathway_report_generation`
- `report_generation`

An unknown step returns failure. A single `--step` invocation does not execute
prerequisites and does not call the full-run `validate_inputs()` method.

### 2.2 Output-directory contract

- A simple name creates
  `<output.base_dir>/tWhEEL_on_<timestamp>_<name>`.
- No name creates `<output.base_dir>/tWhEEL_on_<timestamp>`.
- A complete name beginning with `tWhEEL_on_` is reused below the configured
  output root.
- If `output.base_dir` already names that complete run directory, it is not
  nested a second time.
- `output.base_dir` must be a non-empty string. A leading `~` is expanded.
- The configured `output.intermediate` value is replaced at runtime.

Directory reuse is not checkpoint detection. A full run executes every entry
in `execution.steps` again.

## 3. Configuration contract

### 3.1 Configuration sections by workflow

Every complete pipeline run uses `input`, `output`, and `execution` mappings.
Additional mappings are required by the steps that read them:

| Step or profile | Required mapping |
| --- | --- |
| `gene_prediction` | `busco`; `transdecoder.docker_image` only when TransDecoder mode is effective |
| `annotation` or `sgd_annotation` | `fasta_search` |
| `guide_design` | `crispr` |
| `go_analysis` | `go_analysis` |
| `report_generation` | `report` |
| `pathway_analysis` | `pathway_analysis` |

`pathway_report` is optional; `pathway_report_generation` defaults its HTML
name to `pathway_report.html`. `database` is optional and is read only by the
standard report generator. A pathway-only profile does not read
`go_analysis`, `report`, or `database`.

The portable configuration surface consists of five standalone templates:

| File | Profile |
| --- | --- |
| `config/default.yaml` | Minimal de novo GO-KFI, with no supplied protein or genome |
| `config/go_kfi_existing_protein.yaml` | Existing-protein GO-KFI without a genome |
| `config/pw_kfi_existing_protein.yaml` | Existing-protein SGD/WikiPathways PW-KFI without a genome |
| `config/go_kfi_existing_protein_genome.yaml` | Existing-protein GO-KFI with a genome off-target reference |
| `config/pw_kfi_existing_protein_genome.yaml` | Existing-protein SGD/WikiPathways PW-KFI with a genome off-target reference |

The no-genome profiles are limited to preliminary knockout-feasibility
screening. The genome-backed profiles are the recommended starting point for
any use beyond that preliminary screen.

For a full orchestrated run, including `pathway_report_generation` in
`execution.steps` also makes
`pathway_analysis.wikipathways_gpml_dir` an input-validation requirement.

### 3.2 Input keys

| Key | Contract |
| --- | --- |
| `input.transcriptome` | Required for every full run. Used in both protein input modes. |
| `input.protein_sequences` | Optional existing protein FASTA. |
| `input.protein_transcript_mapping` | Required by annotation modules when existing protein mode is effective. |
| `input.genome_assembly` | Runtime-optional off-target reference; transcriptome is the fallback, but that mode is limited to preliminary feasibility screening. |
| `input.protein_db_for_goanalysis` | Required by full-run validation when `annotation`, `go_analysis`, or `report_generation` is configured. |
| `input.sgd_db` | Required when `sgd_annotation` is configured. Also used for exact pathway cross-reference resolution. |

The existing-protein mapping is a headerless TSV. Blank lines and lines whose
first character is `#` are skipped. The first two columns are:

```text
transcript_id<TAB>protein_id
```

Extra columns are ignored. The loader stores `protein_id -> transcript_id`.
The mapping protein IDs must exactly cover the prepared protein FASTA IDs, and
every mapped transcript ID must exist as a non-empty prepared transcript FASTA
record. Additional transcript FASTA records are allowed. Any mapping format or
identifier inconsistency stops annotation before the similarity search.

### 3.3 Effective protein input mode

`src/modules/input_mode.py` is the shared source of truth:

- existing protein mode is selected when `input.protein_sequences` is non-empty,
  the file exists, and `transdecoder.skip_if_protein_provided` is true;
- `skip_if_protein_provided` defaults to true if omitted; and
- all other cases select TransDecoder mode.

The transcriptome is required in both modes.

### 3.4 Pathway configuration

`pathway_analysis` requires:

```yaml
pathway_analysis:
  wikipathways_gpml_dir: "/path/to/gpml"
  target_gene_databases:
    - database_name: "SGD"
      id_pattern: "^S\\d{9}$"
      annotation_key: "sgd_id"
  minimum_pathway_coverage: 0.5
  skip_low_coverage_pathways: true
```

Every target database entry requires a unique, non-empty `database_name`, a
valid `id_pattern`, and an `annotation_key`. GPML database names match exactly,
and identifiers must fully match the configured regular expression. Only
`*.gpml` files directly inside the configured directory are discovered.

### 3.5 Effective and legacy settings

- `report.generate_html` and `report.generate_summary` are effective.
- `report.guide_chunk_rows` is effective and defaults to `5000`.
- `pathway_report.output_name` is effective and defaults to
  `pathway_report.html`.
- `database.force_rebuild` controls rebuilding the minimal GO cache.
- `execution.verbose` controls the logging level.
- `execution.parallel` is currently ignored; execution is sequential.
- `docker.max_memory` and `docker.use_docker` are currently ignored.
- `report.web_port`, `pathway_report.web_port`, and
  `pathway_report.generate_html` are currently ignored. No web server is
  started.

### 3.6 Docker image and dependency contract

TransDecoder, BUSCO, FASTA3, and CrisFlash image tags are configuration values.
SeqKit is not configurable: `GenePrediction.get_sequence_stats()` uses the
hard-coded `quay.io/biocontainers/seqkit:2.8.0--h9ee0642_0` image. GO helper
image tags are configuration values, but missing images are built from the
repository-local contexts with plain `docker build`.

The configured values are compatibility pins rather than an automatically
updated software set. As reviewed on 2026-07-17:

- BUSCO `v5.2.2_cv1` predates the current BUSCO 6.1.0 stable release;
- SeqKit 2.8.0 predates the current upstream 2.13.0 release;
- TransDecoder 5.7.1 remains the upstream GitHub latest release, while Bioconda
  also publishes 6.0.0; and
- the configured FASTA3 version remains in the current 36.3.8 release family,
  although newer package build revisions exist.

These observations do not authorize automatic replacement: BUSCO and
TransDecoder upgrades can change datasets, output layout, or biological
results and require representative regression analyses.

`add_go_columns`, `compute_go_kfi`, and `run_alignment_and_mask` currently use
`python:3.9-slim`; `count_go_parent_terms` uses `python:3.11-slim`. Python 3.9
reached end of life on 2025-10-31. Base images are not digest-pinned, and most
helper Python dependencies are unpinned or use lower bounds. Consequently, a
later rebuild of the same local tWhEEL helper tag is not guaranteed to reproduce
the same image contents. No minimum Docker Engine or Docker Desktop version is
currently established by automated integration testing.

## 4. Module contracts

### 4.1 `GenePrediction`

File: `src/modules/gene_prediction.py`

Responsibilities:

- normalize the transcriptome to `intermediate/transcriptome.fa`;
- use an existing protein FASTA or run TransDecoder;
- publish proteins as `intermediate/query.pep`;
- collect SeqKit statistics;
- run BUSCO in protein mode; and
- write `gene_prediction_results.txt`.

The TransDecoder minimum protein length defaults to 100 amino acids and must be
a positive integer. TransDecoder, SeqKit, and BUSCO run in Docker. BUSCO is
invoked with `--platform linux/amd64`.

### 4.2 `Annotation` and `Sgd_annotation`

Files:

- `src/modules/annotation.py`
- `src/modules/sgd_annotation.py`

Both modules:

- require `query.pep`;
- require the mapping TSV in effective existing protein mode;
- copy or decompress their configured reference FASTA;
- run `ggsearch36` with one best hit and `-m10` output;
- parse best-hit scores and full descriptions;
- classify hits by shared `subject_id`;
- extract Single and Multi transcript sequences;
- align Multi groups and create consensus sequences; and
- write the common annotation artifacts.

The standard module uses `protein_db_for_goanalysis`; the SGD module uses
`sgd_db` and preserves SGD/Ensembl identifiers for pathway matching.

Classification is exact:

- `Single`: one query row for a best-hit `subject_id`;
- `Multi`: more than one query row for that `subject_id`.

The generated `Single_*` or `Multi_*` identifier is the downstream gene-level
unit. The classification does not independently establish biological gene
redundancy.

FASTA searches use namespace-specific validated caches:

- `ggsearch_results_annotation.txt` and its `.meta.json`;
- `ggsearch_results_sgd.txt` and its `.meta.json`.

Cache identity includes the query and database real paths, sizes, modification
times, namespace, and search settings. Results and metadata are published only
after a successful non-empty search.

Multi-group alignment builds or reuses a content-tagged local Docker image,
then runs a bounded number of single-threaded MAFFT processes. Validated
consensus results are cached by ordered input content, algorithm schema, and
immutable image identity. Every expected group must produce exactly one
non-empty consensus FASTA with the expected identifier. Workspace generation
locks and staged publication prevent concurrent or failed runs from exposing a
mixed output generation.

### 4.3 `GuideDesign`

File: `src/modules/guide_design.py`

Responsibilities:

- copy or decompress the genome assembly for the CrisFlash `-g` off-target
  reference, or fall back to the transcriptome;
- convert supported IUPAC degenerate bases to `N`;
- prepare transcriptome-derived Single sequences and validated Multi consensus
  sequences for the CrisFlash `-s` candidate-target input;
- run CrisFlash separately for the available Single/Multi groups;
- retain site-specific targets; and
- publish one normalized gene/guide generation.

`normalized_genes.tsv` contains exactly one row per `unique_id`, including the
deterministic representative alignment, all query IDs, annotation row count,
guide count, and explicit database identifiers. `normalized_guides.tsv`
contains one row per distinct typed guide tuple. A disk-backed SQLite index
de-duplicates guides while the generation is built, so peak Python memory does
not scale with the full guide set.

Each immutable generation also contains a compact, one-row-per-gene
`annotation_with_guides_final.tsv` for compatibility. Contig-level provenance
remains in `annotation_classified.txt`; the legacy compatibility table is no
longer the authoritative downstream representation. SHA-256 metadata covers
the normalized artifacts, and `normalized_gene_guide_current.json` is replaced
only after a complete generation has been written. A failed publication leaves
the previous current generation selected. Source paths, sizes, modification
times, and SHA-256 hashes bind the generation to its annotation and guide BED
inputs; stale current data fails closed.

CrisFlash target field four has the form
`<sequence>/<MM0>/<MM1-2>`. A site-specific target is retained only when
`MM0 == 1` and `MM1-2 == 0`. Other CrisFlash scores are passed through without
a tWhEEL-enforced range.

A supplied genome improves only the off-target reference. Candidate sequences
remain transcriptome-derived and are not reconciled to genomic exon coordinates
or splice structure. A candidate may therefore cross an exon-exon boundary or
otherwise lack a contiguous genomic target. Experimental guide selection
requires independent genome-coordinate and exon-structure validation.

### 4.4 `GOAnalysis`

File: `src/modules/go_analysis.py`

Responsibilities:

- build missing GO helper images;
- add BP, MF, CC, and combined GO columns through the UniProt REST API;
- compute GO-KFI by generated `unique_id`;
- extract fully targetable GO terms;
- compute parent-term summaries;
- create a compact GO-to-gene reverse lookup; and
- write `go_analysis_results.txt`.

New runs send `normalized_genes.tsv` directly to GO annotation and use its
validated non-negative `guide_count`; they do not reparse serialized guide
payloads. Runs without a normalized current pointer retain the legacy
contig-table fallback. In either form, each `unique_id` is counted once:

Normalized runs require the GO annotation and GO-KFI helper images to declare
normalized schema version `1`. Cached helpers without the label are rebuilt
from the local Docker contexts before execution.

```text
GO-KFI = unique genes with retained guides / unique genes in the GO term
```

Classification values are `Fully`, `Partially`, and `Not`.

GO Slim summaries assign each fully targetable term to one `is_a` category.
Selection uses the minimum ancestor edge distance and breaks equal-distance
ties by ascending GO ID. Summary rows use descending count and ascending GO ID.

Although standard annotation accepts custom FASTA identifiers, GO lookup
requires best-hit subject identifiers that can be interpreted as UniProtKB
accessions.

### 4.5 `PathwayAnalysis`

File: `src/modules/pathway_analysis.py`

Responsibilities:

- parse GPML 2013a `GeneProduct` nodes;
- retain only configured database/identifier pairs;
- match exact annotation identifiers and exact SGD FASTA cross-references;
- cross-validate explicit annotation Ensembl and SGD identities, quarantining
  conflicting subjects from every pathway index;
- collapse duplicate nodes and exact cross-database aliases into canonical gene
  units;
- fail closed for ambiguous or incomplete cross-database identity;
- enforce pathway coverage policy;
- compute PW-KFI; and
- write result and audit artifacts.

Pathway analysis prefers the current normalized generation and reconstructs
raw guide strings only at its compatibility boundary. If no current pointer
exists it reads the legacy table. A present but invalid pointer fails closed
instead of silently mixing generations.

Non-evaluable GPML nodes are excluded from the denominator. Exact canonical
matched and unmatched units are included:

```text
PW-KFI = canonical pathway genes with retained guides
         / evaluable canonical pathway genes
```

Pathways with incomplete cross-database canonicalization are skipped. Low
coverage pathways are skipped only when `skip_low_coverage_pathways` is true.

### 4.6 `ReportGenerator`

File: `src/modules/report_generator.py`

Responsibilities:

- collect standard workflow results;
- publish the report-source snapshot available for the run;
- calculate report statistics;
- generate `summary.txt` when enabled;
- maintain a minimal, integrity-checked SQLite GO metadata cache;
- generate `web_data/summary.json`;
- generate static assets;
- stream immutable browser payloads under
  `web_data/report_payloads/<namespace>-<generation>/`; and
- atomically publish `report.html` when enabled after its payload generation is
  complete; and
- display a prominent guide-design scope warning before report statistics.

`report.html` is a rendering shell and does not embed the analysis payload.
Browser data is published as classic `core.js` and `guides/*.js` scripts so the
report works when opened directly with `file://`, without module loading or
runtime `fetch`. Core data loads with the shell. Guide chunks load only when the
Guide section is opened. Python writes guide chunks by streaming source rows,
publishes the immutable payload generation, and only then atomically replaces
the HTML shell.

Published payload generations remain immutable and are not automatically
deleted. This preserves delayed guide-chunk reads in report tabs that were
opened before a rerun. Unreferenced generations may be removed only when no
reader still uses them. A generation without a matching marker under
`report_payloads/.committed/` is an orphan from a failed HTML publication and
is removed after the next successful commit.

All external strings are context-escaped before HTML rendering. Generated data
scripts use script-safe serialization, and dynamic report tables never place
unescaped external strings into `innerHTML` or inline event-handler source.

The HTML and text summary identify a no-genome run as preliminary-only. When a
genome was configured, they state that it was used for off-target searching but
that candidates remain transcriptome-derived and may cross splice junctions.

QuickGO supplies missing GO labels and definitions. Native analysis tables are
not loaded into SQLite.

### 4.7 `PathwayReportGenerator`

File: `src/modules/pathway_report_generator.py`

Responsibilities:

- require and reconcile `pathway_kfi_table.csv` with
  `pathway_gene_mapping.tsv`;
- reject a canonical mapping unit whose non-empty annotation Ensembl ID is not
  one of its Ensembl pathway xrefs;
- publish those files and referenced GPML files in the pathway snapshot;
- reuse the standard report shell and assets;
- render pathway SVG content from GPML graphics/interactions;
- publish pathway-specific guide exports in the external browser payload;
- generate `web_data/pathway_summary.json`; and
- atomically generate one rendering shell named by
  `pathway_report.output_name` after publishing its payload generation.

Because the pathway report reuses the standard report shell, it displays the
same no-genome preliminary warning or genome-backed transcriptome-design
limitation before its statistics.

GPML labels and identifiers are escaped in generated SVG attributes and text.
The browser validates payload-provided SVG against the report's element and
attribute allowlist before inserting it into the document.

The Pathway HTML is a browser-rendering shell. Its immutable
`web_data/report_payloads/<namespace>-<generation>/core.js` supplies the runtime
projection of summary metadata/statistics, displayed pathway rows, and generated
SVG and pathway-guide data. The shared standard-report guide table remains in
deferred `guides/*.js` chunks and loads only when the Guide section is opened.
Raw `matched_genes` and pathway mapping rows remain regeneration sources rather
than browser payloads.
`web_data/pathway_summary.json` retains the complete machine-readable summary
artifact independently of the browser projection and is not required for
runtime rendering.

The pathway report does not create a pathway SQLite database.

## 5. Artifact contracts

### 5.1 Standard GO artifacts

Key intermediate files include:

- `transcriptome.fa`
- `query.pep`
- `annotation_table.txt`
- `annotation_classified.txt`
- `normalized_gene_guide_current.json`
- `.normalized_gene_guide_generations/<uuid>/normalized_genes.tsv`
- `.normalized_gene_guide_generations/<uuid>/normalized_guides.tsv`
- `.normalized_gene_guide_generations/<uuid>/normalized_gene_guide_manifest.json`
- `.normalized_gene_guide_generations/<uuid>/annotation_with_guides_final.tsv`
- `annotation_with_guides_go.tsv`
- `go_kfi_table.csv`
- `go_kfi_stat.csv`
- `go_kfi_table_fully_targetable.csv`
- `go_reverse_lookup.tsv`
- `go_parent_counts.csv`

Key run/report files include:

- `twheel.log`
- `summary.txt` when enabled
- `report.html` when enabled
- `static/`
- `web_data/summary.json`
- `web_data/report_payloads/<namespace>-<generation>/core.js`
- `web_data/report_payloads/<namespace>-<generation>/guides/*.js`
- `web_data/report_payloads/<namespace>-<generation>/manifest.json`
- `web_data/report_payloads/.committed/<namespace>-<generation>`
- `web_data/twheel_data.db`
- `web_data/report_generation_results.txt`

### 5.2 Pathway artifacts

Pathway analysis always writes non-empty header-bearing contract files even
when no pathway is reportable:

- `pathway_kfi_table.csv`
- `pathway_kfi_fully_targetable.csv`
- `pathway_gene_mapping.tsv`
- `pathway_gene_match_audit.tsv`
- `pathway_annotation_identifier_audit.tsv`
- `pathway_analysis_results.txt`

The report adds:

- `<pathway_report.output_name>`
- `web_data/pathway_summary.json`
- shared `static/` assets
- `web_data/report_payloads/<namespace>-<generation>/core.js`
- shared-guide chunks in
  `web_data/report_payloads/<namespace>-<generation>/guides/*.js`
- `web_data/report_payloads/<namespace>-<generation>/manifest.json`
- `web_data/report_payloads/.committed/<namespace>-<generation>`

After a successful HTML commit, `summary.json` or `pathway_summary.json`
records a `browser_payload` pointer to the generation manifest. The pointer is
metadata; the HTML script path remains the browser's runtime entry point.

Both report HTML files require the adjacent `static/` and `web_data/`
directories when moved or archived. Their classic payload scripts support
direct `file://` opening without a web server.

## 6. Report-source storage

File: `src/modules/report_source_store.py`

Report sources are stored under:

```text
web_data/report_sources/
├── manifest.json
└── snapshots/<namespace>-<uuid>/...
```

The manifest records relative paths, sizes, and SHA-256 hashes. A namespace is
published by creating an independent copy or filesystem clone, hashing it,
renaming the completed snapshot, and atomically replacing the manifest under a
file lock. Unreferenced snapshot directories are pruned.

The standard namespace publishes the report inputs that exist. When all three
normalized gene/guide artifacts are present, it snapshots them and omits the
larger legacy compatibility table. The normalized gene table is cached once;
the guide table is streamed for statistics and into bounded browser payload
chunks without retaining a second all-guide projection. Report-source snapshots
are inputs for regeneration and are not loaded by the browser. Older runs
continue to stage and read the legacy table. The pathway namespace requires the
canonical table/mapping pair and additionally stages referenced GPML files.

Generator classes may load the last verified namespace when their configured
`intermediate` directory does not exist. This is not currently exposed as a
reliable CLI restore path because `tWhEELPipeline` recreates `intermediate`
during initialization.

## 7. SQLite cache schema

`web_data/twheel_data.db` is a GO display cache, not the analysis database.
The current schema contains only:

```sql
CREATE TABLE report_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
) WITHOUT ROWID;

CREATE TABLE go_terms (
    go_id TEXT PRIMARY KEY,
    go_name TEXT NOT NULL,
    go_definition TEXT NOT NULL
) WITHOUT ROWID;
```

The cache is rebuilt atomically when its schema, metadata, GO identifier set,
or content integrity is invalid. `database.force_rebuild` requests a refresh,
while usable cached labels remain fallback data if an API refresh fails.

## 8. CLI metadata

The implemented options are:

- `--config`, `-c`
- `--step`, `-s`
- `--name`, `-n`
- `--list-steps`
- `--version`, `-v`

`src/twheel_metadata.py` is the canonical runtime metadata source. It defines
`__version__`, the eight public step names, their workflow assignments, and their
help summaries. The parser, `--help`, `--list-steps`, and `--step` choices all
consume this shared definition. The current version is 1.3.1, and
`--version` prints `tWhEEL version 1.3.1`.

Help, version, and step-listing requests complete before YAML or analysis-class
imports. This keeps the CLI metadata inspectable without a fully installed
pipeline runtime.

An unknown `--step` value is rejected by `argparse` before configuration loading
or output-directory creation. A valid single-step invocation still bypasses
full-run input validation and never executes prerequisites automatically.

## 9. Development and extension rules

- Run from the repository root because local Docker build contexts are relative
  to it. Runtime data paths are translated from an absolute `intermediate`
  workspace and do not depend on the current working directory.
- Use Python 3.12 for the documented environment.
- Keep code comments, docstrings, configuration comments, and manuals in
  English.
- Update `src/twheel_metadata.py` once when preparing a release, and create a
  matching Git tag. Do not add a second hard-coded runtime version.
- Add a module class with a `run()` method, import it in `src/twheel.py`, add an
  explicit `run_step()` branch, add one `StepDefinition` to the canonical CLI
  metadata, define its config and artifact contracts, and add tests.
- Modules are independent classes; there is no `BaseModule` abstraction.
- Use `logging` for module progress and failures.
- Preserve native large-result formats and publish only report inputs needed by
  the report-source store.
- Use pytest for unit and contract tests.

The test suite covers parsers, input-mode consistency, cache invalidation,
guide filtering, GO/PW-KFI counting, exact pathway identifier reconciliation,
report-source integrity, SQLite cache integrity, and report generation. Docker
and remote-service end-to-end execution remains a separate integration concern.

## 10. Known implementation limitations

- The pathway annotation profile is SGD/*S. cerevisiae*-centered.
- The two annotation profiles overwrite common filenames and should use
  separate runs.
- Execution is sequential; `execution.parallel` is unused.
- Docker cannot be disabled through the current configuration.
- Local helper image builds use plain `docker build` and require daemon access
  for the invoking user; tWhEEL never invokes `sudo` itself.
- BUSCO is forced to the `linux/amd64` platform.
- Three helper/alignment Dockerfiles still use the end-of-life Python 3.9 base.
- Base-image digests and helper Python dependency sets are not fully locked, so
  rebuilding an existing helper tag is not byte-for-byte reproducible.
- The default BUSCO and hard-coded SeqKit images are legacy compatibility pins,
  pending explicit migration and regression testing.
- Report snapshot fallback is tested at the generator-class level but is not a
  working public CLI restore flow after `intermediate` deletion.
- Even genome-backed runs can contain transcriptome-derived candidate guides
  that cross splice junctions; tWhEEL does not yet perform genomic target-site
  reconciliation.
- Portable configuration templates contain explicit absolute-path placeholders
  and are not runnable until those paths and the BUSCO lineage are replaced.

## 11. External services and tools

- TransDecoder
- BUSCO
- SeqKit
- FASTA/ggsearch36
- CrisFlash
- MAFFT
- UniProt REST API
- QuickGO API
- Gene Ontology OBO/GO Slim resources
- WikiPathways GPML

## Document information

- **Status:** Current implementation reference
- **Reviewed:** 2026-07-17
- **Version authority:** `src/twheel_metadata.py`; release Git tags are expected
  to match its `__version__` value
