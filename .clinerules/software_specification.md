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
| `gene_prediction` | `transdecoder`, `busco` |
| `annotation` or `sgd_annotation` | `fasta_search` |
| `guide_design` | `crispr` |
| `go_analysis` | `go_analysis` |
| `report_generation` | `report` |
| `pathway_analysis` | `pathway_analysis` |

`pathway_report` is optional; `pathway_report_generation` defaults its HTML
name to `pathway_report.html`. `database` is optional and is read only by the
standard report generator. A pathway-only profile does not read
`go_analysis`, `report`, or `database`.

For a full orchestrated run, including `pathway_report_generation` in
`execution.steps` also makes
`pathway_analysis.wikipathways_gpml_dir` an input-validation requirement.

### 3.2 Input keys

| Key | Contract |
| --- | --- |
| `input.transcriptome` | Required for every full run. Used in both protein input modes. |
| `input.protein_sequences` | Optional existing protein FASTA. |
| `input.protein_transcript_mapping` | Required by annotation modules when existing protein mode is effective. |
| `input.genome_assembly` | Optional off-target reference; transcriptome is the fallback. |
| `input.protein_db_for_goanalysis` | Required by full-run validation when `annotation`, `go_analysis`, or `report_generation` is configured. |
| `input.sgd_db` | Required when `sgd_annotation` is configured. Also used for exact pathway cross-reference resolution. |

The existing-protein mapping is a headerless TSV. Blank lines and lines whose
first character is `#` are skipped. The first two columns are:

```text
transcript_id<TAB>protein_id
```

Extra columns are ignored. The loader stores `protein_id -> transcript_id`.

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
- `pathway_report.output_name` is effective and defaults to
  `pathway_report.html`.
- `database.force_rebuild` controls rebuilding the minimal GO cache.
- `execution.verbose` controls the logging level.
- `execution.parallel` is currently ignored; execution is sequential.
- `docker.max_memory` and `docker.use_docker` are currently ignored.
- `report.web_port`, `pathway_report.web_port`, and
  `pathway_report.generate_html` are currently ignored. No web server is
  started.

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

Multi-group alignment builds and runs the local Docker image
`run_alignment_and_mask/v1.0`. Every expected group must produce exactly one
non-empty consensus FASTA with the expected identifier.

### 4.3 `GuideDesign`

File: `src/modules/guide_design.py`

Responsibilities:

- copy or decompress the genome assembly, or fall back to the transcriptome;
- convert supported IUPAC degenerate bases to `N`;
- prepare Single sequences and validated Multi consensus sequences;
- run CrisFlash separately for the available Single/Multi groups;
- retain site-specific targets; and
- merge guide data into `annotation_with_guides_final.tsv`.

CrisFlash target field four has the form
`<sequence>/<MM0>/<MM1-2>`. A site-specific target is retained only when
`MM0 == 1` and `MM1-2 == 0`. Other CrisFlash scores are passed through without
a tWhEEL-enforced range.

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

Annotation rows are collapsed by `unique_id`; GO terms and guide entries are
de-duplicated before counting:

```text
GO-KFI = unique genes with retained guides / unique genes in the GO term
```

Classification values are `Fully`, `Partially`, and `Not`.

Although standard annotation accepts custom FASTA identifiers, GO lookup
requires best-hit subject identifiers that can be interpreted as UniProtKB
accessions.

### 4.5 `PathwayAnalysis`

File: `src/modules/pathway_analysis.py`

Responsibilities:

- parse GPML 2013a `GeneProduct` nodes;
- retain only configured database/identifier pairs;
- match exact annotation identifiers and exact SGD FASTA cross-references;
- collapse duplicate nodes and exact cross-database aliases into canonical gene
  units;
- fail closed for ambiguous or incomplete cross-database identity;
- enforce pathway coverage policy;
- compute PW-KFI; and
- write result and audit artifacts.

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
- generate static assets; and
- generate `report.html` when enabled.

QuickGO supplies missing GO labels and definitions. Native analysis tables are
not loaded into SQLite.

### 4.7 `PathwayReportGenerator`

File: `src/modules/pathway_report_generator.py`

Responsibilities:

- require and reconcile `pathway_kfi_table.csv` with
  `pathway_gene_mapping.tsv`;
- publish those files and referenced GPML files in the pathway snapshot;
- reuse the standard report shell and assets;
- render pathway SVG content from GPML graphics/interactions;
- embed pathway-specific guide exports;
- generate `web_data/pathway_summary.json`; and
- generate one HTML file named by `pathway_report.output_name`.

The pathway report does not create a pathway SQLite database.

## 5. Artifact contracts

### 5.1 Standard GO artifacts

Key intermediate files include:

- `transcriptome.fa`
- `query.pep`
- `annotation_table.txt`
- `annotation_classified.txt`
- `annotation_with_guides_final.tsv`
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

The standard namespace publishes the report inputs that exist. The pathway
namespace requires the canonical table/mapping pair and additionally stages
referenced GPML files.

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

Two metadata limitations remain in the current source:

- `--list-steps` lists only the five GO-profile steps even though dispatch also
  accepts `sgd_annotation` and the two pathway steps.
- `--version` returns a hard-coded `1.0.0` that is not synchronized with Git
  tags.

## 9. Development and extension rules

- Run from the repository root because several Docker build contexts and data
  mounts are relative to the current working directory.
- Use Python 3.12 for the documented environment.
- Keep code comments, docstrings, configuration comments, and manuals in
  English.
- Add a module class with a `run()` method, import it in `src/twheel.py`, add an
  explicit `run_step()` branch, add the step to CLI metadata, define its config
  and artifact contracts, and add tests.
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
- Local helper image builds currently invoke `sudo docker build`.
- BUSCO is forced to the `linux/amd64` platform.
- Report snapshot fallback is tested at the generator-class level but is not a
  working public CLI restore flow after `intermediate` deletion.
- CLI step and version metadata are stale as described above.
- Dated configuration files reference Git-ignored `LOCAL_DATA` and are local run
  records rather than portable examples.

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
- **Reviewed:** 2026-07-16
- **Version authority:** Source behavior and Git history; the hard-coded CLI
  version is currently stale
