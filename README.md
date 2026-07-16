# tWhEEL

**Transcriptome Workflow for Hierarchical Exploration and Expression Levels**

tWhEEL is a Docker-oriented transcriptome analysis pipeline for protein
prediction, similarity-based annotation, SpCas9 guide design, and knockout
feasibility assessment. The current implementation has two report workflows:

```text
GO workflow
Transcriptome -> Gene prediction -> Protein annotation -> Guide design
              -> GO analysis -> GO report

Pathway workflow
Transcriptome -> Gene prediction -> SGD annotation -> Guide design
              -> Pathway analysis -> Pathway report
```

The pathway workflow is currently centered on *Saccharomyces cerevisiae*
because its annotation step uses an SGD reference database. Run the GO and
pathway workflows in separate output directories: `annotation` and
`sgd_annotation` intentionally write the same intermediate filenames.

## Requirements

- Linux or macOS
- Python 3.12 (the currently used and documented environment)
- A running Docker daemon and the Docker CLI
- Network access when Docker images, BUSCO data, GO ontology files, or QuickGO
  records are not already cached
- Active network access to `https://rest.uniprot.org` during every
  `go_analysis` annotation run; these UniProt responses are not cached across
  runs
- RAM and storage appropriate for the input data and selected Docker images

The annotation implementations build a content-tagged local multi-gene
alignment image only when the matching build context is not already available.
GO analysis still invokes `sudo docker build` when one of its three configured
helper images is missing. The user running tWhEEL must therefore be able to
execute the required Docker commands in the local setup.

## Installation

```bash
git clone https://github.com/KazukiNakamae/tWhEEL.git
cd tWhEEL

conda create -n twheel python=3.12 -y
conda activate twheel
pip install -r requirements.txt
```

The GO helper images may be built before the first run:

```bash
docker build -t kazukinakamae/add-go-columns:v1.2 \
  src/docker/add_go_columns/
docker build -t kazukinakamae/compute_go_kfi:v1.2 \
  src/docker/compute_go_kfi/
docker build -t kazukinakamae/count_go_parent_terms:v1.1 \
  src/docker/count_go_parent_terms/
```

The multi-gene alignment image is built automatically. Its local tag includes
a hash of the Dockerfile, requirements, and alignment script, so unchanged
runs reuse the existing image instead of rebuilding it.

Other configured images (TransDecoder, BUSCO, SeqKit, FASTA3, and CrisFlash)
are pulled by Docker on first use if they are not available locally.

Run commands from the repository root because local Docker build contexts are
resolved from it. Runtime analysis files are mounted from the absolute
`intermediate` directory, so `output.base_dir` may be relative or absolute and
may be outside the repository. Docker Desktop users must allow file sharing for
an absolute output location that is not already shared with Docker.

## Quick start: GO workflow

Copy the standard template and edit the paths before running it. The unedited
`config/default.yaml` cannot run because its required input paths are empty.

```bash
cp config/default.yaml config/my_analysis.yaml
```

At minimum, set:

```yaml
input:
  transcriptome: "/absolute/path/to/transcriptome.fasta"
  protein_db_for_goanalysis: "/absolute/path/to/uniprot_sprot.fasta"

output:
  base_dir: "outputs"
```

Then run the steps listed in `execution.steps`:

```bash
python src/twheel.py \
  --config config/my_analysis.yaml \
  --name my_analysis
```

This creates a directory such as
`outputs/tWhEEL_on_YYYYMMDD_HHMMSS_my_analysis/`.

### Protein database requirement

`protein_db_for_goanalysis` is the implemented key; `uniprot_db` is not read by
the pipeline. The annotation parser accepts opaque/custom FASTA identifiers,
but the GO helper extracts a UniProtKB accession from each best-hit
`subject_id` and queries the UniProt REST API. A complete GO workflow therefore
requires best-hit identifiers that UniProt can resolve. A custom database with
unresolvable identifiers may be usable for annotation alone but will not
produce usable GO annotations.

The GO annotation helper requires an active network connection at runtime. It
tries each failed UniProt batch up to three times, then stops without
publishing a partially annotated table. If the error says that UniProt cannot
be reached, ensure both the host and Docker can access
`https://rest.uniprot.org`, connect to the internet if needed, and rerun only
`go_analysis` in the same run directory:

```bash
python src/twheel.py \
  --config config/my_analysis.yaml \
  --step go_analysis \
  --name tWhEEL_on_YYYYMMDD_HHMMSS_my_analysis
```

## Existing protein mode

Provide a protein FASTA and enable `skip_if_protein_provided` to skip
TransDecoder. The transcriptome is still required for sequence extraction and
guide design.

```yaml
input:
  transcriptome: "/absolute/path/to/transcripts.fasta"
  protein_sequences: "/absolute/path/to/proteins.fasta"
  protein_transcript_mapping: "/absolute/path/to/mapping.tsv"
  protein_db_for_goanalysis: "/absolute/path/to/uniprot_sprot.fasta"

transdecoder:
  skip_if_protein_provided: true
```

The mapping file is a headerless TSV. Lines beginning with `#` and blank lines
are ignored; extra columns are ignored.

```text
transcript_id<TAB>protein_id
```

The first two columns are read in that order, and the implementation constructs
a `protein_id -> transcript_id` lookup. Header rows are not allowed.

Before a FASTA search starts, the mapping is checked against the prepared
`intermediate/query.pep` and `intermediate/transcriptome.fa` files. Mapping
protein IDs must exactly match all query protein IDs, and every mapped
transcript ID must name a non-empty transcriptome FASTA record. Additional
transcriptome records are allowed. Malformed rows, empty IDs, conflicting
mappings for one protein ID, duplicate FASTA IDs, missing mappings, and missing
mapped transcripts stop annotation with an error. These failures are not
interpreted as genes without targetable guides.

Existing protein mode is selected only when the configured protein file exists
and `skip_if_protein_provided` is true (the option defaults to true when
omitted). Otherwise, the pipeline uses TransDecoder mode. See
`config/example_existing_protein.yaml` for a complete template.

## Pathway workflow

The pathway workflow replaces `annotation` with `sgd_annotation` and requires a
local WikiPathways GPML directory. A minimal configuration outline is:

```yaml
input:
  transcriptome: "/absolute/path/to/transcripts.fasta"
  protein_sequences: "/absolute/path/to/proteins.fasta"
  protein_transcript_mapping: "/absolute/path/to/mapping.tsv"
  sgd_db: "/absolute/path/to/sgd_reference.fasta"
  genome_assembly: "/absolute/path/to/genome.fasta"  # optional

pathway_analysis:
  wikipathways_gpml_dir: "/absolute/path/to/gpml_directory"
  target_gene_databases:
    - database_name: "SGD"
      id_pattern: "^S\\d{9}$"
      annotation_key: "sgd_id"
    - database_name: "Ensembl"
      id_pattern: "^Y[A-Z]{2}\\d{3}[CW]$"
      annotation_key: "ensembl_id"
  minimum_pathway_coverage: 0.5
  skip_low_coverage_pathways: true

pathway_report:
  output_name: "pathway_report.html"

execution:
  steps:
    - gene_prediction
    - sgd_annotation
    - guide_design
    - pathway_analysis
    - pathway_report_generation
  verbose: true
```

`database_name` must match the GPML `Xref` database name exactly, and
`id_pattern` is applied as a full regular-expression match. The dated pathway
configuration in `config/input_for_PathwayAnalysis_20250917.yaml` is a local
run record: it references `LOCAL_DATA`, which is excluded from Git and is not
available in a fresh clone.

When the SGD reference provides an exact Ensembl-to-SGD cross-reference, each
annotation record that contains both an Ensembl ID and an explicit SGD identity
must agree with that cross-reference. A conflicting record is quarantined from
all pathway indexes, contributes no guide to PW-KFI, and is recorded as
`conflicting_annotation_cross_reference` in both pathway audit artifacts. Empty
legacy Ensembl fields and custom subjects without an explicit SGD identity
remain supported. Pathway report generation also rejects stale mapping
artifacts whose non-empty annotation Ensembl ID disagrees with their Ensembl
pathway xref.

## Command-line usage

```bash
python src/twheel.py --help
python src/twheel.py --config config/my_analysis.yaml --name my_analysis
python src/twheel.py --config config/my_analysis.yaml \
  --step annotation \
  --name tWhEEL_on_YYYYMMDD_HHMMSS_my_analysis
```

Supplying a simple `--name` creates a new timestamped run. Supplying a complete
name beginning with `tWhEEL_on_` reuses that run directory. Reuse is not an
automatic resume mechanism: tWhEEL does not detect completed steps. Without
`--step`, every configured step runs again in order. A single-step run must
reuse the correct directory and its prerequisite intermediate files. Within
`annotation` and `sgd_annotation`, validated Multi-gene consensus results are
an exception: they are cached per input group and reused on subsequent runs.

The implemented step names and their main contracts are:

| Step | Main prerequisite | Main output |
| --- | --- | --- |
| `gene_prediction` | Transcriptome; optional existing proteins | `query.pep`, BUSCO output |
| `annotation` | `query.pep`, transcriptome, GO protein DB | `annotation_classified.txt` |
| `sgd_annotation` | `query.pep`, transcriptome, SGD DB | `annotation_classified.txt` |
| `guide_design` | Classified Single/Multi sequences | Versioned normalized gene/guide generation |
| `go_analysis` | Normalized genes, or a legacy annotation/guide table | GO-KFI tables and GO lookup files |
| `pathway_analysis` | Normalized genes/guides or a legacy table, plus GPML files | PW-KFI, mappings, and audit tables |
| `report_generation` | GO workflow artifacts | `report.html`, optional `summary.txt` |
| `pathway_report_generation` | Canonical pathway table and mapping | Configured pathway HTML report |

`--step` dispatches all eight names above. The current `--list-steps` output is
legacy CLI metadata and lists only the five standard GO-workflow steps. The
current `--version` output is also a legacy `1.0.0` string and is not synchronized
with the latest Git tag or commits.

## Input behavior

- `input.transcriptome` is required for a full pipeline run in both protein
  input modes.
- Transcriptome, protein, genome, and annotation database inputs may be plain
  FASTA or gzip-compressed FASTA; format handling is based on a `.gz` suffix.
- The implementation does not validate a fixed FASTA extension or an ACGTN-only
  alphabet. During off-target reference preparation, supported IUPAC degenerate
  bases are converted to `N` for CrisFlash compatibility.
- `input.genome_assembly` is optional. When absent, guide design uses the
  transcriptome as its off-target reference and logs that this is less
  comprehensive.
- YAML duplicate keys are rejected rather than silently overwritten.

## Configuration behavior

`config/default.yaml` covers the GO workflow.
`config/input_for_PathwayAnalysis_20250917.yaml` shows the additional pathway
sections, but its input paths are local and must be replaced.

The following current behaviors are easy to miss:

- `output.base_dir` is the output root. At startup, the pipeline replaces
  `output.intermediate` with `<run>/intermediate`; the YAML value is not an
  independent output location.
- Steps always run sequentially. `execution.parallel` is currently not read.
- `multi_gene_alignment.workers` bounds concurrent MAFFT processes inside the
  alignment container. `multi_gene_alignment.memory` is the total container
  memory limit, and `multi_gene_alignment.cache_enabled` controls the validated
  per-group consensus cache.
- The alignment cache is durable and has no automatic eviction. When no run is
  active, removing `intermediate/multi_genes_alignment_cache` safely reclaims
  its space; the next annotation run recomputes missing groups.
- Guide design publishes immutable normalized generations. The small
  `normalized_gene_guide_current.json` pointer selects the complete generation;
  a failed rerun does not replace the previous pointer. Source hashes make a
  rerun of annotation invalidate stale normalized data until `guide_design` is
  rerun.
- Normalized GO input requires helper images that declare schema support. The
  bundled `v1.2` Dockerfiles carry this label; cached legacy/custom helpers
  without it are rebuilt from the local helper contexts.
- `report.generate_html` and `report.generate_summary` control standard report
  outputs.
- `report.guide_chunk_rows` sets the positive number of shared guide rows per
  deferred browser payload script and defaults to `5000`.
- `pathway_report.output_name` controls the pathway HTML filename.
- `report.web_port`, `pathway_report.web_port`,
  `pathway_report.generate_html`, `docker.max_memory`, and
  `docker.use_docker` are retained in older configurations but are not read by
  the current Python implementation.
- `database.force_rebuild` controls rebuilding the minimal GO metadata cache.

## Output layout

Depending on the selected workflow and report settings, a run has this shape:

```text
<output.base_dir>/tWhEEL_on_<timestamp>[_<name>]/
├── twheel.log
├── intermediate/
│   ├── transcriptome.fa
│   ├── query.pep
│   ├── annotation_classified.txt
│   ├── multi_genes_alignment_cache/             # validated MAFFT group cache
│   ├── normalized_gene_guide_current.json
│   ├── .normalized_gene_guide_generations/
│   │   └── <uuid>/
│   │       ├── normalized_genes.tsv              # one row per unique_id
│   │       ├── normalized_guides.tsv             # one row per typed guide
│   │       ├── normalized_gene_guide_manifest.json
│   │       └── annotation_with_guides_final.tsv  # compact legacy view
│   ├── go_kfi_table.csv                         # GO workflow
│   ├── pathway_kfi_table.csv                    # pathway workflow
│   ├── pathway_gene_mapping.tsv                 # pathway workflow
│   ├── pathway_gene_match_audit.tsv             # pathway workflow
│   └── pathway_annotation_identifier_audit.tsv  # pathway workflow
├── report.html                    # standard report when enabled
├── summary.txt                    # standard text summary when enabled
├── pathway_report.html            # default pathway report name
├── static/                        # CSS, JavaScript, logo, report images
└── web_data/
    ├── summary.json               # GO workflow
    ├── pathway_summary.json       # pathway workflow
    ├── twheel_data.db             # GO label/definition cache only
    ├── report_generation_results.txt
    ├── report_payloads/
    │   ├── .committed/<namespace>-<generation>  # HTML commit marker
    │   └── <namespace>-<generation>/
    │       ├── core.js                 # immutable core browser data
    │       ├── guides/*.js             # deferred guide chunks
    │       └── manifest.json           # sizes and SHA-256 checksums
    └── report_sources/
        ├── manifest.json
        └── snapshots/
```

Most analysis results remain native TSV, CSV, FASTA, BED, text, or GPML files;
they are not copied into a result database. `web_data/twheel_data.db` contains
only report metadata and cached GO term labels/definitions. The pathway report
does not create a separate SQLite result database.

Report generators atomically publish versioned, checksum-verified snapshots of
the report inputs available to them. For normalized runs, the standard report
snapshot contains `normalized_genes.tsv`, `normalized_guides.tsv`, and their
manifest, and deliberately omits the legacy compatibility table. These three
files are validated together. The gene table is cached once; guides are
streamed without full normalized-row materialization for statistics and are
written as bounded guide payload chunks. Older runs retain the legacy-table
fallback. The pathway report requires a complete canonical KFI-table/mapping
pair. Report-source snapshots under `web_data/report_sources/` are regeneration
inputs; they are not browser payloads. The generator classes can fall back to a
published source snapshot when their configured `intermediate` directory does
not exist. The public CLI currently recreates `intermediate` during
initialization, so this fallback is an internal generator contract, not a CLI
backup/restore procedure.

`report.html` and `pathway_report.html` are browser-rendering shells and do not
contain the analysis payload inline. Each shell loads immutable data from
`web_data/report_payloads/<namespace>-<generation>/core.js` and, when needed,
`guides/*.js`. These files use classic scripts instead of module scripts or
runtime `fetch` calls, so the reports continue to work when opened directly with
`file://`. Core data is loaded when the shell opens; guide chunks are deferred
until the Guide section is opened. Python streams guide rows into bounded chunks,
publishes the complete payload generation first, and then atomically replaces
the HTML shell so a report never points at a partial payload generation.
For the Pathway report, the pathway table, SVGs, and pathway-specific guide
exports live in its external `core.js`; the shared standard-report guide table
uses the same deferred chunk contract as `report.html`.

Published payload generations are not removed automatically. A report tab may
read its core payload immediately and request guide chunks from that same
generation much later. After all report tabs are closed, unreferenced older
generation directories may be removed manually; keep the generations named by
the current HTML script paths and summary JSON pointers.
Generations that never reached an HTML commit have no corresponding marker in
`report_payloads/.committed/` and are removed after the next successful
publication.

When moving, archiving, or saving a report, keep both the adjacent `static/` and
`web_data/` directories with the HTML file. `web_data/summary.json` and
`web_data/pathway_summary.json` remain complete machine-readable output
artifacts; browser rendering uses the versioned JavaScript payload generation
rather than reading either summary JSON file at runtime. Raw `matched_genes` and
pathway mapping rows remain source artifacts rather than browser payloads.
On successful publication, each summary JSON records a `browser_payload`
pointer with the schema version, namespace, generation ID, and manifest path.

Report inputs such as annotation text, FASTA identifiers, configuration paths,
GO labels, and GPML labels are treated as untrusted display text. Generated data
scripts use script-safe serialization, dynamic tables use text nodes or explicit
HTML escaping, and generated pathway SVG is validated before insertion. These
protections preserve the original strings in browser payloads and display them
as text while preventing them from becoming executable HTML or JavaScript.

## Result semantics

### Single and Multi classifications

- **Single**: exactly one query sequence has a given best-hit `subject_id`.
- **Multi**: two or more query sequences share that best-hit `subject_id`.

These are similarity-hit groups, not an independent proof of biological gene or
transcript redundancy. Each group receives a generated `unique_id`, which is
the downstream gene-level counting unit. `annotation_classified.txt` retains
the contig-level provenance. Guide design then writes one normalized gene row
per `unique_id` and one normalized guide row per distinct typed guide tuple.
The generation-local legacy TSV is a compact one-row-per-gene compatibility
view, not the primary source for new GO, pathway, or report processing.

### Site-specific guides

CrisFlash serializes the guide sequence and mismatch counts in the fourth BED
field as `<sequence>/<MM0>/<MM1-2>`. The current filter retains a guide only when
`MM0 == 1` and `MM1-2 == 0`. Any additional score is passed through from
CrisFlash; tWhEEL does not enforce a documented 0-1 range.

### GO-KFI

For each GO term, the normalized gene row supplies a validated `guide_count`,
and each `unique_id` is counted once. Legacy runs still collapse annotation
rows and de-duplicate serialized guide targets before applying the same rule:

```text
GO-KFI = unique genes with at least one retained guide
         / unique genes annotated to the GO term
```

Scores are classified as Fully (`1.0`), Partially (`0 < score < 1.0`), or Not
targetable (`0.0`).

The GO Slim parent summary assigns each fully targetable GO term to exactly one
`is_a` ancestor. The input term itself is preferred, otherwise the nearest GO
Slim ancestor is selected by edge distance. Equal-distance candidates are
resolved by ascending GO ID. Equal-count summary rows are also ordered by GO
ID, so a fixed ontology and input produce the same assignment and table order.

### PW-KFI

Pathway gene nodes are filtered by configured database/identifier contracts and
collapsed to exact canonical gene units. Non-evaluable nodes are excluded from
the denominator; exact canonical matched and unmatched units are included:

```text
PW-KFI = canonical pathway genes with at least one retained guide
         / evaluable canonical pathway genes
```

Pathways with unresolved cross-database canonical identity are skipped. Low
annotation coverage is either reported or skipped according to
`minimum_pathway_coverage` and `skip_low_coverage_pathways`. The two pathway
audit TSV files record matching, exclusion, and identifier issues. Conflicting
Ensembl/SGD identities are fail-closed and never transfer guides between genes.

## Testing

With the development dependencies from `requirements.txt` installed:

```bash
python -m pytest -q
```

The unit suite exercises parsing, input-mode resolution, guide filtering,
GO/PW-KFI counting, pathway identifier reconciliation, report generation, and
report-source integrity. It does not replace an end-to-end run of all external
Docker images and network services.

## Troubleshooting

### A single step cannot find its input

Use the complete existing run name with `--name`. Omitting it creates a new,
empty timestamped run directory.

### Docker build or permission failure

Confirm that the daemon is running and that `docker` works for the current
user. Multi-gene alignment uses a user-writable temporary buildx state and does
not invoke `sudo`; GO helper-image fallback builds still use `sudo docker
build`. Linux group-management commands such as `usermod` do not apply to
Docker Desktop on macOS.

### Workspace lock timeout

Annotation, guide-design handoff, and standalone Multi-gene alignment use the
atomic directory `intermediate/.annotation-workspace.lock.d` to prevent mixed
generations across host processes and Docker Desktop containers. If a process
is force-killed, the directory intentionally remains and the next run stops
safely. Confirm that no tWhEEL process or container is still using the run,
then remove that lock directory and rerun `annotation` or `sgd_annotation`
before any downstream step. If container cleanup could not be confirmed, the
previous complete consensus set is retained in a hidden
`.multi-gene-alignment-backup-*` directory for recovery.

### Protein-transcript mapping validation fails

Use a headerless `transcript_id<TAB>protein_id` TSV. Its protein IDs must match
the first header tokens in the supplied protein FASTA exactly, and each mapped
transcript must match a non-empty first header token in the transcriptome
FASTA. Do not reverse the two columns. Correct the input files and rerun the
`annotation` or `sgd_annotation` step; mapping inconsistencies are deliberately
not converted into untargetable genes. If a supplied protein intentionally has
no corresponding transcript record, remove it from the protein FASTA instead
of using a partial mapping.

### UniProt API cannot be reached

If a connection or timeout error occurs, ensure that both the host and Docker
can access `https://rest.uniprot.org`. Connect to the internet if needed, then
rerun `go_analysis` with the same configuration file and the complete existing
`tWhEEL_on_...` run name, as shown in the command above. The failed run does not
publish a partially annotated table.

### GO annotations are empty

Check that the best-hit subject identifiers resolve as UniProtKB accessions and
that the UniProt REST API is reachable. A generic custom FASTA identifier is
accepted by annotation but is not automatically convertible to a UniProt
record. A valid UniProtKB accession can also legitimately have no GO
annotation.

### Reports open without styling or images

Open the HTML from its run directory and retain the neighboring `static/` and
`web_data/` directories. The classic payload scripts support direct `file://`
opening; no Flask server is started by the current pipeline.

## Citation

The repository does not currently provide a finalized publication citation.
Until one is available, cite the repository URL together with the Git tag or
commit used and the access date.

## License and support

tWhEEL is distributed under the MIT License; see [LICENSE](LICENSE).

Report reproducible problems through the repository's GitHub Issues page and
include the command, configuration (with sensitive paths removed),
`twheel.log`, and the failing step.

## Related resources

- [TransDecoder](https://github.com/TransDecoder/TransDecoder)
- [BUSCO](https://busco.ezlab.org/)
- [FASTA programs](https://fasta.bioch.virginia.edu/)
- [CrisFlash](https://github.com/tsuryuchka/CRISPRFlash)
- [Gene Ontology](https://geneontology.org/)
- [UniProt](https://www.uniprot.org/)
- [WikiPathways](https://www.wikipathways.org/)
