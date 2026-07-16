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

The current annotation implementations invoke `sudo docker build` for the
local multi-gene alignment image. GO analysis also invokes `sudo docker build`
when one of its three configured helper images is missing. The user running
tWhEEL must therefore be able to execute those commands in the local Docker
setup.

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
docker build -t kazukinakamae/add-go-columns:v1.1 \
  src/docker/add_go_columns/
docker build -t kazukinakamae/compute_go_kfi:v1.1 \
  src/docker/compute_go_kfi/
docker build -t kazukinakamae/count_go_parent_terms:v1.1 \
  src/docker/count_go_parent_terms/
```

Multi-gene alignment uses the local tag below. The current annotation code
rebuilds it during each applicable annotation run even if it already exists.

```bash
docker build -t run_alignment_and_mask/v1.0 \
  src/docker/run_alignment_and_mask/
```

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
reuse the correct directory and its prerequisite intermediate files.

The implemented step names and their main contracts are:

| Step | Main prerequisite | Main output |
| --- | --- | --- |
| `gene_prediction` | Transcriptome; optional existing proteins | `query.pep`, BUSCO output |
| `annotation` | `query.pep`, transcriptome, GO protein DB | `annotation_classified.txt` |
| `sgd_annotation` | `query.pep`, transcriptome, SGD DB | `annotation_classified.txt` |
| `guide_design` | Classified Single/Multi sequences | `annotation_with_guides_final.tsv` |
| `go_analysis` | Final annotation/guide table | GO-KFI tables and GO lookup files |
| `pathway_analysis` | Final annotation/guide table and GPML files | PW-KFI, mappings, and audit tables |
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
- `report.generate_html` and `report.generate_summary` control standard report
  outputs.
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
│   ├── annotation_with_guides_final.tsv
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
    └── report_sources/
        ├── manifest.json
        └── snapshots/
```

Most analysis results remain native TSV, CSV, FASTA, BED, text, or GPML files;
they are not copied into a result database. `web_data/twheel_data.db` contains
only report metadata and cached GO term labels/definitions. The pathway report
does not create a separate SQLite result database.

Report generators atomically publish versioned, checksum-verified snapshots of
the report inputs available to them. The standard report snapshot may contain
only the files present in that run; the pathway report requires a complete
canonical KFI-table/mapping pair. The generator classes can fall back to the
published snapshot when their configured `intermediate` directory does not
exist. The public CLI currently recreates `intermediate` during initialization,
so this fallback is an internal generator contract, not a CLI backup/restore
procedure.

The HTML files are static reports. Keep the adjacent `static/` directory with
the HTML file when moving a report.

The Pathway HTML keeps its interactive data inline so it also works when opened
directly from the filesystem. Its embedded payload contains summary statistics,
the displayed pathway rows, pathway guide data, and generated SVG data once
each. Raw `matched_genes` and pathway mapping rows are not copied into the HTML.
The full machine-readable `web_data/pathway_summary.json` remains a separate
output artifact and is not fetched by the report at runtime.

Report inputs such as annotation text, FASTA identifiers, configuration paths,
GO labels, and GPML labels are treated as untrusted display text. Embedded JSON
uses script-safe Unicode escaping, dynamic tables use text nodes or explicit
HTML escaping, and generated pathway SVG is validated before insertion. These
protections preserve the original strings in embedded data and display them as
text while preventing them from becoming executable HTML or JavaScript.

## Result semantics

### Single and Multi classifications

- **Single**: exactly one query sequence has a given best-hit `subject_id`.
- **Multi**: two or more query sequences share that best-hit `subject_id`.

These are similarity-hit groups, not an independent proof of biological gene or
transcript redundancy. Each group receives a generated `unique_id`, which is
the downstream gene-level counting unit. The final annotation TSV retains its
contig-level rows as provenance, but report gene counts, report guide counts,
GO-KFI, and PW-KFI collapse those rows by `unique_id` and count each distinct
guide entry once within the resulting gene group.

### Site-specific guides

CrisFlash serializes the guide sequence and mismatch counts in the fourth BED
field as `<sequence>/<MM0>/<MM1-2>`. The current filter retains a guide only when
`MM0 == 1` and `MM1-2 == 0`. Any additional score is passed through from
CrisFlash; tWhEEL does not enforce a documented 0-1 range.

### GO-KFI

For each GO term, annotation rows are first collapsed by `unique_id`, guide
targets are de-duplicated, and each gene is counted once:

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

Confirm that the daemon is running and that both `docker` and the current
implementation's `sudo docker build` calls work in your environment. Linux
group-management commands such as `usermod` do not apply to Docker Desktop on
macOS.

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

Open the HTML from its run directory and retain the neighboring `static/`
directory. No Flask server is started by the current pipeline.

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
