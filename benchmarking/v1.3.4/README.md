# tWhEEL v1.3.4 publication benchmark rerun

This directory prepares a reproducible restart at `guide_design`, the first
pipeline step affected by v1.3.4. The runner preserves every v1.3.3 output and
creates a separate v1.3.4 run directory for each analysis.

## Publication scope

`publication_runs.tsv` includes the eight TSA GO-KFI analyses used in Table 1,
the CEN.PK113-7D GO-KFI and PW-KFI analyses, and the two GRCz12tu RefSeq
comparisons used in the manuscript: NP-only and NP+XP. The other RefSeq runs
remain pending and are deliberately excluded.

Supplementary Table S2 is not part of this rerun plan. The NP+XP run remains in
scope only because its aggregate GO-KFI result is discussed in the manuscript.

## Restart boundary

The runner clones these validated upstream artifacts from each v1.3.3 run:

- gene-prediction and annotation result summaries;
- the prepared transcriptome and protein FASTA files;
- the classified annotation and Single/Multi group FASTA files; and
- the aligned Multi consensus directory.

macOS copy-on-write cloning (`cp -cR`) keeps the new workspace independent
without duplicating the physical storage immediately. Source and cloned
artifacts are checksum-verified before the temporary seed directory is renamed
to its final run name. The source output is never modified.
An existing seed is reused only after its manifest metadata and every cloned
artifact still match the current source, configuration, and guide-design code.

From that seed, the runner executes:

1. `guide_design` for every run;
2. `go_analysis` and `report_generation` for GO-KFI runs; or
3. `pathway_analysis` and `pathway_report_generation` for the PW-KFI run.

## Commands

Run these commands from the repository root:

```bash
bash benchmarking/v1.3.4/rerun_from_guide_design.sh --list
bash benchmarking/v1.3.4/rerun_from_guide_design.sh --preflight-only
bash benchmarking/v1.3.4/rerun_from_guide_design.sh --prepare
```

`--preflight-only` is the default and performs no writes. `--prepare` creates
the independently named v1.3.4 workspaces but does not start analysis.

Start with the affected Arabidopsis TSA run:

```bash
bash benchmarking/v1.3.4/rerun_from_guide_design.sh \
  --execute --run tsa_gkif01
```

After reviewing that result, execute or resume all manifest entries:

```bash
bash benchmarking/v1.3.4/rerun_from_guide_design.sh --execute
bash benchmarking/v1.3.4/verify_publication_reruns.sh
```

Set `TWHEEL_V134_PYTHON` only when a different validated Python environment is
required. Use `--skip-docker-check` only for file-level preparation on a host
where Docker is intentionally unavailable; do not use it for an analysis run.

## Release provenance

Each seeded output records the source run, configuration checksum, current Git
commit, and the SHA-256 checksum of `guide_design.py`. The canonical runtime
version is now 1.3.4. A Git commit, `v1.3.4` tag, and immutable software archive
should be created only after code review and before the formal publication
rerun; this preparation does not create them automatically.
