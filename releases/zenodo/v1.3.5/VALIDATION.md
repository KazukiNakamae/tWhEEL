# Archive preparation validation

Verified on 2026-10-26. Preparation checks are followed by publication checks
below; the original source ZIP remains unchanged.

| Check | Observed result |
| --- | --- |
| Remote v1.3.5 tag | Remote tag object and peeled commit match local values recorded in `README.md`. |
| Source version | Tagged `src/twheel_metadata.py` declares `1.3.5`. |
| Source ZIP | 96 files, 1,038,561 bytes; every file's Git blob hash matches the pinned tree. |
| ZIP SHA-256 | `2c47af2d1b625f5c5ae84cee5d9659d39d29b40a3d60b0294cc228e74c735c90` |
| Independent rebuild | Source ZIP rebuilt to a second directory is byte-identical. |
| Existing destination | Build refuses overwrite; original ZIP remains unchanged. |
| CFF format | Passes the official CFF 1.2.0 JSON Schema with format checking. |
| Metadata agreement | Title, version, both creators in the supplied order, affiliations, ORCIDs, and license agree between CFF and Zenodo JSON. |
| Creator provenance | Kazuki Nakamae followed by Hidemasa Bono, with affiliations and ORCIDs supplied by the owner on 2026-10-26. Both ORCID checksums are valid; identity ownership was not independently verified. |
| DOI fields | The version DOI `10.5281/zenodo.23660613` is registered and resolves to the published record. |
| GO-KFI report regressions | 8 passed in 1.37 s, using code extracted from the source ZIP. |
| Git whitespace check | `git diff --check` passed. |

The report tests ran under the existing Python 3.12 environment at
`.handoff_runtime/twheel_env-osx-64`, with Node.js available for the generated
JavaScript assertions. The test command, run from a temporary extraction, was:

```bash
/path/to/python3.12 -B -m pytest -q -p no:cacheprovider \
  tests/test_report_go_kfi_chart.py
```

The official schema was retrieved from
<https://raw.githubusercontent.com/citation-file-format/citation-file-format/1.2.0/schema.json>.
Validation used PyYAML and jsonschema 4.26.0 with
`Draft7Validator(schema, format_checker=FormatChecker())`. The jsonschema
dependencies were installed only in a temporary directory; the existing
analysis environment was not changed.

The eight regression tests execute the real generated chart methods in a
minimal JavaScript DOM harness. They are not a real-browser inspection,
Docker end-to-end analysis, rerun of the publication benchmark, or biological
validation. No such additional validation was performed for archive preparation.

## Publication verification

- Published record: <https://zenodo.org/records/23660613>.
- `PUBLISHED_RECORD.json` was exported from the published record page;
  it reports `status: published` and `is_published: true`.
- Version DOI: `10.5281/zenodo.23660613`; concept DOI:
  `10.5281/zenodo.23660612`. The version DOI was followed in the browser and
  resolved to the published record without a preview parameter.
- Publication date: 2026-10-26. Resource type Software; access Open; license MIT.
- Both creators, their order, their supplied ORCIDs, and the two distinct
  affiliations for Hidemasa Bono were verified in the published metadata.
- The record contains one source ZIP, 1,038,561 bytes. Its published MD5 is
  `aa82a0112de31104a5a6adf08446d1fd`, matching the local archive.
- The ZIP was downloaded from the published record and its SHA-256 matches
  `2c47af2d1b625f5c5ae84cee5d9659d39d29b40a3d60b0294cc228e74c735c90`.
- Software Heritage archival remains unconfirmed; the observed external
  resources panel listed OpenAIRE but no Software Heritage archive.

## GitHub release verification

The public GitHub Release was published on 2026-10-26 at
<https://github.com/KazukiNakamae/tWhEEL/releases/tag/v1.3.5>. It is neither a
draft nor a prerelease. GitHub reports the attached ZIP SHA-256 as
`2c47af2d1b625f5c5ae84cee5d9659d39d29b40a3d60b0294cc228e74c735c90`,
matching Zenodo. A separate SHA-256 checksum asset accompanies the ZIP.
