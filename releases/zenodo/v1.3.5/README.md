# Zenodo preparation for tWhEEL v1.3.5

## Status and identity

The v1.3.5 source archive was published on Zenodo on 2026-10-26:

- Record: <https://zenodo.org/records/23660613>
- Version DOI: <https://doi.org/10.5281/zenodo.23660613>
- All-versions DOI: <https://doi.org/10.5281/zenodo.23660612>
- Published metadata export: `PUBLISHED_RECORD.json`.

The version DOI resolves to the published record. Its source ZIP was downloaded
and verified byte-for-byte by SHA-256 against the prepared archive. Only the
source ZIP was uploaded; the preparation files remain local supporting records.
Software Heritage archival has not yet been confirmed.

The source is pinned to tag `v1.3.5`, commit
`77d7a7c9aa14ef1d5fc8e709293ab86999226d8a`, tree
`d79595e0665aedaee187e25333dde274363b4006`. The tag/source date is
2026-09-27; the Zenodo publication date is 2026-10-26.

On 2026-10-26, `git ls-remote` verified the remote tag object
`83f5d5ff29bec8f56c2ad6d96d8e6e40f155e7c3` and its peeled commit against
the local tag. The GitHub connector also returned this commit. During initial preparation,
GitHub Release availability was unverified and remote main was older than the
correction. The repository is now public, main includes v1.3.5, and the
[GitHub Release](https://github.com/KazukiNakamae/tWhEEL/releases/tag/v1.3.5)
was published on 2026-10-26. Its attached source ZIP has the same SHA-256 as the
Zenodo archive. Continue to use the pinned tag for source preservation.

The software creators, in the order supplied by the owner on 2026-10-26, are:

1. **Kazuki Nakamae** — Genome Editing Innovation Center, Hiroshima University,
   Higashi-Hiroshima, 739-0046, Japan.
   ORCID: <https://orcid.org/0000-0002-4469-664X>.
2. **Hidemasa Bono** — Genome Editing Innovation Center, Hiroshima University,
   Higashi-Hiroshima, 739-0046, Japan; Graduate School of Integrated Sciences for
   Life, Hiroshima University, Higashi-Hiroshima, 739-0046, Japan.
   ORCID: <https://orcid.org/0000-0003-4413-0651>.

Both metadata files preserve this order and the supplied affiliations/ORCIDs.
The two affiliations for Hidemasa Bono share one semicolon-separated affiliation
field in CFF and the legacy Zenodo JSON format. When entering the Zenodo form,
add them as two affiliations if the form supports separate entries.
No paper DOI, grant, or Zenodo DOI has been invented.

## Build and verify

From the repository root:

```bash
bash releases/zenodo/v1.3.5/build_archive.sh
cd dist/zenodo/v1.3.5
shasum -a 256 -c SHA256SUMS
```

Only Git, Bash, and Python 3's standard library are needed. The build refuses
to overwrite an existing destination. Pass a new output directory as its first
argument to rebuild. The ZIP comes from the pinned Git commit; local edits and
untracked data cannot enter it. The builder verifies the tagged runtime version
and every ZIP member against the Git tree and emits per-file SHA-256 checksums
in `PROVENANCE.json`.

The source ZIP preserves all tracked files exactly. It therefore retains old
tracked `__pycache__` files and `outputs/test_report_fix/static` assets. This is
an exact source snapshot, not a cleaned or rebuilt distribution. No container
images, reference databases, or untracked analysis results are bundled.

`CITATION.cff`, `.zenodo.json`, these instructions, and release notes describe
the archive but are separate from the historical source ZIP. These files are
kept in this version-specific preparation directory and were added after the
v1.3.5 tag. Root-level copies of `CITATION.cff` and `.zenodo.json` provide
GitHub citation and integration metadata. Do not move the old tag to include them.

## Manual registration reference (completed for this version)

The following steps document the workflow. Do not create another v1.3.5 deposit;
use the published record above. For future releases, use its **New version**
action to retain the existing concept DOI.

1. Sign in at <https://zenodo.org> and inspect your records and GitHub
   integration for an existing tWhEEL deposit. If a software record already
   exists, use its **New version** action; avoid creating an unrelated record
   and a second concept DOI. An exact v1.3.5 record should be reused.
2. Create a new upload only when no record exists. Upload
   `tWhEEL-1.3.5-source.zip`. For Software Heritage ingestion, Zenodo requires
   a single compressed source file in the record. `PROVENANCE.json`,
   `SHA256SUMS`, `CITATION.cff`, and `RELEASE_NOTES.md` may be attached as
   uncompressed supporting files. Do not upload a second ZIP or tarball.
3. Set resource type to **Software**, version to **1.3.5**, language to
   **English**, access to **Open**, and license to **MIT**. Copy the title,
   description, creator names/order, keywords, and related identifiers from
   `.zenodo.json`. Uploading that JSON as a normal file does not populate the
   form; the automatic `.zenodo.json` convention applies to GitHub integration.
4. Enter both creators in the supplied order, with the affiliations and ORCIDs
   above. Add grants, communities, or a related publication only when their
   exact values and applicability are established. Keep the software release
   date (2026-09-27) distinct from the
   actual Zenodo publication date. Do not invent a publication DOI.
5. Review the rendered metadata and uploaded ZIP checksum. Reserve a DOI if
   needed for review. A reserved DOI is not evidence of a published archive.
   Publish only after reviewing the complete deposit.
6. Record the published record URL, version-specific DOI, concept DOI, and
   publication date. Verify downloads against `SHA256SUMS`. Check the record's
   **External resources > Archived in** status separately; Software Heritage
   archival can follow later.
7. Add the version DOI to this release's CFF `doi` field, cite that DOI for
   v1.3.5, and link the published record from repository documentation. Use the
   concept DOI only for citations intended to cover all versions. Keep a copy
   of the deposited metadata and upload receipt with the DOI.

## GitHub integration alternative

If Zenodo has not yet ingested the release, enable this repository in Zenodo
and check whether the existing tag can be selected for a GitHub Release.
Verify the release target resolves to the pinned commit before publishing.
Do not create another release if the version has already been archived.

Adding metadata to the current branch does not retroactively add it to the
v1.3.5 tag. For this existing tag, use manual deposit metadata or edit the
generated Zenodo draft after ingestion, before publishing. For future versions,
place updated `CITATION.cff` and `.zenodo.json` at the repository root before
tagging. For a new release, update the version and release date and remove the
previous version DOI from both metadata files until the new DOI is known; do not
reuse the v1.3.5 DOI for a different release. If both files exist, Zenodo uses
`.zenodo.json`; keep their creators, title,
license, and version synchronized. GitHub's citation display uses the root CFF.

## Verification limits

Archive integrity checks establish source identity, not successful Docker
execution or experimental knockout performance. The v1.3.5 regression suite
can verify report behavior separately. Publication benchmarks and external
tools must retain their own evidence and version provenance. This preparation
does not freeze external Docker tags or online annotation responses.

## Official references

- [Upload software manually](https://help.zenodo.org/docs/github/archive-software/manual-upload/)
- [Archive a release from GitHub](https://help.zenodo.org/docs/github/archive-software/github-upload/)
- [Zenodo JSON metadata](https://help.zenodo.org/docs/github/describe-software/zenodo-json/)
- [CITATION.cff metadata](https://help.zenodo.org/docs/github/describe-software/citation-file/)
