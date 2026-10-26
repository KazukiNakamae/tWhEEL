#!/usr/bin/env bash
# Build an exact tagged-source ZIP without using or modifying the working tree.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(git -C "$script_dir" rev-parse --show-toplevel)"
python3 - "$repo_dir" "$script_dir" "${1:-$repo_dir/dist/zenodo/v1.3.5}" <<'PY'
import ast
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

repo, metadata_dir, destination = map(Path, sys.argv[1:])
destination = destination.resolve()
tag = "v1.3.5"
commit = "77d7a7c9aa14ef1d5fc8e709293ab86999226d8a"
prefix = "tWhEEL-1.3.5/"

def git(*args):
    return subprocess.check_output(["git", "-C", str(repo), *args])

def sha256(data):
    return hashlib.sha256(data).hexdigest()

if git("rev-parse", f"refs/tags/{tag}^{{commit}}").decode().strip() != commit:
    raise SystemExit("The v1.3.5 tag does not match the pinned release commit.")
tree = git("rev-parse", f"{commit}^{{tree}}").decode().strip()
module = ast.parse(git("show", f"{commit}:src/twheel_metadata.py"))
version = next(ast.literal_eval(node.value) for node in module.body
               if isinstance(node, ast.Assign)
               and any(isinstance(t, ast.Name) and t.id == "__version__"
                       for t in node.targets))
if version != "1.3.5":
    raise SystemExit("Tagged runtime version is not 1.3.5.")
metadata = json.loads((metadata_dir / ".zenodo.json").read_text())
if metadata["version"] != version or metadata["upload_type"] != "software":
    raise SystemExit("Zenodo metadata does not describe this software version.")
if destination.exists():
    raise SystemExit(f"Refusing to overwrite existing destination: {destination}")
if destination == metadata_dir.resolve() or destination in metadata_dir.resolve().parents:
    raise SystemExit("Destination must not contain the preparation source files.")
destination.parent.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix=".twheel-zenodo-", dir=destination.parent) as tmp:
    staging = Path(tmp)
    archive = staging / "tWhEEL-1.3.5-source.zip"
    subprocess.run(["git", "-C", str(repo), "archive", "--format=zip",
                    f"--prefix={prefix}", f"--output={archive}", commit], check=True)
    entries = {}
    for row in git("ls-tree", "-rz", commit).split(b"\0"):
        if not row:
            continue
        attributes, raw_path = row.split(b"\t", 1)
        mode, kind, oid = attributes.decode().split()
        if kind != "blob" or mode not in {"100644", "100755"}:
            raise SystemExit("Unsupported Git tree entry; inspect before archiving.")
        entries[raw_path.decode()] = (mode, oid)
    manifest = []
    with zipfile.ZipFile(archive) as zipped:
        files = [item for item in zipped.infolist() if not item.is_dir()]
        if len(files) != len(entries) or {item.filename for item in files} != {
                prefix + name for name in entries}:
            raise SystemExit("Archive members differ from the pinned Git tree.")
        for item in files:
            path = item.filename[len(prefix):]
            data = zipped.read(item)
            mode, expected_blob = entries[path]
            actual_blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
            if actual_blob != expected_blob:
                raise SystemExit(f"Archive content differs from Git: {path}")
            manifest.append({"path": path, "mode": mode, "git_blob": expected_blob,
                             "size_bytes": len(data), "sha256": sha256(data)})
    for name in ("CITATION.cff", ".zenodo.json", "RELEASE_NOTES.md", "README.md", "VALIDATION.md", "PUBLISHED_RECORD.json"):
        shutil.copyfile(metadata_dir / name, staging / name)
    provenance = {
        "software": "tWhEEL", "version": version, "tag": tag,
        "commit": commit, "tree": tree,
        "tag_object": git("rev-parse", f"refs/tags/{tag}").decode().strip(),
        "source_date": "2026-09-27",
        "repository": "https://github.com/KazukiNakamae/tWhEEL",
        "archive": archive.name, "archive_sha256": sha256(archive.read_bytes()),
        "file_count": len(manifest),
        "content_verification": "Every file matches the pinned Git blob; no exclusions or additions.",
        "metadata_scope": "Separate preparation files created after the source tag; not in source ZIP.",
        "publication_status": "Published on Zenodo on 2026-10-26; downloaded source ZIP verified.",
        "zenodo_record": "https://zenodo.org/records/23660613",
        "version_doi": "10.5281/zenodo.23660613",
        "concept_doi": "10.5281/zenodo.23660612",
        "creator_status": "Creator order, affiliations, and ORCIDs supplied by the owner on 2026-10-26.",
        "creators": metadata["creators"],
        "files": sorted(manifest, key=lambda item: item["path"]),
    }
    (staging / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    checksums = [f"{sha256(path.read_bytes())}  {path.name}\n"
                 for path in sorted(staging.iterdir())]
    (staging / "SHA256SUMS").write_text("".join(checksums))
    staging.rename(destination)
print(f"Prepared {destination}")
print(f"Verified {len(manifest)} files at {commit}")
print(f"Source ZIP SHA-256: {provenance['archive_sha256']}")
PY
