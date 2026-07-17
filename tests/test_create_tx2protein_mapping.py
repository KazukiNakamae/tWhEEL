import contextlib
import gzip
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "docker"
    / "create_tx2protein_mapping"
    / "create_tx2protein_mapping.py"
)
SPEC = importlib.util.spec_from_file_location(
    "create_tx2protein_mapping", SCRIPT_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

MAPPING_SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "modules"
    / "transcript_mapping.py"
)
MAPPING_SPEC = importlib.util.spec_from_file_location(
    "transcript_mapping", MAPPING_SCRIPT_PATH
)
TRANSCRIPT_MAPPING = importlib.util.module_from_spec(MAPPING_SPEC)
MAPPING_SPEC.loader.exec_module(TRANSCRIPT_MAPPING)


REFSEQ_GFF = (
    "##gff-version 3\n"
    "#!processor NCBI annotwriter\n"
    "#!annotation-source NCBI RefSeq test\n"
    "chr1\tBestRefSeq\tCDS\t1\t30\t.\t+\t0\t"
    "ID=cds-NP_000001.1;Parent=rna-NM_000001.1;"
    "protein_id=NP_000001.1\n"
    "chr1\tBestRefSeq\tCDS\t31\t60\t.\t+\t0\t"
    "ID=cds-NP_000001.1;Parent=rna-NM_000001.1;"
    "protein_id=NP_000001.1\n"
    "chr1\tGnomon\tCDS\t100\t150\t.\t+\t0\t"
    "ID=cds-XP_000002.2;Parent=rna-XM_000002.3;"
    "protein_id=XP_000002.2\n"
    "chr1\tRefSeq\tCDS\t200\t230\t.\t+\t0\t"
    "ID=cds-rearrangement;Parent=rna-VDJ_segment\n"
)
REFSEQ_TRANSCRIPTS = """>NM_000001.1 curated transcript
ATGAAATAG
>XM_000002.3 predicted transcript
ATGCCCTAG
>XR_000003.1 non-coding transcript
AAAA
"""
REFSEQ_PROTEINS = """>NP_000001.1 curated protein
MK
>XP_000002.2 predicted protein
MP
"""


class CreateTx2ProteinMappingTests(unittest.TestCase):
    def _write(self, path, text, compressed=False):
        if compressed:
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(text)
        else:
            path.write_text(text, encoding="utf-8")

    def _refseq_inputs(self, root, compressed=False):
        suffix = ".gz" if compressed else ""
        gff = root / f"annotation.gff3{suffix}"
        transcripts = root / f"rna.fna{suffix}"
        proteins = root / f"protein.faa{suffix}"
        self._write(gff, REFSEQ_GFF, compressed)
        self._write(transcripts, REFSEQ_TRANSCRIPTS, compressed)
        self._write(proteins, REFSEQ_PROTEINS, compressed)
        return gff, transcripts, proteins

    def _run(self, arguments):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = MODULE.main([str(value) for value in arguments])
        return status, stdout.getvalue(), stderr.getvalue()

    def test_refseq_include_xp_creates_exact_mapping_and_audit(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff, transcripts, proteins = self._refseq_inputs(root)
            output = root / "tx2protein.tsv"

            status, stdout, stderr = self._run(
                [
                    "--source", "auto",
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(0, status, stderr)
            self.assertEqual(
                "NM_000001.1\tNP_000001.1\nXM_000002.3\tXP_000002.2\n",
                output.read_text(encoding="utf-8"),
            )
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual("ok", audit["status"])
            self.assertEqual("refseq", audit["source_detected"])
            self.assertEqual(2, audit["counts"]["output_mapping_rows"])
            self.assertEqual(0, audit["counts"]["excluded_xp_proteins"])
            self.assertEqual(
                1, audit["counts"]["cds_rows_without_protein_id_attribute"]
            )
            self.assertIn("Mappings: 2", stdout)
            mapping = TRANSCRIPT_MAPPING.load_mapping_file(output)
            validation = TRANSCRIPT_MAPPING.validate_mapping_against_fastas(
                mapping, proteins, transcripts
            )
            self.assertEqual(2, validation["protein_count"])

    def test_refseq_exclude_xp_creates_matching_filtered_proteins(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff, transcripts, proteins = self._refseq_inputs(root)
            output = root / "refseq_tx2protein.tsv"

            status, _, stderr = self._run(
                [
                    "--source", "refseq",
                    "--refseq-xp", "exclude",
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(0, status, stderr)
            self.assertEqual(
                "NM_000001.1\tNP_000001.1\n",
                output.read_text(encoding="utf-8"),
            )
            filtered = root / "refseq_tx2protein.matched.faa"
            filtered_text = filtered.read_text(encoding="utf-8")
            self.assertIn(">NP_000001.1", filtered_text)
            self.assertNotIn(">XP_000002.2", filtered_text)
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual(1, audit["counts"]["excluded_xp_proteins"])
            self.assertEqual(str(filtered), audit["outputs"]["filtered_proteins"])
            mapping = TRANSCRIPT_MAPPING.load_mapping_file(output)
            validation = TRANSCRIPT_MAPPING.validate_mapping_against_fastas(
                mapping, filtered, transcripts
            )
            self.assertEqual(1, validation["protein_count"])

    def test_ensembl_reconciles_unique_version_suffixes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff = root / "annotation.gff3"
            transcripts = root / "cdna.fa"
            proteins = root / "pep.fa"
            output = root / "tx2protein.tsv"
            self._write(
                gff,
                "##gff-version 3\n"
                "chr1\tensembl\tCDS\t1\t90\t.\t+\t0\t"
                "ID=CDS:ENSP000001;Parent=transcript:ENST000001;"
                "protein_id=ENSP000001\n"
                "chr1\tensembl\tCDS\t100\t190\t.\t+\t0\t"
                "ID=CDS:ENSP000001;Parent=transcript:ENST000001;"
                "protein_id=ENSP000001\n",
            )
            self._write(transcripts, ">ENST000001.2 cdna\nATGAAATAG\n")
            self._write(
                proteins,
                ">ENSP000001.5 pep\nMK\n",
            )

            status, _, stderr = self._run(
                [
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(0, status, stderr)
            self.assertEqual(
                "ENST000001.2\tENSP000001.5\n",
                output.read_text(encoding="utf-8"),
            )
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual("ensembl", audit["source_detected"])
            self.assertGreaterEqual(
                audit["counts"]["version_normalized_pairs"], 1
            )

    def test_ambiguous_protein_mapping_fails_without_mapping_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff = root / "annotation.gff3"
            transcripts = root / "rna.fa"
            proteins = root / "protein.fa"
            output = root / "tx2protein.tsv"
            self._write(
                gff,
                "##gff-version 3\n#!processor NCBI annotwriter\n"
                "chr1\tRefSeq\tCDS\t1\t9\t.\t+\t0\t"
                "Parent=rna-NM_1.1;protein_id=NP_1.1\n"
                "chr1\tRefSeq\tCDS\t20\t30\t.\t+\t0\t"
                "Parent=rna-NM_2.1;protein_id=NP_1.1\n",
            )
            self._write(transcripts, ">NM_1.1\nATG\n>NM_2.1\nATG\n")
            self._write(proteins, ">NP_1.1\nM\n")

            status, _, stderr = self._run(
                [
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(2, status)
            self.assertFalse(output.exists())
            self.assertIn("mapping validation failed", stderr)
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertIn(
                "ambiguous_protein_mappings",
                {problem["code"] for problem in audit["problems"]},
            )

    def test_missing_protein_mapping_fails(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff, transcripts, proteins = self._refseq_inputs(root)
            with proteins.open("a", encoding="utf-8") as handle:
                handle.write(">NP_999999.1 absent from GFF\nMZ\n")
            output = root / "tx2protein.tsv"

            status, _, _ = self._run(
                [
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(2, status)
            self.assertFalse(output.exists())
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            problem = next(
                value for value in audit["problems"]
                if value["code"] == "unmapped_proteins"
            )
            self.assertEqual(["NP_999999.1"], problem["examples"])

    def test_unmapped_protein_can_be_explicitly_excluded(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff, transcripts, proteins = self._refseq_inputs(root)
            with gff.open("a", encoding="utf-8") as handle:
                handle.write(
                    "chrM\tRefSeq\tCDS\t1\t30\t.\t+\t0\t"
                    "ID=cds-NP_059331.1;Parent=rna-ND1;"
                    "protein_id=NP_059331.1\n"
                )
            with proteins.open("a", encoding="utf-8") as handle:
                handle.write(">NP_059331.1 mitochondrial protein\nMZ\n")
            output = root / "tx2protein.tsv"

            status, _, stderr = self._run(
                [
                    "--unmapped-proteins", "exclude",
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(0, status, stderr)
            matched = root / "tx2protein.matched.faa"
            matched_text = matched.read_text(encoding="utf-8")
            self.assertNotIn(">NP_059331.1", matched_text)
            self.assertIn(">NP_000001.1", matched_text)
            self.assertIn(">XP_000002.2", matched_text)
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual(1, audit["counts"]["excluded_unmapped_proteins"])
            self.assertEqual(
                ["NP_059331.1"],
                audit["examples"]["excluded_unmapped_proteins"],
            )

    def test_ambiguous_versionless_transcript_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff = root / "annotation.gff3"
            transcripts = root / "cdna.fa"
            proteins = root / "pep.fa"
            output = root / "tx2protein.tsv"
            self._write(
                gff,
                "##gff-version 3\n"
                "chr1\tensembl\tCDS\t1\t90\t.\t+\t0\t"
                "ID=CDS:ENSP1;Parent=transcript:ENST1;protein_id=ENSP1\n",
            )
            self._write(
                transcripts,
                ">ENST1.1 first version\nATG\n"
                ">ENST1.2 second version\nATG\n",
            )
            self._write(proteins, ">ENSP1\nM\n")

            status, _, _ = self._run(
                [
                    "--unmapped-proteins", "exclude",
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(2, status)
            self.assertFalse(output.exists())
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertIn(
                "ambiguous_transcript_identifiers",
                {problem["code"] for problem in audit["problems"]},
            )

    def test_conflicting_explicit_versions_are_not_normalized(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff = root / "annotation.gff3"
            transcripts = root / "cdna.fa"
            proteins = root / "pep.fa"
            output = root / "tx2protein.tsv"
            self._write(
                gff,
                "##gff-version 3\n"
                "chr1\tensembl\tCDS\t1\t90\t.\t+\t0\t"
                "ID=CDS:ENSP1.1;Parent=transcript:ENST1.1;"
                "protein_id=ENSP1.1\n",
            )
            self._write(transcripts, ">ENST1.2\nATG\n")
            self._write(proteins, ">ENSP1.2\nM\n")

            status, _, _ = self._run(
                [
                    "--unmapped-proteins", "exclude",
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(2, status)
            self.assertFalse(output.exists())
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertIn(
                "conflicting_protein_versions",
                {problem["code"] for problem in audit["problems"]},
            )

    def test_gzip_inputs_are_supported(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff, transcripts, proteins = self._refseq_inputs(
                root, compressed=True
            )
            output = root / "tx2protein.tsv"

            status, _, stderr = self._run(
                [
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(0, status, stderr)
            self.assertEqual(2, len(output.read_text(encoding="utf-8").splitlines()))

    def test_explicit_source_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            gff, transcripts, proteins = self._refseq_inputs(root)
            output = root / "tx2protein.tsv"

            status, _, _ = self._run(
                [
                    "--source", "ensembl",
                    "--annotation", gff,
                    "--transcripts", transcripts,
                    "--proteins", proteins,
                    "--output", output,
                ]
            )

            self.assertEqual(2, status)
            audit = json.loads(
                Path(f"{output}.audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual("source_mismatch", audit["problems"][0]["code"])


if __name__ == "__main__":
    unittest.main()
