from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
VALIDATOR = REPOSITORY / "bin" / "check_samplesheet.py"


def run_validator(root: Path, header: str, rows: list[str], group_column: str = ""):
    for sample in ("sample_a", "sample_b", "sample_c"):
        for mate in (1, 2):
            (root / f"{sample}_R{mate}.fastq.gz").write_bytes(b"synthetic\n")
    source = root / "samples.csv"
    output = root / "validated.csv"
    metadata = root / "sample_metadata.tsv"
    source.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    arguments = [
        sys.executable,
        str(VALIDATOR),
        "--input",
        str(source),
        "--output",
        str(output),
        "--metadata-output",
        str(metadata),
    ]
    if group_column:
        arguments.extend(["--group-column", group_column])
    result = subprocess.run(arguments, text=True, capture_output=True, check=False)
    return result, output, metadata


class CheckSamplesheetTests(unittest.TestCase):
    def test_classic_three_column_samplesheet_remains_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result, output, metadata = run_validator(
                root,
                "sample,fastq_1,fastq_2",
                ["sample_a,sample_a_R1.fastq.gz,sample_a_R2.fastq.gz"],
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["group"], "")
            with metadata.open(encoding="utf-8") as handle:
                metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(metadata_rows[0]["sample_id"], "sample_a")
            self.assertEqual(metadata_rows[0]["group"], "")

    def test_configured_group_is_required_and_preserved_for_any_number_of_categories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result, output, metadata = run_validator(
                root,
                "sample,fastq_1,fastq_2,condition",
                [
                    "sample_a,sample_a_R1.fastq.gz,sample_a_R2.fastq.gz,Treatment_A",
                    "sample_b,sample_b_R1.fastq.gz,sample_b_R2.fastq.gz,Treatment_B",
                    "sample_c,sample_c_R1.fastq.gz,sample_c_R2.fastq.gz,Treatment_C",
                ],
                "condition",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(
                [row["group"] for row in rows],
                ["Treatment_A", "Treatment_B", "Treatment_C"],
            )
            with metadata.open(encoding="utf-8") as handle:
                metadata_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [row["group"] for row in metadata_rows],
                ["Treatment_A", "Treatment_B", "Treatment_C"],
            )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            missing, output, metadata = run_validator(
                root,
                "sample,fastq_1,fastq_2,condition",
                ["sample_a,sample_a_R1.fastq.gz,sample_a_R2.fastq.gz,"],
                "condition",
            )
            self.assertEqual(missing.returncode, 1)
            self.assertIn("group column 'condition' is empty", missing.stderr)
            self.assertFalse(output.exists())
            self.assertFalse(metadata.exists())

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            multiple, output, metadata = run_validator(
                root,
                "sample,fastq_1,fastq_2,condition",
                [
                    "sample_a,sample_a_R1.fastq.gz,sample_a_R2.fastq.gz,"
                    "Treatment_A,unexpected_second_value"
                ],
                "condition",
            )
            self.assertEqual(multiple.returncode, 1)
            self.assertIn("unexpected extra fields", multiple.stderr)
            self.assertFalse(output.exists())
            self.assertFalse(metadata.exists())


if __name__ == "__main__":
    unittest.main()
