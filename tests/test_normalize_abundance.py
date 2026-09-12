from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
NORMALIZER = REPOSITORY / "bin" / "normalize_coverm_abundance.py"
OUTPUT_FIELDS = (
    "sample",
    "group",
    "mag_id",
    "relative_abundance_percent",
    "mean_coverage",
    "covered_fraction",
    "genome_length",
)


def coverm_header(*samples: str) -> str:
    fields = ["Genome"]
    for sample in samples:
        fields.extend(
            [
                f"{sample} Relative Abundance (%)",
                f"{sample} Mean",
                f"{sample} Covered Fraction",
                f"{sample} Length",
            ]
        )
    return "\t".join(fields)


def run_normalizer(
    root: Path,
    content: str,
    groups: dict[str, str] | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    source = root / "coverm.tsv"
    metadata = root / "sample_metadata.tsv"
    output = root / "abundance.long.tsv"
    source.write_text(content, encoding="utf-8")
    sample_names: list[str] = []
    for field in content.splitlines()[0].split("\t")[1:]:
        for suffix in (
            " Relative Abundance (%)",
            " Mean",
            " Covered Fraction",
            " Length",
        ):
            if field.endswith(suffix):
                sample = field[: -len(suffix)]
                if sample and sample not in sample_names:
                    sample_names.append(sample)
                break
    if not sample_names:
        sample_names = ["metadata_only_sample"]
    selected_groups = groups or {sample: "" for sample in sample_names}
    metadata.write_text(
        "sample_id\tgroup\n"
        + "".join(f"{sample}\t{selected_groups.get(sample, '')}\n" for sample in sample_names),
        encoding="utf-8",
    )
    result = subprocess.run(
        [
            sys.executable,
            str(NORMALIZER),
            "--input",
            str(source),
            "--sample-metadata",
            str(metadata),
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    return result, output


class NormalizeCovermAbundanceTests(unittest.TestCase):
    def test_two_samples_by_two_mags_excludes_unmapped_without_changing_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            content = (
                coverm_header("sample_b", "sample_a")
                + "\n"
                + "MAG_2\t2.5\t3.25\t0.75\t2048\t4.5\t6.25\t0.8\t2048\n"
                + "unmapped\t70\t0\t0\t1\t60\t0\t0\t1\n"
                + "MAG_1\t7.125\t8.5\t1\t4096\t9.25\t10.75\t0.125\t4096\n"
            )
            result, output = run_normalizer(
                root,
                content,
                {"sample_a": "Control", "sample_b": "PD"},
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames, list(OUTPUT_FIELDS))
                rows = list(reader)

            self.assertEqual(
                rows,
                [
                    {
                        "sample": "sample_a",
                        "group": "Control",
                        "mag_id": "MAG_2",
                        "relative_abundance_percent": "4.5",
                        "mean_coverage": "6.25",
                        "covered_fraction": "0.8",
                        "genome_length": "2048",
                    },
                    {
                        "sample": "sample_b",
                        "group": "PD",
                        "mag_id": "MAG_2",
                        "relative_abundance_percent": "2.5",
                        "mean_coverage": "3.25",
                        "covered_fraction": "0.75",
                        "genome_length": "2048",
                    },
                    {
                        "sample": "sample_a",
                        "group": "Control",
                        "mag_id": "MAG_1",
                        "relative_abundance_percent": "9.25",
                        "mean_coverage": "10.75",
                        "covered_fraction": "0.125",
                        "genome_length": "4096",
                    },
                    {
                        "sample": "sample_b",
                        "group": "PD",
                        "mag_id": "MAG_1",
                        "relative_abundance_percent": "7.125",
                        "mean_coverage": "8.5",
                        "covered_fraction": "1",
                        "genome_length": "4096",
                    },
                ],
            )
            self.assertNotIn("unmapped", {row["mag_id"] for row in rows})

    def test_group_metadata_accepts_more_than_two_categories_without_changing_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            content = (
                coverm_header("sample_a", "sample_b", "sample_c")
                + "\nMAG_1\t1\t2\t0.5\t1000\t3\t4\t0.6\t1000\t5\t6\t0.7\t1000\n"
            )
            result, output = run_normalizer(
                root,
                content,
                {"sample_a": "Treatment_A", "sample_b": "Treatment_B", "sample_c": "Treatment_C"},
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {row["sample"]: row["group"] for row in rows},
                {
                    "sample_a": "Treatment_A",
                    "sample_b": "Treatment_B",
                    "sample_c": "Treatment_C",
                },
            )
            self.assertEqual(
                [row["relative_abundance_percent"] for row in rows],
                ["1", "3", "5"],
            )

    def test_unmapped_filter_is_case_sensitive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            content = (
                coverm_header("sample")
                + "\n"
                + "unmapped\t90\t0\t0\t1\n"
                + "Unmapped\t10\t1.5\t0.5\t1000\n"
            )
            result, output = run_normalizer(root, content)

            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["mag_id"], "Unmapped")

    def test_rejects_no_samples_and_incomplete_sample_metrics(self) -> None:
        cases = (
            (
                "no sample metrics",
                "Genome\nMAG_1\n",
                "contains no sample metric columns",
            ),
            (
                "missing metric",
                "Genome\tsample Relative Abundance (%)\tsample Mean\t"
                "sample Covered Fraction\nMAG_1\t10\t1\t0.5\n",
                "lacks metrics: genome_length",
            ),
        )
        for label, content, expected_error in cases:
            with self.subTest(case=label), tempfile.TemporaryDirectory() as temporary:
                result, output = run_normalizer(Path(temporary), content)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected_error, result.stderr)
                self.assertFalse(output.exists())

    def test_rejects_only_unmapped_duplicate_and_empty_mag_ids(self) -> None:
        header = coverm_header("sample") + "\n"
        cases = (
            (
                "only unmapped",
                header + "unmapped\t100\t0\t0\t1\n",
                "contains no MAG rows",
            ),
            (
                "duplicate MAG",
                header
                + "MAG_1\t10\t1\t0.5\t1000\n"
                + "MAG_1\t20\t2\t0.75\t1000\n",
                "duplicate genome ID: MAG_1",
            ),
            (
                "empty MAG",
                header + "\t10\t1\t0.5\t1000\n",
                "empty genome ID",
            ),
        )
        for label, content, expected_error in cases:
            with self.subTest(case=label), tempfile.TemporaryDirectory() as temporary:
                result, output = run_normalizer(Path(temporary), content)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected_error, result.stderr)
                self.assertFalse(output.exists())

    def test_rejects_non_numeric_non_finite_and_out_of_range_metrics(self) -> None:
        header = coverm_header("sample") + "\n"
        invalid_rows = (
            ("non-numeric", "MAG_1\t10\tnot-a-number\t0.5\t1000", "Non-numeric mean_coverage"),
            ("NaN", "MAG_1\t10\tNaN\t0.5\t1000", "Non-finite mean_coverage"),
            ("infinite", "MAG_1\tInfinity\t1\t0.5\t1000", "Non-finite relative_abundance_percent"),
            ("negative abundance", "MAG_1\t-0.01\t1\t0.5\t1000", "Out-of-range relative_abundance_percent"),
            ("abundance over 100", "MAG_1\t100.01\t1\t0.5\t1000", "Out-of-range relative_abundance_percent"),
            ("negative coverage", "MAG_1\t10\t-0.01\t0.5\t1000", "Out-of-range mean_coverage"),
            ("negative fraction", "MAG_1\t10\t1\t-0.01\t1000", "Out-of-range covered_fraction"),
            ("fraction over one", "MAG_1\t10\t1\t1.01\t1000", "Out-of-range covered_fraction"),
            ("zero length", "MAG_1\t10\t1\t0.5\t0", "Out-of-range genome_length"),
        )
        for label, row, expected_error in invalid_rows:
            with self.subTest(case=label), tempfile.TemporaryDirectory() as temporary:
                result, output = run_normalizer(Path(temporary), header + row + "\n")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected_error, result.stderr)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
