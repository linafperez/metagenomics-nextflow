from __future__ import annotations

import os
import re
import shutil
import subprocess
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
AWK = os.environ.get("AWK") or shutil.which("awk")


def coverage_program(source: str | None = None) -> str:
    """Extract the MAG AWK and decode its real Groovy string escapes only."""
    if source is None:
        source = (REPOSITORY / "modules/core/coverm/main.nf").read_text(encoding="utf-8")
    genome_process = source.split("process COVERM_GENOME {", 1)[1]
    match = re.search(
        r"(?ms)^\s*awk '\n(.*?)^\s*' \"\$\{prefix\}\.mag_abundance\.raw\.tsv\"",
        genome_process,
    )
    if match is None:
        raise AssertionError("COVERM_GENOME MAG abundance AWK program not found")
    escapes = {"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"', "$": "$"}
    return re.sub(r'\\([ntr\\"$])', lambda item: escapes[item[1]], match[1])


def header(*samples: str) -> str:
    return "\t".join([
        "Genome",
        *(f"{sample} {metric}" for sample in samples for metric in (
            "Relative Abundance (%)", "Mean", "Covered Bases", "Length",
        )),
    ])


@unittest.skipUnless(AWK, "awk is required (or set AWK to its executable path)")
class CovermMagAbundanceAwkTests(unittest.TestCase):
    def run_fragment(self, content: str) -> list[list[str]]:
        result = subprocess.run(
            [AWK, coverage_program()], input=content, capture_output=True, text=True, timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        return [line.split("\t") for line in result.stdout.splitlines()]

    def test_unmapped_na_and_numeric_mags_with_one_two_or_eight_samples(self) -> None:
        for count in (1, 2, 8):
            with self.subTest(samples=count):
                samples = [f"sample_{idx}" for idx in range(count)]
                unmapped = ["unmapped", *(["60.0", "NA", "NA", "NA"] * count)]
                mag = ["MAG_TEST_0001"]
                for idx in range(count):
                    mag.extend([f"{idx + 1}.0", "2.5", str(900 - idx * 100), "1000"])
                zero = ["MAG_ZERO", *(["0.0", "0.0", "0", "1000"] * count)]
                raw_rows = [unmapped, mag, zero]
                content = header(*samples) + "\n" + "".join("\t".join(row) + "\n" for row in raw_rows)
                result = self.run_fragment(content)
                self.assertEqual(result[0], header(*samples).replace("Covered Bases", "Covered Fraction").split("\t"))
                self.assertEqual(len(result), 4)
                self.assertEqual(result[1], unmapped)
                for before, after in zip(raw_rows, result[1:]):
                    self.assertEqual(len(after), 1 + 4 * count)
                    self.assertEqual(len(before), len(after))
                    for column in range(len(before)):
                        if column not in range(3, len(before), 4):
                            self.assertEqual(after[column], before[column])
                self.assertEqual(result[2][3::4], [f"0.{9 - idx}" for idx in range(count)])
                self.assertEqual(result[3][3::4], ["0"] * count)

    def test_missing_or_nonpositive_length_keeps_na_and_valid_zero_stays_zero(self) -> None:
        cases = (
            ("NA", "NA", "NA"),
            ("NA", "1000", "NA"),
            ("", "1000", "NA"),
            ("900", "NA", "NA"),
            ("900", "", "NA"),
            ("900", "0", "NA"),
            ("0", "0", "NA"),
            ("900", "-100", "NA"),
            ("900", "not-a-number", "NA"),
            ("0", "1000", "0"),
            ("0.0", "1000.0", "0"),
            ("9e2", "1e3", "0.9"),
        )
        raw_rows = [[f"MAG_{idx}", "1.0", "2.5", covered, length] for idx, (covered, length, _) in enumerate(cases)]
        content = header("sample") + "\n" + "".join("\t".join(row) + "\n" for row in raw_rows)
        result = self.run_fragment(content)
        self.assertEqual(len(result), len(cases) + 1)
        for before, after, (_, _, fraction) in zip(raw_rows, result[1:], cases):
            self.assertEqual(after, [*before[:3], fraction, before[4]])

    def test_hpc_mag_values_keep_covered_bases_divided_by_length(self) -> None:
        raw = ["MAG_MEGAHIT_0003", "0.2167008", "2.996999", "1785520", "1976380"]
        result = self.run_fragment(header("sample") + "\n" + "\t".join(raw) + "\n")
        self.assertEqual(result[1][:3], raw[:3])
        self.assertEqual(result[1][4], raw[4])
        self.assertAlmostEqual(float(result[1][3]), 1785520 / 1976380, places=6)


if __name__ == "__main__":
    unittest.main()
