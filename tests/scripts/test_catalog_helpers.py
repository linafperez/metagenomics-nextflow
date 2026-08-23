#!/usr/bin/env python3
"""Unit tests for catalog construction helper scripts."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).parents[2]
BIN_DIR = PROJECT_DIR / "bin"


def run_script(name: str, *arguments: object) -> None:
    command = [sys.executable, str(BIN_DIR / name), *(str(value) for value in arguments)]
    subprocess.run(command, check=True, capture_output=True, text=True)


def write_fasta(path: Path, identifier: str, length: int) -> None:
    path.write_text(f">{identifier}\n{'ACGT' * (length // 4)}\n", encoding="utf-8")


class CatalogHelperTests(unittest.TestCase):
    def test_normalize_coverm_abundance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wide = root / "wide.tsv"
            wide.write_text(
                "Genome\tsample_A Relative Abundance (%)\tsample_A Mean\t"
                "sample_A Covered Fraction\tsample_A Length\n"
                "MAG_001\t50\t12.5\t0.8\t2400\n",
                encoding="utf-8",
            )
            long_table = root / "long.tsv"
            run_script(
                "normalize_coverm_abundance.py",
                "--input",
                wide,
                "--output",
                long_table,
            )
            with long_table.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["sample"], "sample_A")
            self.assertEqual(rows[0]["covered_fraction"], "0.8")

    def test_filter_and_strict_quality_selection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "assembly.fa"
            assembly.write_text(
                f">short\n{'A' * 100}\n>long\n{'C' * 1600}\n", encoding="utf-8"
            )
            run_script(
                "filter_fasta_by_length.py",
                "--input",
                assembly,
                "--output",
                root / "filtered.fa",
                "--stats",
                root / "stats.tsv",
                "--min-length",
                1500,
            )
            self.assertIn(">long", (root / "filtered.fa").read_text(encoding="utf-8"))
            self.assertNotIn(">short", (root / "filtered.fa").read_text(encoding="utf-8"))

            bins = []
            for identifier in ("pass", "equal_completeness", "equal_contamination"):
                path = root / f"{identifier}.fa"
                write_fasta(path, identifier, 1600)
                bins.append(path)
            quality = root / "quality.tsv"
            quality.write_text(
                "Name\tCompleteness\tContamination\n"
                "pass\t91\t4\n"
                "equal_completeness\t90\t1\n"
                "equal_contamination\t99\t5\n",
                encoding="utf-8",
            )
            run_script(
                "select_high_quality_mags.py",
                "--bins",
                *bins,
                "--quality",
                quality,
                "--output-dir",
                root / "selected",
                "--selected-table",
                root / "selected.tsv",
                "--genome-info",
                root / "genomeInfo.csv",
                "--completeness",
                90,
                "--contamination",
                5,
                "--assembler",
                "megahit",
            )
            self.assertEqual([path.name for path in (root / "selected").glob("*.fa")], ["pass.fa"])

    def test_combine_and_finalize_catalog(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            megahit = root / "megahit_bin.fa"
            spades = root / "spades_bin.fa"
            write_fasta(megahit, "megahit_bin", 1600)
            write_fasta(spades, "spades_bin", 1600)
            megahit_quality = root / "megahit.tsv"
            spades_quality = root / "spades.tsv"
            megahit_quality.write_text(
                "Name\tCompleteness\tContamination\nmegahit_bin\t96\t1\n",
                encoding="utf-8",
            )
            spades_quality.write_text(
                "Name\tCompleteness\tContamination\nspades_bin\t95\t2\n",
                encoding="utf-8",
            )
            combined = root / "combined"
            provenance = root / "provenance.tsv"
            quality = root / "quality.tsv"
            run_script(
                "combine_mag_catalogs.py",
                "--megahit-bins",
                megahit,
                "--spades-bins",
                spades,
                "--megahit-quality",
                megahit_quality,
                "--spades-quality",
                spades_quality,
                "--output-dir",
                combined,
                "--provenance",
                provenance,
                "--quality-table",
                quality,
                "--genome-info",
                root / "genomeInfo.csv",
            )
            self.assertEqual(
                sorted(path.name for path in combined.glob("*.fa")),
                ["MAG_MEGAHIT_0001.fa", "MAG_SPADES_0001.fa"],
            )

            representative = combined / "MAG_MEGAHIT_0001.fa"
            run_script(
                "finalize_mag_catalog.py",
                "--representatives",
                representative,
                "--provenance",
                provenance,
                "--quality",
                quality,
                "--output-dir",
                root / "final",
                "--output-provenance",
                root / "final.provenance.tsv",
                "--output-quality",
                root / "final.quality.tsv",
            )
            with (root / "final.quality.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["mag_id"], "MAG_MEGAHIT_0001")


if __name__ == "__main__":
    unittest.main()
