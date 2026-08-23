#!/usr/bin/env python3
"""Convert a wide CoverM genome table to analysis-ready long format."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


METRICS = {
    "Relative Abundance (%)": "relative_abundance_percent",
    "Mean": "mean_coverage",
    "Covered Fraction": "covered_fraction",
    "Length": "genome_length",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_columns(fieldnames: list[str]) -> tuple[str, dict[str, dict[str, str]]]:
    genome_column = fieldnames[0]
    samples: dict[str, dict[str, str]] = {}
    for column in fieldnames[1:]:
        matched = False
        for suffix, metric in METRICS.items():
            marker = f" {suffix}"
            if column.endswith(marker):
                sample = column[: -len(marker)]
                samples.setdefault(sample, {})[metric] = column
                matched = True
                break
        if not matched:
            raise ValueError(f"Unsupported CoverM abundance column: {column}")

    expected = set(METRICS.values())
    for sample, columns in samples.items():
        missing = expected - set(columns)
        if missing:
            raise ValueError(f"Sample '{sample}' lacks metrics: {', '.join(sorted(missing))}")
    return genome_column, samples


def main() -> int:
    args = parse_args()
    with args.input.open(newline="", encoding="utf-8-sig") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("CoverM abundance table has no header")
        genome_column, samples = parse_columns(reader.fieldnames)
        rows = list(reader)

    with args.output.open("w", newline="", encoding="utf-8") as output_handle:
        fieldnames = (
            "sample",
            "mag_id",
            "relative_abundance_percent",
            "mean_coverage",
            "covered_fraction",
            "genome_length",
        )
        writer = csv.DictWriter(
            output_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            for sample in sorted(samples):
                columns = samples[sample]
                writer.writerow(
                    {
                        "sample": sample,
                        "mag_id": row[genome_column],
                        **{metric: row[column] for metric, column in columns.items()},
                    }
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
