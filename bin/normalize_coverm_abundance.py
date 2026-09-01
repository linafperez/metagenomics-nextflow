#!/usr/bin/env python3
"""Convert a wide CoverM genome table to analysis-ready long format."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path


METRICS = {
    "Relative Abundance (%)": "relative_abundance_percent",
    "Mean": "mean_coverage",
    "Covered Fraction": "covered_fraction",
    "Length": "genome_length",
}
OUTPUT_FIELDS = (
    "sample",
    "mag_id",
    "relative_abundance_percent",
    "mean_coverage",
    "covered_fraction",
    "genome_length",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def parse_columns(fieldnames: list[str]) -> tuple[str, dict[str, dict[str, str]]]:
    if len(set(fieldnames)) != len(fieldnames):
        raise ValueError("CoverM abundance table has duplicate columns")
    genome_column = fieldnames[0]
    if not genome_column.strip():
        raise ValueError("CoverM abundance table has an empty genome column")
    samples: dict[str, dict[str, str]] = {}
    for column in fieldnames[1:]:
        matched = False
        for suffix, metric in METRICS.items():
            marker = f" {suffix}"
            if column.endswith(marker):
                sample = column[: -len(marker)]
                if not sample:
                    raise ValueError(f"CoverM abundance column has no sample name: {column}")
                sample_columns = samples.setdefault(sample, {})
                if metric in sample_columns:
                    raise ValueError(
                        f"Sample '{sample}' has duplicate metric column: {suffix}"
                    )
                sample_columns[metric] = column
                matched = True
                break
        if not matched:
            raise ValueError(f"Unsupported CoverM abundance column: {column}")

    if not samples:
        raise ValueError("CoverM abundance table contains no sample metric columns")
    expected = set(METRICS.values())
    for sample, columns in samples.items():
        missing = expected - set(columns)
        if missing:
            raise ValueError(f"Sample '{sample}' lacks metrics: {', '.join(sorted(missing))}")
    return genome_column, samples


def validate_metric(value: str | None, metric: str, sample: str, mag_id: str) -> None:
    raw = (value or "").strip()
    try:
        number = Decimal(raw)
    except InvalidOperation as exc:
        raise ValueError(
            f"Non-numeric {metric} for sample '{sample}', MAG '{mag_id}': {raw!r}"
        ) from exc
    if not number.is_finite():
        raise ValueError(
            f"Non-finite {metric} for sample '{sample}', MAG '{mag_id}': {raw!r}"
        )
    if metric == "genome_length":
        valid = number > 0
    elif metric == "covered_fraction":
        valid = Decimal(0) <= number <= Decimal(1)
    elif metric == "relative_abundance_percent":
        valid = Decimal(0) <= number <= Decimal(100)
    else:
        valid = number >= 0
    if not valid:
        raise ValueError(
            f"Out-of-range {metric} for sample '{sample}', MAG '{mag_id}': {raw!r}"
        )


def main() -> int:
    args = parse_args()
    with args.input.open(newline="", encoding="utf-8-sig") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("CoverM abundance table has no header")
        genome_column, samples = parse_columns(reader.fieldnames)
        rows = list(reader)

    output_rows: list[dict[str, str]] = []
    seen_mags: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        if None in row:
            raise ValueError(f"CoverM abundance row {row_number} has extra fields")
        mag_id = (row.get(genome_column) or "").strip()
        # CoverM emits this pseudo-genome for the fraction of reads not mapped
        # to the supplied catalog.  Its percentage remains reflected in the
        # real MAG values, but it is not itself a MAG abundance row.
        if mag_id == "unmapped":
            continue
        if not mag_id:
            raise ValueError(f"CoverM abundance row {row_number} has an empty genome ID")
        if mag_id in seen_mags:
            raise ValueError(f"CoverM abundance table has duplicate genome ID: {mag_id}")
        seen_mags.add(mag_id)
        for sample in sorted(samples):
            columns = samples[sample]
            for metric, column in columns.items():
                validate_metric(row.get(column), metric, sample, mag_id)
            output_rows.append(
                {
                    "sample": sample,
                    "mag_id": mag_id,
                    **{metric: row[column] for metric, column in columns.items()},
                }
            )

    if not seen_mags:
        raise ValueError("CoverM abundance table contains no MAG rows")

    with args.output.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.DictWriter(
            output_handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
