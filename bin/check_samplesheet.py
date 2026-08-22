#!/usr/bin/env python3
"""Validate and normalize the paired-end FASTQ samplesheet."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


EXPECTED_HEADER = ["sample", "fastq_1", "fastq_2"]
FASTQ_EXTENSIONS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
SAMPLE_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a paired-end FASTQ samplesheet."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def resolve_fastq(raw_path: str, samplesheet_dir: Path) -> Path:
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = samplesheet_dir / path
    return path.resolve()


def validate_samplesheet(samplesheet: Path) -> list[dict[str, str]]:
    errors: list[str] = []
    normalized_rows: list[dict[str, str]] = []
    seen_samples: set[str] = set()

    with samplesheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames != EXPECTED_HEADER:
            raise ValueError(
                "Samplesheet header must be exactly: " + ",".join(EXPECTED_HEADER)
            )

        for line_number, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            fastq_1_value = (row.get("fastq_1") or "").strip()
            fastq_2_value = (row.get("fastq_2") or "").strip()

            if not sample:
                errors.append(f"Line {line_number}: sample is empty")
                continue

            if not SAMPLE_PATTERN.fullmatch(sample):
                errors.append(
                    f"Line {line_number}: sample '{sample}' contains unsupported characters"
                )

            if sample in seen_samples:
                errors.append(f"Line {line_number}: duplicate sample ID '{sample}'")
            seen_samples.add(sample)

            if not fastq_1_value or not fastq_2_value:
                errors.append(
                    f"Line {line_number}: both fastq_1 and fastq_2 are required"
                )
                continue

            fastq_1 = resolve_fastq(fastq_1_value, samplesheet.parent)
            fastq_2 = resolve_fastq(fastq_2_value, samplesheet.parent)

            for field_name, fastq in (("fastq_1", fastq_1), ("fastq_2", fastq_2)):
                if not str(fastq).lower().endswith(FASTQ_EXTENSIONS):
                    errors.append(
                        f"Line {line_number}: {field_name} has an unsupported extension: {fastq}"
                    )
                if not fastq.is_file():
                    errors.append(
                        f"Line {line_number}: {field_name} does not exist or is not a file: {fastq}"
                    )

            if fastq_1 == fastq_2:
                errors.append(
                    f"Line {line_number}: fastq_1 and fastq_2 resolve to the same file"
                )

            normalized_rows.append(
                {
                    "sample": sample,
                    "fastq_1": str(fastq_1),
                    "fastq_2": str(fastq_2),
                }
            )

    if not normalized_rows:
        errors.append("Samplesheet contains no sample rows")

    if errors:
        raise ValueError("\n".join(errors))

    return normalized_rows


def main() -> int:
    args = parse_args()
    samplesheet = args.input.resolve()

    try:
        rows = validate_samplesheet(samplesheet)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPECTED_HEADER)
        writer.writeheader()
        writer.writerows(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
