#!/usr/bin/env python3
"""Validate and normalize the paired-end FASTQ samplesheet."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


REQUIRED_HEADER = ["sample", "fastq_1", "fastq_2"]
NORMALIZED_HEADER = [*REQUIRED_HEADER, "group"]
METADATA_HEADER = [
    "sample_id",
    "biosample_accession",
    "group",
    "run_accessions",
    "project_accession",
    "sample_order",
]
FASTQ_EXTENSIONS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
SAMPLE_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a paired-end FASTQ samplesheet."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--group-column", default="")
    return parser.parse_args()


def resolve_fastq(raw_path: str, samplesheet_dir: Path) -> Path:
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = samplesheet_dir / path
    return path.resolve()


def validate_samplesheet(
    samplesheet: Path, group_column: str = ""
) -> list[dict[str, str]]:
    errors: list[str] = []
    normalized_rows: list[dict[str, str]] = []
    seen_samples: set[str] = set()

    with samplesheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames or reader.fieldnames[:3] != REQUIRED_HEADER:
            raise ValueError(
                "Samplesheet header must begin with: " + ",".join(REQUIRED_HEADER)
            )
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Samplesheet contains duplicate column names")
        if group_column and group_column not in reader.fieldnames:
            raise ValueError(
                f"Configured group column {group_column!r} is absent from the samplesheet"
            )
        if group_column in REQUIRED_HEADER:
            raise ValueError(f"{group_column!r} cannot be used as the group column")

        for line_number, row in enumerate(reader, start=2):
            if None in row:
                errors.append(
                    f"Line {line_number}: row contains unexpected extra fields"
                )
                continue
            sample = (row.get("sample") or "").strip()
            fastq_1_value = (row.get("fastq_1") or "").strip()
            fastq_2_value = (row.get("fastq_2") or "").strip()
            group = (row.get(group_column) or "").strip() if group_column else ""

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

            if group_column and not group:
                errors.append(
                    f"Line {line_number}: group column {group_column!r} is empty"
                )
            if any(character in group for character in "\t\r\n"):
                errors.append(
                    f"Line {line_number}: group contains a tab or line break"
                )

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
                    "group": group,
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
        rows = validate_samplesheet(samplesheet, args.group_column)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_HEADER)
        writer.writeheader()
        writer.writerows(rows)

    with args.metadata_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=METADATA_HEADER,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for order, row in enumerate(rows, start=1):
            writer.writerow(
                {
                    "sample_id": row["sample"],
                    "biosample_accession": "",
                    "group": row["group"],
                    "run_accessions": "",
                    "project_accession": "",
                    "sample_order": str(order),
                }
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
