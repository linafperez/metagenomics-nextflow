#!/usr/bin/env python3
"""Materialize final representatives and subset catalog metadata."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives", nargs="+", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--quality", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-provenance", required=True, type=Path)
    parser.add_argument("--output-quality", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected_names = {path.name for path in args.representatives}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for representative in sorted(args.representatives, key=lambda path: path.name):
        shutil.copyfile(representative, args.output_dir / representative.name)

    def subset(source: Path, destination: Path) -> int:
        with source.open(newline="", encoding="utf-8-sig") as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            if not reader.fieldnames or "genome" not in reader.fieldnames:
                raise ValueError(f"Metadata table lacks a genome column: {source}")
            rows = [row for row in reader if row["genome"] in selected_names]
        with destination.open("w", newline="", encoding="utf-8") as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=reader.fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        return len(rows)

    provenance_count = subset(args.provenance, args.output_provenance)
    quality_count = subset(args.quality, args.output_quality)
    if provenance_count != len(selected_names) or quality_count != len(selected_names):
        raise ValueError(
            "Final catalog metadata does not match the representative FASTA set: "
            f"FASTA={len(selected_names)}, provenance={provenance_count}, quality={quality_count}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
