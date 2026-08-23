#!/usr/bin/env python3
"""Merge assembler-specific MAG catalogs with collision-free identifiers."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--megahit-bins", nargs="+", required=True, type=Path)
    parser.add_argument("--spades-bins", nargs="+", required=True, type=Path)
    parser.add_argument("--megahit-quality", required=True, type=Path)
    parser.add_argument("--spades-quality", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--quality-table", required=True, type=Path)
    parser.add_argument("--genome-info", required=True, type=Path)
    return parser.parse_args()


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def fasta_stem(path_or_name: str) -> str:
    name = Path(path_or_name).name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def load_quality(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Quality table has no header: {path}")
        fields = {normalized(field): field for field in reader.fieldnames}
        name_column = fields.get("name") or fields.get("genome") or fields.get("magid")
        completeness_column = fields.get("completeness")
        contamination_column = fields.get("contamination")
        if not name_column or not completeness_column or not contamination_column:
            raise ValueError(f"Unsupported quality table header: {path}")
        return {
            fasta_stem(row[name_column]): (
                float(row[completeness_column]),
                float(row[contamination_column]),
            )
            for row in reader
        }


def main() -> int:
    args = parse_args()
    quality_by_branch = {
        "megahit": load_quality(args.megahit_quality),
        "spades": load_quality(args.spades_quality),
    }
    inputs = {
        "megahit": sorted(args.megahit_bins, key=lambda path: path.name),
        "spades": sorted(args.spades_bins, key=lambda path: path.name),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, str | float]] = []
    for assembler in ("megahit", "spades"):
        for index, source in enumerate(inputs[assembler], start=1):
            original_id = fasta_stem(source.name)
            if original_id not in quality_by_branch[assembler]:
                raise ValueError(
                    f"No clean CheckM2 quality row for {assembler} MAG '{source.name}'"
                )
            completeness, contamination = quality_by_branch[assembler][original_id]
            mag_id = f"MAG_{assembler.upper()}_{index:04d}"
            output_name = f"{mag_id}.fa"
            shutil.copyfile(source, args.output_dir / output_name)
            records.append(
                {
                    "mag_id": mag_id,
                    "genome": output_name,
                    "assembler": assembler,
                    "original_mag_id": original_id,
                    "source_file": source.name,
                    "completeness": completeness,
                    "contamination": contamination,
                }
            )

    if not records:
        raise ValueError("No MAGs were provided for final catalog construction")

    with args.provenance.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
            fieldnames=("mag_id", "genome", "assembler", "original_mag_id", "source_file"),
        )
        writer.writeheader()
        writer.writerows(records)

    with args.quality_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
            fieldnames=(
                "mag_id",
                "genome",
                "assembler",
                "original_mag_id",
                "completeness",
                "contamination",
            ),
        )
        writer.writeheader()
        writer.writerows(records)

    with args.genome_info.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("genome", "completeness", "contamination"))
        for record in records:
            writer.writerow(
                (record["genome"], record["completeness"], record["contamination"])
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
