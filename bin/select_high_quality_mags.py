#!/usr/bin/env python3
"""Select high-quality MAGs from a CheckM2 quality report."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bins", nargs="+", required=True, type=Path)
    parser.add_argument("--quality", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--selected-table", required=True, type=Path)
    parser.add_argument("--genome-info", required=True, type=Path)
    parser.add_argument("--completeness", required=True, type=float)
    parser.add_argument("--contamination", required=True, type=float)
    parser.add_argument("--assembler", required=True)
    return parser.parse_args()


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def find_column(fieldnames: list[str], *candidates: str) -> str:
    lookup = {normalized(field): field for field in fieldnames}
    for candidate in candidates:
        if normalized(candidate) in lookup:
            return lookup[normalized(candidate)]
    raise ValueError(f"Missing required quality column: {', '.join(candidates)}")


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
            raise ValueError(f"Quality report has no header: {path}")
        name_column = find_column(reader.fieldnames, "Name", "genome", "mag_id")
        completeness_column = find_column(reader.fieldnames, "Completeness")
        contamination_column = find_column(reader.fieldnames, "Contamination")
        quality: dict[str, tuple[float, float]] = {}
        for row in reader:
            name = fasta_stem(row[name_column])
            quality[name] = (
                float(row[completeness_column]),
                float(row[contamination_column]),
            )
    return quality


def main() -> int:
    args = parse_args()
    quality = load_quality(args.quality)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected: list[tuple[str, str, float, float]] = []

    for bin_path in sorted(args.bins, key=lambda path: path.name):
        original_id = fasta_stem(bin_path.name)
        if original_id not in quality:
            raise ValueError(f"No CheckM2 row found for bin '{bin_path.name}'")
        completeness, contamination = quality[original_id]
        if completeness <= args.completeness or contamination >= args.contamination:
            continue
        output_name = f"{original_id}.fa"
        destination = args.output_dir / output_name
        if destination.exists():
            raise ValueError(f"Selected MAG filename collision: {output_name}")
        shutil.copyfile(bin_path, destination)
        selected.append((output_name, original_id, completeness, contamination))

    if not selected:
        raise ValueError(
            "No MAGs satisfy completeness > "
            f"{args.completeness} and contamination < {args.contamination}"
        )

    with args.selected_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("genome", "original_mag_id", "assembler", "completeness", "contamination")
        )
        for genome, original_id, completeness, contamination in selected:
            writer.writerow(
                (genome, original_id, args.assembler, completeness, contamination)
            )

    with args.genome_info.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("genome", "completeness", "contamination"))
        for genome, _original_id, completeness, contamination in selected:
            writer.writerow((genome, completeness, contamination))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
