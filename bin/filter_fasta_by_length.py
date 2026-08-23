#!/usr/bin/env python3
"""Filter FASTA records by minimum sequence length."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--min-length", required=True, type=int)
    return parser.parse_args()


def fasta_records(path: Path):
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line
                sequence = []
            elif header is None:
                raise ValueError(f"Sequence found before FASTA header in {path}")
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def main() -> int:
    args = parse_args()
    if args.min_length < 1:
        raise ValueError("Minimum length must be positive")

    total_records = 0
    total_bases = 0
    kept_records = 0
    kept_bases = 0

    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        for header, sequence in fasta_records(args.input):
            total_records += 1
            total_bases += len(sequence)
            if len(sequence) < args.min_length:
                continue
            kept_records += 1
            kept_bases += len(sequence)
            output.write(header + "\n")
            for offset in range(0, len(sequence), 80):
                output.write(sequence[offset : offset + 80] + "\n")

    if kept_records == 0:
        args.output.unlink(missing_ok=True)
        raise ValueError(
            f"No contigs in {args.input} meet the {args.min_length} bp threshold"
        )

    with args.stats.open("w", encoding="utf-8", newline="\n") as stats:
        stats.write("metric\tvalue\n")
        stats.write(f"minimum_length\t{args.min_length}\n")
        stats.write(f"input_contigs\t{total_records}\n")
        stats.write(f"input_bases\t{total_bases}\n")
        stats.write(f"retained_contigs\t{kept_records}\n")
        stats.write(f"retained_bases\t{kept_bases}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
