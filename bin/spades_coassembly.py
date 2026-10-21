#!/usr/bin/env python3
"""Merge condition SPAdes FASTAs with group-aware contig provenance."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


def group_id(group: str) -> str:
    # Keep identical to SpadesCoassembly.groupId in lib/SpadesCoassembly.groovy.
    return "group_" + re.sub(r"[^A-Za-z0-9._-]+", "_", group)[:80]


def validate_groups(groups: list[str]) -> None:
    if not groups:
        raise ValueError("SPAdes condition coassembly requires at least one sample/group")
    seen: dict[str, str] = {}
    for group in groups:
        if not group.strip() or any(char in group for char in "\t\r\n"):
            raise ValueError("SPAdes condition coassembly requires non-empty groups without tabs or line breaks")
        key = group_id(group).lower()
        if key in seen and seen[key] != group:
            raise ValueError(f"SPAdes group filename collision: {seen[key]!r} and {group!r}")
        seen[key] = group


def merge(args: argparse.Namespace) -> None:
    assemblies = json.loads(args.manifest.read_text(encoding="utf-8"))
    validate_groups([item["group"] for item in assemblies])
    if len({item["group"] for item in assemblies}) != len(assemblies):
        raise ValueError("Merge requires exactly one SPAdes assembly per group")
    combined_ids: set[str] = set()
    with args.output.open("w", encoding="utf-8", newline="\n") as fasta, args.provenance.open(
        "w", encoding="utf-8", newline=""
    ) as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(["combined_contig_id", "original_contig_id", "group", "group_id", "coassembly_id"])
        for item in sorted(assemblies, key=lambda entry: entry["group"]):
            safe = group_id(item["group"])
            if item["group_id"] != safe:
                raise ValueError("Assembly group identifier disagrees with deterministic sanitization")
            count = 0
            sequence_seen = False
            with Path(item["path"]).open(encoding="utf-8") as source:
                for line in source:
                    if line.startswith(">"):
                        if count and not sequence_seen:
                            raise ValueError(f"Empty FASTA contig in {item['path']}")
                        header = line[1:].rstrip("\r\n")
                        parts = header.split(maxsplit=1)
                        if not parts:
                            raise ValueError(f"Empty FASTA header in {item['path']}")
                        original = parts[0]
                        combined = f"{safe}__{original}"
                        if combined in combined_ids:
                            raise ValueError(f"Duplicate combined contig ID: {combined}")
                        combined_ids.add(combined)
                        fasta.write(f">{combined}" + (f" {parts[1]}" if len(parts) == 2 else "") + "\n")
                        writer.writerow([combined, original, item["group"], safe, item["coassembly_id"]])
                        count += 1
                        sequence_seen = False
                    elif line.strip():
                        if not count:
                            raise ValueError(f"Sequence before FASTA header in {item['path']}")
                        fasta.write(line.rstrip("\r\n") + "\n")
                        sequence_seen = True
            if not count or not sequence_seen:
                raise ValueError(f"Empty assembly/contig in {item['path']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    combine = commands.add_parser("merge")
    combine.add_argument("--manifest", type=Path, required=True)
    combine.add_argument("--output", type=Path, required=True)
    combine.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    try:
        merge(args)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        print(f"ERROR: SPAdes coassembly: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
