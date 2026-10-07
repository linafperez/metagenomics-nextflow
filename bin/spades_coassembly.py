#!/usr/bin/env python3
"""Validate condition metadata, bind strategy state, and merge SPAdes FASTAs."""

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


def validate_input(args: argparse.Namespace) -> None:
    if not args.group_column:
        raise ValueError("--spades-coassembly-mode condition requires --group-column")
    id_column = "sample" if args.kind == "local" else "biosample_accession"
    if args.group_column == id_column or (args.kind == "local" and args.group_column in {"fastq_1", "fastq_2"}):
        raise ValueError("An input identity/read column cannot also be the group column")
    with args.input.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="," if args.kind == "local" else "\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError("Condition input contains duplicate column names")
        if id_column not in fields or args.group_column not in fields:
            raise ValueError(f"Condition input requires columns {id_column!r} and {args.group_column!r}")
        groups: list[str] = []
        samples: set[str] = set()
        for line, row in enumerate(reader, 2):
            if None in row:
                raise ValueError(f"Condition input line {line} contains unexpected extra fields")
            sample = (row.get(id_column) or "").strip()
            if args.kind == "sra":
                sample = sample.upper()
            if not sample or sample in samples:
                raise ValueError(f"Condition input line {line} has empty or duplicate/ambiguous sample ID {sample!r}")
            samples.add(sample)
            # Same whitespace normalization as the existing input validators.
            groups.append((row.get(args.group_column) or "").strip())
        validate_groups(groups)


def bind_strategy(args: argparse.Namespace) -> None:
    """Fail closed when a results root belongs to a different assembly strategy."""
    root = args.results_dir.resolve()
    info = root / "pipeline_info"
    record = info / "spades_coassembly_strategy.json"
    expected = {"schema_version": 1, "spades_coassembly_mode": args.mode,
                "group_column": args.group_column if args.mode == "condition" else ""}
    if record.exists():
        if json.loads(record.read_text(encoding="utf-8")) != expected:
            raise ValueError("Results root uses a different SPAdes strategy/group column; use a fresh --outdir for the comparison")
        return
    success = info / "sra" / "sra_global_success.json"
    if success.exists():
        old = json.loads(success.read_text(encoding="utf-8"))
        if old.get("spades_coassembly_mode", "global") != args.mode:
            raise ValueError("Existing SRA success marker uses a different SPAdes strategy; use a fresh --outdir")
    if args.mode == "condition" and any(root.glob("02_mag_construction/**/*")):
        raise ValueError("Existing scientific results have no strategy record; legacy results are global. Use a fresh --outdir for condition mode")
    info.mkdir(parents=True, exist_ok=True)
    # Exclusive creation, under the launcher's results lock; never overwrite state.
    with record.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(expected, handle, indent=2, sort_keys=True)
        handle.write("\n")


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
    validate = commands.add_parser("validate-input")
    validate.add_argument("--input", type=Path, required=True)
    validate.add_argument("--kind", choices=["local", "sra"], required=True)
    validate.add_argument("--group-column", required=True)
    bind = commands.add_parser("bind-strategy")
    bind.add_argument("--results-dir", type=Path, required=True)
    bind.add_argument("--mode", choices=["global", "condition"], required=True)
    bind.add_argument("--group-column", default="")
    combine = commands.add_parser("merge")
    combine.add_argument("--manifest", type=Path, required=True)
    combine.add_argument("--output", type=Path, required=True)
    combine.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()
    try:
        {"validate-input": validate_input, "bind-strategy": bind_strategy, "merge": merge}[args.command](args)
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        print(f"ERROR: SPAdes coassembly: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
