#!/usr/bin/env python3
"""Integrate eggNOG-mapper and InterProScan annotations per predicted protein."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO


MISSING_VALUES = {"", "-", "NA", "N/A", "none", "None", "null"}

OUTPUT_COLUMNS = (
    "mag_id",
    "protein_id",
    "gene_id",
    "sequence_length_aa",
    "contig",
    "feature_start",
    "feature_end",
    "strand",
    "preferred_name",
    "description",
    "seed_ortholog",
    "eggnog_ogs",
    "max_annotation_level",
    "cog_category",
    "go_terms",
    "ec_numbers",
    "kegg_ko",
    "kegg_pathways",
    "kegg_modules",
    "kegg_reactions",
    "cazy_families",
    "pfam_accessions",
    "interpro_accessions",
    "interpro_descriptions",
    "interpro_member_databases",
    "interpro_signature_accessions",
    "interpro_signature_descriptions",
    "interpro_pathways",
    "annotation_sources",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteins", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--eggnog", required=True, type=Path)
    parser.add_argument("--interpro", required=True, type=Path)
    parser.add_argument("--mag-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lstrip("#").lower())


def present(value: str | None) -> str:
    if value is None:
        return ""
    stripped = value.strip()
    return "" if stripped in MISSING_VALUES else stripped


def split_values(value: str | None, separator_pattern: str = r"[,;|]") -> set[str]:
    cleaned = present(value)
    if not cleaned:
        return set()
    return {
        token.strip()
        for token in re.split(separator_pattern, cleaned)
        if present(token)
    }


def joined(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if present(value)}))


def first_value(row: dict[str, str], *candidates: str) -> str:
    lookup = {normalized(key): value for key, value in row.items()}
    for candidate in candidates:
        value = present(lookup.get(normalized(candidate)))
        if value:
            return value
    return ""


def parse_fasta(path: Path) -> tuple[list[str], dict[str, int]]:
    identifiers: list[str] = []
    lengths: dict[str, int] = {}
    current: str | None = None

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                if identifier in lengths:
                    raise ValueError(f"Duplicate protein identifier '{identifier}' in {path}")
                identifiers.append(identifier)
                lengths[identifier] = 0
                current = identifier
                continue
            if current is None:
                raise ValueError(f"Sequence precedes the first FASTA header in {path}")
            lengths[current] += len(line.replace(" ", ""))

    if not identifiers:
        raise ValueError(f"Protein FASTA contains no records: {path}")
    return identifiers, lengths


def parse_attributes(raw_attributes: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for raw_part in raw_attributes.split(";"):
        part = raw_part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", maxsplit=1)
        elif re.search(r"\s", part):
            key, value = part.split(maxsplit=1)
        else:
            continue
        attributes[key.strip()] = value.strip().strip('"')
    return attributes


def parse_gff(path: Path) -> dict[str, dict[str, str]]:
    features: dict[str, dict[str, str]] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"Expected nine GFF columns at {path}:{line_number}")
            if fields[2].lower() not in {"cds", "gene"}:
                continue
            attributes = parse_attributes(fields[8])
            feature = {
                "gene_id": present(
                    attributes.get("gene_id")
                    or attributes.get("ID")
                    or attributes.get("locus_tag")
                ),
                "contig": fields[0],
                "start": fields[3],
                "end": fields[4],
                "strand": fields[6],
            }
            aliases = {
                attributes.get("protein_id"),
                attributes.get("ID"),
                attributes.get("gene_id"),
                attributes.get("locus_tag"),
                attributes.get("Name"),
            }
            for alias in aliases:
                if present(alias):
                    features[present(alias)] = feature
    return features


def parse_eggnog(path: Path) -> dict[str, dict[str, str]]:
    annotations: dict[str, dict[str, str]] = {}
    header: list[str] | None = None

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                candidate = line.lstrip("#").split("\t")
                if candidate and normalized(candidate[0]) in {"query", "queryname"}:
                    header = candidate
                continue
            if header is None:
                raise ValueError(f"eggNOG annotations have no header before {path}:{line_number}")
            values = line.split("\t")
            values.extend([""] * (len(header) - len(values)))
            row = dict(zip(header, values, strict=False))
            query = first_value(row, "query", "query_name")
            if not query:
                raise ValueError(f"eggNOG row has no query identifier at {path}:{line_number}")
            annotations[query] = row

    if header is None:
        raise ValueError(f"eggNOG annotations contain no recognized header: {path}")
    return annotations


def parse_interpro(path: Path) -> dict[str, dict[str, set[str]]]:
    annotations: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 13:
                raise ValueError(
                    f"Expected at least 13 InterProScan TSV columns at {path}:{line_number}"
                )
            protein_id = present(fields[0])
            if not protein_id:
                raise ValueError(f"InterProScan row has no protein identifier at {path}:{line_number}")
            record = annotations[protein_id]
            analysis = present(fields[3])
            signature_accession = present(fields[4])
            signature_description = present(fields[5])
            interpro_accession = present(fields[11])
            interpro_description = present(fields[12])
            record["analyses"].update({analysis} if analysis else set())
            record["signatures"].update(
                {signature_accession} if signature_accession else set()
            )
            record["signature_descriptions"].update(
                {signature_description} if signature_description else set()
            )
            record["interpro_accessions"].update(
                {interpro_accession} if interpro_accession else set()
            )
            record["interpro_descriptions"].update(
                {interpro_description} if interpro_description else set()
            )
            if len(fields) > 13:
                record["go_terms"].update(split_values(fields[13]))
            if len(fields) > 14:
                record["pathways"].update(split_values(fields[14]))
            if analysis.lower() == "pfam" or signature_accession.upper().startswith("PF"):
                record["pfam"].update(
                    {signature_accession} if signature_accession else set()
                )
    return annotations


def validate_annotation_ids(
    protein_ids: set[str], source_name: str, annotation_ids: Iterable[str]
) -> None:
    unknown = sorted(set(annotation_ids) - protein_ids)
    if unknown:
        preview = ", ".join(unknown[:5])
        raise ValueError(
            f"{source_name} contains identifiers absent from the protein FASTA: {preview}"
        )


def build_rows(
    mag_id: str,
    identifiers: list[str],
    lengths: dict[str, int],
    features: dict[str, dict[str, str]],
    eggnog: dict[str, dict[str, str]],
    interpro: dict[str, dict[str, set[str]]],
) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for protein_id in identifiers:
        egg = eggnog.get(protein_id, {})
        ipr = interpro.get(protein_id, {})
        feature = features.get(protein_id, {})
        egg_go = split_values(first_value(egg, "GOs", "GO"), r"[,;]")
        egg_pfam = split_values(first_value(egg, "PFAMs", "PFAM"), r"[,;]")
        sources = []
        if egg:
            sources.append("eggNOG-mapper")
        if ipr:
            sources.append("InterProScan")

        rows.append(
            {
                "mag_id": mag_id,
                "protein_id": protein_id,
                "gene_id": feature.get("gene_id", ""),
                "sequence_length_aa": lengths[protein_id],
                "contig": feature.get("contig", ""),
                "feature_start": feature.get("start", ""),
                "feature_end": feature.get("end", ""),
                "strand": feature.get("strand", ""),
                "preferred_name": first_value(egg, "Preferred_name", "preferred_name"),
                "description": first_value(egg, "Description"),
                "seed_ortholog": first_value(egg, "seed_ortholog"),
                "eggnog_ogs": joined(
                    split_values(
                        first_value(egg, "eggNOG_OGs", "eggnog_ogs"), r"[,;]"
                    )
                ),
                "max_annotation_level": first_value(egg, "max_annot_lvl"),
                "cog_category": joined(
                    split_values(first_value(egg, "COG_category", "cog"), r"[,;]")
                ),
                "go_terms": joined(egg_go | ipr.get("go_terms", set())),
                "ec_numbers": joined(
                    split_values(first_value(egg, "EC"), r"[,;]")
                ),
                "kegg_ko": joined(
                    split_values(first_value(egg, "KEGG_ko"), r"[,;]")
                ),
                "kegg_pathways": joined(
                    split_values(first_value(egg, "KEGG_Pathway"), r"[,;]")
                ),
                "kegg_modules": joined(
                    split_values(first_value(egg, "KEGG_Module"), r"[,;]")
                ),
                "kegg_reactions": joined(
                    split_values(first_value(egg, "KEGG_Reaction"), r"[,;]")
                ),
                "cazy_families": joined(
                    split_values(first_value(egg, "CAZy"), r"[,;]")
                ),
                "pfam_accessions": joined(egg_pfam | ipr.get("pfam", set())),
                "interpro_accessions": joined(ipr.get("interpro_accessions", set())),
                "interpro_descriptions": joined(ipr.get("interpro_descriptions", set())),
                "interpro_member_databases": joined(ipr.get("analyses", set())),
                "interpro_signature_accessions": joined(ipr.get("signatures", set())),
                "interpro_signature_descriptions": joined(
                    ipr.get("signature_descriptions", set())
                ),
                "interpro_pathways": joined(ipr.get("pathways", set())),
                "annotation_sources": ";".join(sources),
            }
        )
    return rows


def write_table(path: Path, rows: list[dict[str, str | int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=OUTPUT_COLUMNS,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, mag_id: str, rows: list[dict[str, str | int]]) -> None:
    summary = {
        "mag_id": mag_id,
        "protein_count": len(rows),
        "eggnog_annotated_proteins": sum(
            "eggNOG-mapper" in str(row["annotation_sources"]) for row in rows
        ),
        "interpro_annotated_proteins": sum(
            "InterProScan" in str(row["annotation_sources"]) for row in rows
        ),
        "proteins_with_go": sum(bool(row["go_terms"]) for row in rows),
        "proteins_with_kegg": sum(
            any(
                row[field]
                for field in (
                    "kegg_ko",
                    "kegg_pathways",
                    "kegg_modules",
                    "kegg_reactions",
                )
            )
            for row in rows
        ),
        "proteins_with_ec": sum(bool(row["ec_numbers"]) for row in rows),
        "proteins_with_cazy": sum(bool(row["cazy_families"]) for row in rows),
        "proteins_with_pfam": sum(bool(row["pfam_accessions"]) for row in rows),
        "proteins_with_interpro": sum(bool(row["interpro_accessions"]) for row in rows),
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    args = parse_args()
    identifiers, lengths = parse_fasta(args.proteins)
    protein_ids = set(identifiers)
    features = parse_gff(args.gff)
    eggnog = parse_eggnog(args.eggnog)
    interpro = parse_interpro(args.interpro)
    validate_annotation_ids(protein_ids, "eggNOG-mapper", eggnog)
    validate_annotation_ids(protein_ids, "InterProScan", interpro)
    rows = build_rows(
        args.mag_id,
        identifiers,
        lengths,
        features,
        eggnog,
        interpro,
    )
    write_table(args.output, rows)
    write_summary(args.summary, args.mag_id, rows)
    print(f"Integrated functional annotations for {len(rows)} proteins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
