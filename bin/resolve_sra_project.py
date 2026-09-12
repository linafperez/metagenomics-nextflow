#!/usr/bin/env python3
"""Resolve and freeze an NCBI BioProject RunInfo manifest.

The normal mode performs one E-utilities history search followed by one RunInfo
fetch.  ``--runinfo-file`` provides a completely offline path for testing and
reproducible re-resolution.  No sequence data are downloaded by this program.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA_VERSION = 2
PROJECT_PATTERN = re.compile(r"^PRJ(?:NA|EB|DB)[0-9]+$")
RUN_PATTERN = re.compile(r"^[SED]RR[0-9]+$")
EXPERIMENT_PATTERN = re.compile(r"^[SED]RX[0-9]+$")
BIOSAMPLE_PATTERN = re.compile(r"^SAM(?:N|EA|D)[0-9]+$")
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

RUNINFO_FILENAME = "sra_project_runinfo.csv"
RUN_MANIFEST_FILENAME = "sra_project_manifest.tsv"
SAMPLE_MANIFEST_FILENAME = "sra_sample_manifest.tsv"
EXCLUSIONS_FILENAME = "sra_project_exclusions.tsv"
SUMMARY_FILENAME = "sra_project_summary.json"
REQUESTED_FILENAME = "sra_requested_samples.tsv"
SAMPLE_METADATA_FILENAME = "sample_metadata.tsv"

RUN_FIELDS = (
    "project_accession",
    "sample_order",
    "sample_id",
    "identity_source",
    "biosample_accession",
    "group",
    "selection_file_sha256",
    "experiment_accession",
    "run_order",
    "run_accession",
    "layout",
    "strategy",
    "source",
    "platform",
    "model",
    "public_access",
    "spots",
    "spots_with_mates",
    "bases",
    "size_mb",
    "download_path",
    "eligibility",
    "exclusion_reason",
    "metadata_warnings",
)

SAMPLE_FIELDS = (
    "project_accession",
    "sample_order",
    "sample_id",
    "identity_source",
    "biosample_accession",
    "group",
    "selection_file_sha256",
    "experiment_accession",
    "run_count",
    "run_accessions",
    "layout",
    "strategy",
    "source",
    "platform",
    "model",
    "total_spots",
    "total_spots_with_mates",
    "total_bases",
    "total_size_mb",
    "metadata_warnings",
)

REQUESTED_FIELDS = (
    "selection_order",
    "biosample_accession",
    "group",
)

SAMPLE_METADATA_FIELDS = (
    "sample_id",
    "biosample_accession",
    "group",
    "run_accessions",
    "project_accession",
    "sample_order",
)

ALIASES: dict[str, tuple[str, ...]] = {
    "run": ("Run", "run_accession"),
    "project": ("BioProject", "project_accession"),
    "biosample": ("BioSample", "biosample_accession"),
    "experiment": ("Experiment", "experiment_accession"),
    "layout": ("LibraryLayout", "layout"),
    "strategy": ("LibraryStrategy", "strategy"),
    "source": ("LibrarySource", "source"),
    "platform": ("Platform", "platform"),
    "model": ("Model", "model"),
    "consent": ("Consent", "consent"),
    "access": ("public_access", "Access", "access"),
    "dbgap": ("dbgap_study_accession", "dbGaP"),
    "spots": ("spots",),
    "spots_with_mates": ("spots_with_mates", "spotsWithMates"),
    "bases": ("bases",),
    "size_mb": ("size_MB", "size_mb"),
    "download_path": ("download_path", "DownloadPath"),
}


class ResolverError(RuntimeError):
    """A user-facing resolution or validation error."""


@dataclass
class RunRow:
    """A normalized RunInfo row and its validation findings."""

    source_index: int
    values: dict[str, str]
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def add_reason(self, reason: str) -> None:
        if reason not in self.reasons:
            self.reasons.append(reason)

    def add_warning(self, warning: str) -> None:
        if warning not in self.warnings:
            self.warnings.append(warning)

    def as_dict(self) -> dict[str, str]:
        result = {name: self.values.get(name, "") for name in RUN_FIELDS}
        result["eligibility"] = "excluded" if self.reasons else "eligible"
        result["exclusion_reason"] = ";".join(self.reasons)
        result["metadata_warnings"] = ";".join(self.warnings)
        return result


@dataclass
class Resolution:
    project: str
    rows: list[RunRow]
    errors: list[str]
    allowed_platforms: tuple[str, ...]
    source: str
    requested: list[dict[str, str]] = field(default_factory=list)
    selection_file_sha256: str = ""
    selection_file_name: str = ""
    selection_file_bytes: int = 0
    group_column: str = ""
    project_runinfo_record_count: int = 0
    source_details: dict[str, Any] = field(default_factory=dict)

    @property
    def eligible(self) -> list[RunRow]:
        return [row for row in self.rows if not row.reasons]

    @property
    def excluded(self) -> list[RunRow]:
        return [row for row in self.rows if row.reasons]

    @property
    def valid(self) -> bool:
        return bool(self.requested) and bool(self.eligible) and not self.errors


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve an NCBI BioProject into deterministic run- and sample-level "
            "manifests for paired short-read shotgun metagenomics."
        )
    )
    parser.add_argument("project", nargs="?", help="BioProject accession")
    parser.add_argument("--project", dest="project_option", help="BioProject accession")
    parser.add_argument("--output-dir", type=Path, help="Directory for frozen reports")
    parser.add_argument(
        "--runinfo-file",
        type=Path,
        help="Offline RunInfo CSV/TSV fixture; disables all NCBI requests",
    )
    parser.add_argument(
        "--selection-file",
        type=Path,
        help="TSV containing the explicit biosample_accession selection",
    )
    parser.add_argument(
        "--group-column",
        default="",
        help="Optional selection column to preserve as categorical group metadata",
    )
    parser.add_argument(
        "--validate-existing",
        type=Path,
        metavar="REPORT_DIR",
        help="Validate an already frozen report directory without network access",
    )
    parser.add_argument(
        "--platform",
        "--platforms",
        "--allowed-platforms",
        dest="platforms",
        default="ILLUMINA,BGISEQ",
        help="Comma-separated platform allowlist (default: ILLUMINA,BGISEQ)",
    )
    parser.add_argument("--email", help="Contact email sent to NCBI E-utilities")
    parser.add_argument("--api-key", help="NCBI API key")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument(
        "--write-invalid-and-succeed",
        action="store_true",
        help=(
            "Write invalid reports but exit zero so an orchestrator can publish "
            "them before a separate validation process fails"
        ),
    )
    args = parser.parse_args(argv)

    if args.validate_existing is not None:
        if args.project or args.project_option or args.runinfo_file or args.output_dir:
            parser.error(
                "--validate-existing cannot be combined with project, --runinfo-file, "
                "or --output-dir"
            )
    else:
        selected = [value for value in (args.project, args.project_option) if value]
        if len(selected) != 1:
            parser.error("specify exactly one BioProject accession (positional or --project)")
        args.project = selected[0]
        if args.output_dir is None:
            parser.error("--output-dir is required")
        if args.selection_file is None:
            parser.error("--selection-file is required; BioProject-only discovery is disabled")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.retries < 0:
        parser.error("--retries cannot be negative")
    if args.group_column and not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", args.group_column):
        parser.error("--group-column contains unsupported characters")
    return args


def normalize_project(value: str) -> str:
    project = value.strip().upper()
    if not PROJECT_PATTERN.fullmatch(project):
        raise ResolverError(
            f"invalid BioProject accession {value!r}; expected PRJNA, PRJEB, or PRJDB "
            "followed by digits"
        )
    return project


def parse_platforms(value: str) -> tuple[str, ...]:
    platforms = tuple(
        sorted({part.strip().upper() for part in value.split(",") if part.strip()})
    )
    if not platforms:
        raise ResolverError("the platform allowlist cannot be empty")
    return platforms


def load_selection(
    path: Path, group_column: str
) -> tuple[list[dict[str, str]], bytes, str]:
    """Read and validate the explicit BioSample cohort without reordering it."""

    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ResolverError(f"cannot read BioSample selection file {path}: {exc}") from exc
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ResolverError("BioSample selection file is not valid UTF-8") from exc

    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    if not reader.fieldnames:
        raise ResolverError("BioSample selection file has no header")
    if any(not (name or "").strip() for name in reader.fieldnames):
        raise ResolverError("BioSample selection file contains an empty column name")
    if len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ResolverError("BioSample selection file contains duplicate column names")
    if "biosample_accession" not in reader.fieldnames:
        raise ResolverError(
            "BioSample selection file must contain a biosample_accession column"
        )
    if group_column and group_column not in reader.fieldnames:
        raise ResolverError(
            f"configured group column {group_column!r} is absent from the BioSample selection file"
        )
    if group_column == "biosample_accession":
        raise ResolverError("biosample_accession cannot also be the group column")

    requested: list[dict[str, str]] = []
    seen: set[str] = set()
    errors: list[str] = []
    try:
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                errors.append(f"selection line {line_number} contains unexpected extra fields")
                continue
            accession = (row.get("biosample_accession") or "").strip().upper()
            group = (row.get(group_column) or "").strip() if group_column else ""
            if not BIOSAMPLE_PATTERN.fullmatch(accession):
                errors.append(
                    f"selection line {line_number} has invalid BioSample accession {accession!r}"
                )
                continue
            if accession in seen:
                errors.append(f"selection line {line_number} duplicates BioSample {accession}")
                continue
            seen.add(accession)
            if group_column and not group:
                errors.append(
                    f"selection line {line_number} has an empty {group_column!r} group value"
                )
            if any(character in group for character in "\t\r\n"):
                errors.append(
                    f"selection line {line_number} group contains a tab or line break"
                )
            requested.append(
                {
                    "selection_order": str(len(requested) + 1),
                    "biosample_accession": accession,
                    "group": group,
                }
            )
    except csv.Error as exc:
        errors.append(f"BioSample selection TSV parsing failed: {exc}")
    if not requested:
        errors.append("BioSample selection file contains no samples")
    if errors:
        raise ResolverError("\n".join(errors))
    return requested, raw, hashlib.sha256(raw).hexdigest()


def _http_bytes(
    url: str,
    params: dict[str, str],
    *,
    timeout: float,
    retries: int,
    user_agent: str,
) -> bytes:
    encoded = urllib.parse.urlencode(params).encode("ascii")
    retryable_statuses = {429, 500, 502, 503, 504}
    for attempt in range(retries + 1):
        request = urllib.request.Request(
            url,
            data=encoded,
            headers={
                "Accept": "application/json,text/csv,text/plain,*/*",
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": user_agent,
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code not in retryable_statuses or attempt == retries:
                raise ResolverError(f"NCBI request failed with HTTP {exc.code}: {url}") from exc
        except urllib.error.URLError as exc:
            if attempt == retries:
                raise ResolverError(f"NCBI request failed: {exc.reason}") from exc
        time.sleep(min(2**attempt, 8))
    raise AssertionError("unreachable")


def fetch_runinfo(
    project: str,
    *,
    email: str | None,
    api_key: str | None,
    timeout: float,
    retries: int,
) -> tuple[bytes, dict[str, Any]]:
    """Resolve project membership once via ESearch history and EFetch RunInfo."""

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    user_agent = "metagenomics-nextflow-sra-resolver/1.0"
    if email:
        user_agent += f" ({email})"
    common = {"tool": "metagenomics_nextflow"}
    if email:
        common["email"] = email
    if api_key:
        common["api_key"] = api_key

    search_params = {
        **common,
        "db": "sra",
        "term": f"{project}[BioProject]",
        "retmode": "json",
        "retmax": "0",
        "usehistory": "y",
    }
    search_raw = _http_bytes(
        f"{base}/esearch.fcgi",
        search_params,
        timeout=timeout,
        retries=retries,
        user_agent=user_agent,
    )
    try:
        result = json.loads(search_raw.decode("utf-8"))["esearchresult"]
        count = int(result["count"])
        webenv = str(result["webenv"])
        query_key = str(result["querykey"])
    except (KeyError, TypeError, ValueError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ResolverError("NCBI ESearch returned an unexpected response") from exc

    if count == 0:
        return b"", {
            "eutils_database": "sra",
            "query": f"{project}[BioProject]",
            "esearch_record_count": 0,
        }

    fetch_params = {
        **common,
        "db": "sra",
        "query_key": query_key,
        "WebEnv": webenv,
        "rettype": "runinfo",
        "retmode": "text",
        "retstart": "0",
        "retmax": str(count),
    }
    runinfo = _http_bytes(
        f"{base}/efetch.fcgi",
        fetch_params,
        timeout=timeout,
        retries=retries,
        user_agent=user_agent,
    )
    return runinfo, {
        "eutils_database": "sra",
        "query": f"{project}[BioProject]",
        # SRA database records are experiment packages and can contain more
        # than one run, so this count must not be compared with RunInfo rows.
        "esearch_record_count": count,
    }


def _normalized_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.casefold())


def parse_runinfo(raw: bytes) -> tuple[list[dict[str, str]], list[str]]:
    errors: list[str] = []
    if not raw.strip():
        return [], ["RunInfo response is empty"]
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        return [], [f"RunInfo is not valid UTF-8: {exc}"]

    first_line = text.splitlines()[0] if text.splitlines() else ""
    delimiter = "\t" if first_line.count("\t") > first_line.count(",") else ","
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter=delimiter)
    if not reader.fieldnames:
        return [], ["RunInfo has no header"]
    normalized_headers = [_normalized_header(name or "") for name in reader.fieldnames]
    if len(set(normalized_headers)) != len(normalized_headers):
        errors.append("RunInfo contains duplicate column names after normalization")
    rows: list[dict[str, str]] = []
    try:
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                errors.append(f"RunInfo row {row_number} contains unexpected extra fields")
            clean = {(key or ""): (value or "").strip() for key, value in row.items() if key}
            if any(clean.values()):
                rows.append(clean)
    except csv.Error as exc:
        errors.append(f"RunInfo CSV parsing failed: {exc}")
    if not rows:
        errors.append("RunInfo contains no records")
    return rows, errors


def _lookup(row: dict[str, str], field: str) -> str:
    normalized = {_normalized_header(key): value.strip() for key, value in row.items()}
    for alias in ALIASES[field]:
        value = normalized.get(_normalized_header(alias), "")
        if value:
            return value.strip()
    return ""


def _integer_value(raw: str) -> tuple[str, int | None, str | None]:
    if not raw:
        return "", None, None
    try:
        value = int(raw)
    except ValueError:
        return raw, None, "not_an_integer"
    return str(value), value, None


def _decimal_value(raw: str) -> tuple[str, Decimal | None, str | None]:
    if not raw:
        return "", None, None
    try:
        value = Decimal(raw)
    except InvalidOperation:
        return raw, None, "not_a_number"
    if not value.is_finite():
        return raw, None, "not_a_finite_number"
    normalized = format(value, "f")
    if "." in normalized:
        normalized = normalized.rstrip("0").rstrip(".")
    return normalized or "0", value, None


def _public_access(row: dict[str, str], download_path: str) -> tuple[bool, str | None]:
    explicit = _lookup(row, "access").strip().casefold()
    consent = _lookup(row, "consent").strip().casefold()
    dbgap = _lookup(row, "dbgap").strip()
    false_values = {"false", "no", "0", "controlled", "restricted", "protected"}
    true_values = {"true", "yes", "1", "public", "open"}
    if explicit in false_values:
        return False, "explicit_non_public_access"
    if consent in false_values or dbgap:
        return False, "controlled_or_restricted_access"
    if not download_path:
        return False, "missing_public_download_path"
    if explicit and explicit not in true_values:
        return False, "unrecognized_access_metadata"
    return True, None


def _project_is_present(raw: str, expected: str) -> bool:
    accessions = {match.group(0) for match in re.finditer(r"PRJ(?:NA|EB|DB)[0-9]+", raw.upper())}
    return expected in accessions


def normalize_rows(
    project: str,
    raw_rows: Iterable[dict[str, str]],
    allowed_platforms: tuple[str, ...],
    requested: Sequence[dict[str, str]],
    selection_file_sha256: str,
) -> tuple[list[RunRow], list[str]]:
    requested_by_accession = {
        row["biosample_accession"]: row for row in requested
    }
    rows: list[RunRow] = []
    for index, raw in enumerate(raw_rows, start=1):
        run = _lookup(raw, "run").upper()
        experiment = _lookup(raw, "experiment").upper()
        biosample = _lookup(raw, "biosample").upper()
        if biosample not in requested_by_accession:
            continue
        selection = requested_by_accession[biosample]
        layout = _lookup(raw, "layout").upper()
        strategy = _lookup(raw, "strategy").upper()
        source = _lookup(raw, "source").upper()
        platform = _lookup(raw, "platform").upper()
        model = _lookup(raw, "model")
        download_path = _lookup(raw, "download_path")

        row = RunRow(
            source_index=index,
            values={
                "project_accession": project,
                "sample_order": selection["selection_order"],
                "sample_id": biosample,
                "identity_source": "BioSample",
                "biosample_accession": biosample,
                "group": selection["group"],
                "selection_file_sha256": selection_file_sha256,
                "experiment_accession": experiment,
                "run_order": "",
                "run_accession": run,
                "layout": layout,
                "strategy": strategy,
                "source": source,
                "platform": platform,
                "model": model,
                "public_access": "false",
                "spots": "",
                "spots_with_mates": "",
                "bases": "",
                "size_mb": "",
                "download_path": download_path,
                "eligibility": "",
                "exclusion_reason": "",
                "metadata_warnings": "",
            },
        )

        if not run:
            row.add_reason("missing_run_accession")
        elif not RUN_PATTERN.fullmatch(run):
            row.add_reason("invalid_run_accession")
        if experiment and not EXPERIMENT_PATTERN.fullmatch(experiment):
            row.add_warning("invalid_experiment_accession")
        elif not experiment:
            row.add_warning("missing_experiment_accession")

        raw_project = _lookup(raw, "project")
        if not raw_project:
            row.add_reason("missing_bioproject_accession")
        elif not _project_is_present(raw_project, project):
            row.add_reason("bioproject_accession_mismatch")

        if layout != "PAIRED":
            row.add_reason("library_layout_not_paired" if layout else "missing_library_layout")
        if strategy != "WGS":
            row.add_reason("library_strategy_not_wgs" if strategy else "missing_library_strategy")
        if source != "METAGENOMIC":
            row.add_reason("library_source_not_metagenomic" if source else "missing_library_source")
        if not platform:
            row.add_reason("missing_platform")
        elif platform not in allowed_platforms:
            row.add_reason("unsupported_platform")

        is_public, access_reason = _public_access(raw, download_path)
        row.values["public_access"] = "true" if is_public else "false"
        if access_reason:
            row.add_reason(access_reason)

        spots_text, spots, spots_error = _integer_value(_lookup(raw, "spots"))
        row.values["spots"] = spots_text
        if spots_error:
            row.add_reason("invalid_spots")
        elif spots is None:
            row.add_warning("spots_unavailable")
        elif spots <= 0:
            row.add_reason("non_positive_spots")

        mates_text, mates, mates_error = _integer_value(_lookup(raw, "spots_with_mates"))
        row.values["spots_with_mates"] = mates_text
        if mates_error:
            row.add_reason("invalid_spots_with_mates")
        elif mates is None:
            row.add_warning("spots_with_mates_unavailable")
        elif mates < 0:
            row.add_reason("negative_spots_with_mates")
        elif spots is not None and mates > spots:
            row.add_reason("spots_with_mates_exceeds_spots")
        elif spots is not None and mates != spots:
            # RunInfo often reports fewer paired spots than total spots even
            # for PAIRED libraries. Treat that aggregate as advisory: the
            # acquisition helper still validates the materialized mates in
            # full and refuses unequal/non-matching FASTQ pairs.
            row.add_warning("spots_with_mates_differs_from_spots")

        bases_text, bases, bases_error = _integer_value(_lookup(raw, "bases"))
        row.values["bases"] = bases_text
        if bases_error:
            row.add_reason("invalid_bases")
        elif bases is None:
            row.add_warning("bases_unavailable")
        elif bases <= 0:
            row.add_reason("non_positive_bases")

        size_text, size_mb, size_error = _decimal_value(_lookup(raw, "size_mb"))
        row.values["size_mb"] = size_text
        if size_error:
            row.add_warning("invalid_size_mb")
        elif size_mb is None:
            row.add_warning("size_mb_unavailable")
        elif size_mb < 0:
            row.add_warning("negative_size_mb")

        if not model:
            row.add_warning("model_unavailable")
        rows.append(row)

    errors: list[str] = []
    run_counts = Counter(row.values["run_accession"] for row in rows if row.values["run_accession"])
    for row in rows:
        if run_counts[row.values["run_accession"]] > 1:
            row.add_reason("duplicate_run_accession")
    for run, count in sorted(run_counts.items()):
        if count > 1:
            errors.append(
                f"selected RunInfo metadata contains duplicate or contradictory rows for run {run}"
            )

    eligible = sorted(
        (row for row in rows if not row.reasons),
        key=lambda item: (
            int(item.values["sample_order"]),
            item.values["run_accession"],
        ),
    )
    run_orders: defaultdict[str, int] = defaultdict(int)
    for row in eligible:
        sample_id = row.values["sample_id"]
        run_orders[sample_id] += 1
        row.values["run_order"] = str(run_orders[sample_id])

    for selection in requested:
        accession = selection["biosample_accession"]
        matching = [row for row in rows if row.values["biosample_accession"] == accession]
        if not matching:
            errors.append(
                f"requested BioSample {accession} is absent from BioProject RunInfo metadata"
            )
            continue
        # Membership is represented directly by normalized RunInfo fields,
        # avoiding any inference from human-readable sample names.
        if all(
            "bioproject_accession_mismatch" in row.reasons
            or "missing_bioproject_accession" in row.reasons
            for row in matching
        ):
            errors.append(
                f"requested BioSample {accession} does not belong to BioProject {project}"
            )
        if not any(not row.reasons for row in matching):
            errors.append(
                f"requested BioSample {accession} has no eligible paired WGS METAGENOMIC run"
            )
    return rows, errors


def _unique_join(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def _sum_integers(rows: Sequence[dict[str, str]], field_name: str) -> str:
    values = [row[field_name] for row in rows]
    if not values or any(not value for value in values):
        return ""
    return str(sum(int(value) for value in values))


def _sum_decimals(rows: Sequence[dict[str, str]], field_name: str) -> str:
    values = [row[field_name] for row in rows]
    if not values or any(not value for value in values):
        return ""
    result = format(sum((Decimal(value) for value in values), Decimal(0)), "f")
    if "." in result:
        result = result.rstrip("0").rstrip(".")
    return result or "0"


def build_sample_rows(run_rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    grouped: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in run_rows:
        grouped[row["sample_id"]].append(row)
    result: list[dict[str, str]] = []
    for sample_id in sorted(grouped, key=lambda value: int(grouped[value][0]["sample_order"])):
        rows = sorted(grouped[sample_id], key=lambda row: int(row["run_order"]))
        warnings: set[str] = set()
        for row in rows:
            warnings.update(filter(None, row["metadata_warnings"].split(";")))
        result.append(
            {
                "project_accession": rows[0]["project_accession"],
                "sample_order": rows[0]["sample_order"],
                "sample_id": sample_id,
                "identity_source": _unique_join(row["identity_source"] for row in rows),
                "biosample_accession": _unique_join(row["biosample_accession"] for row in rows),
                "group": _unique_join(row["group"] for row in rows),
                "selection_file_sha256": _unique_join(
                    row["selection_file_sha256"] for row in rows
                ),
                "experiment_accession": _unique_join(row["experiment_accession"] for row in rows),
                "run_count": str(len(rows)),
                "run_accessions": ";".join(row["run_accession"] for row in rows),
                "layout": _unique_join(row["layout"] for row in rows),
                "strategy": _unique_join(row["strategy"] for row in rows),
                "source": _unique_join(row["source"] for row in rows),
                "platform": _unique_join(row["platform"] for row in rows),
                "model": _unique_join(row["model"] for row in rows),
                "total_spots": _sum_integers(rows, "spots"),
                "total_spots_with_mates": _sum_integers(rows, "spots_with_mates"),
                "total_bases": _sum_integers(rows, "bases"),
                "total_size_mb": _sum_decimals(rows, "size_mb"),
                "metadata_warnings": ";".join(sorted(warnings)),
            }
        )
    return result


def build_sample_metadata(
    sample_rows: Sequence[dict[str, str]],
) -> list[dict[str, str]]:
    return [
        {
            "sample_id": row["sample_id"],
            "biosample_accession": row["biosample_accession"],
            "group": row["group"],
            "run_accessions": row["run_accessions"],
            "project_accession": row["project_accession"],
            "sample_order": row["sample_order"],
        }
        for row in sample_rows
    ]


def _tsv_bytes(rows: Sequence[dict[str, str]], fields: Sequence[str]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=list(fields),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def render_reports(resolution: Resolution, raw_runinfo: bytes) -> dict[str, bytes]:
    run_rows = [
        row.as_dict()
        for row in sorted(
            resolution.eligible,
            key=lambda item: (int(item.values["sample_order"]), int(item.values["run_order"])),
        )
    ]
    exclusion_rows = [
        row.as_dict()
        for row in sorted(
            resolution.excluded,
            key=lambda item: (
                item.values["sample_id"],
                item.values["run_accession"],
                item.source_index,
            ),
        )
    ]
    sample_rows = build_sample_rows(run_rows)
    requested_rows = [dict(row) for row in resolution.requested]
    sample_metadata = build_sample_metadata(sample_rows)
    reports = {
        RUNINFO_FILENAME: raw_runinfo,
        RUN_MANIFEST_FILENAME: _tsv_bytes(run_rows, RUN_FIELDS),
        SAMPLE_MANIFEST_FILENAME: _tsv_bytes(sample_rows, SAMPLE_FIELDS),
        EXCLUSIONS_FILENAME: _tsv_bytes(exclusion_rows, RUN_FIELDS),
        REQUESTED_FILENAME: _tsv_bytes(requested_rows, REQUESTED_FIELDS),
        SAMPLE_METADATA_FILENAME: _tsv_bytes(
            sample_metadata, SAMPLE_METADATA_FIELDS
        ),
    }
    reason_counts = Counter(
        reason for row in resolution.excluded for reason in row.reasons
    )
    warning_counts = Counter(warning for row in resolution.rows for warning in row.warnings)
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "project_accession": resolution.project,
        "valid": resolution.valid,
        "resolution_source": resolution.source,
        "group_column": resolution.group_column,
        "allowed_platforms": list(resolution.allowed_platforms),
        "project_runinfo_record_count": resolution.project_runinfo_record_count,
        "record_count": len(resolution.rows),
        "requested_biosample_count": len(requested_rows),
        "resolved_biosample_count": len(sample_rows),
        "eligible_run_count": len(run_rows),
        "excluded_run_count": len(exclusion_rows),
        "sample_count": len(sample_rows),
        "exclusion_reason_counts": dict(sorted(reason_counts.items())),
        "metadata_warning_counts": dict(sorted(warning_counts.items())),
        "validation_errors": resolution.errors,
        "source_details": resolution.source_details,
        "selection_file": {
            "name": resolution.selection_file_name,
            "sha256": resolution.selection_file_sha256,
            "bytes": resolution.selection_file_bytes,
        },
        "files": {
            name: {"sha256": _sha256(content), "bytes": len(content)}
            for name, content in sorted(reports.items())
        },
    }
    reports[SUMMARY_FILENAME] = (
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")
    return reports


def write_reports(output_dir: Path, reports: dict[str, bytes]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename in (
        RUNINFO_FILENAME,
        RUN_MANIFEST_FILENAME,
        SAMPLE_MANIFEST_FILENAME,
        EXCLUSIONS_FILENAME,
        REQUESTED_FILENAME,
        SAMPLE_METADATA_FILENAME,
        SUMMARY_FILENAME,
    ):
        _atomic_write(output_dir / filename, reports[filename])


def _read_tsv(path: Path, expected_fields: Sequence[str]) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(expected_fields):
                raise ResolverError(
                    f"{path.name} has an unexpected header; expected "
                    + ",".join(expected_fields)
                )
            return [dict(row) for row in reader]
    except OSError as exc:
        raise ResolverError(f"cannot read {path}: {exc}") from exc


def validate_existing(
    report_dir: Path,
    selection_file: Path | None = None,
    group_column: str = "",
) -> list[str]:
    errors: list[str] = []
    try:
        summary_raw = (report_dir / SUMMARY_FILENAME).read_bytes()
        summary = json.loads(summary_raw.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return [f"cannot read valid {SUMMARY_FILENAME}: {exc}"]

    if summary.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"unsupported summary schema_version: {summary.get('schema_version')!r}")
    try:
        project = normalize_project(str(summary.get("project_accession", "")))
    except ResolverError as exc:
        errors.append(str(exc))
        project = ""
    if summary.get("valid") is not True:
        errors.append("frozen project summary is marked invalid")

    report_names = (
        RUNINFO_FILENAME,
        RUN_MANIFEST_FILENAME,
        SAMPLE_MANIFEST_FILENAME,
        EXCLUSIONS_FILENAME,
        REQUESTED_FILENAME,
        SAMPLE_METADATA_FILENAME,
    )
    for filename in report_names:
        path = report_dir / filename
        try:
            content = path.read_bytes()
        except OSError as exc:
            errors.append(f"cannot read {filename}: {exc}")
            continue
        expected = summary.get("files", {}).get(filename, {})
        if expected.get("sha256") != _sha256(content):
            errors.append(f"SHA-256 mismatch for {filename}")
        if expected.get("bytes") != len(content):
            errors.append(f"byte-count mismatch for {filename}")

    try:
        run_rows = _read_tsv(report_dir / RUN_MANIFEST_FILENAME, RUN_FIELDS)
        exclusion_rows = _read_tsv(report_dir / EXCLUSIONS_FILENAME, RUN_FIELDS)
        sample_rows = _read_tsv(report_dir / SAMPLE_MANIFEST_FILENAME, SAMPLE_FIELDS)
        requested_rows = _read_tsv(report_dir / REQUESTED_FILENAME, REQUESTED_FIELDS)
        sample_metadata = _read_tsv(
            report_dir / SAMPLE_METADATA_FILENAME, SAMPLE_METADATA_FIELDS
        )
    except ResolverError as exc:
        errors.append(str(exc))
        return errors

    if not requested_rows:
        errors.append("frozen requested-BioSample manifest contains no samples")
    if not run_rows:
        errors.append("frozen run manifest contains no eligible runs")
    selection_description = summary.get("selection_file", {})
    selection_hash = selection_description.get("sha256", "")
    if not re.fullmatch(r"[0-9a-f]{64}", str(selection_hash)):
        errors.append("frozen summary contains an invalid selection-file SHA-256")
    if summary.get("group_column", "") != group_column and selection_file is not None:
        errors.append("frozen group column differs from the requested group column")

    requested_by_accession: dict[str, dict[str, str]] = {}
    for expected_order, row in enumerate(requested_rows, start=1):
        accession = row["biosample_accession"]
        if row["selection_order"] != str(expected_order):
            errors.append("requested BioSamples are not in consecutive selection order")
        if not BIOSAMPLE_PATTERN.fullmatch(accession):
            errors.append(f"requested manifest contains invalid BioSample {accession!r}")
        if accession in requested_by_accession:
            errors.append(f"requested manifest duplicates BioSample {accession}")
        requested_by_accession[accession] = row
        if summary.get("group_column") and not row["group"]:
            errors.append(f"requested BioSample {accession} has an empty group")

    if selection_file is not None:
        try:
            current_requested, current_raw, current_hash = load_selection(
                selection_file, group_column
            )
        except ResolverError as exc:
            errors.append(str(exc))
        else:
            if current_requested != requested_rows:
                errors.append("frozen requested BioSamples differ from the selection file")
            if current_hash != selection_hash:
                errors.append("selection-file SHA-256 differs from frozen SRA state")
            if len(current_raw) != selection_description.get("bytes"):
                errors.append("selection-file byte count differs from frozen SRA state")
    seen_runs: set[str] = set()
    expected_sort: list[tuple[int, int]] = []
    orders_by_sample: defaultdict[str, list[int]] = defaultdict(list)
    for row_number, row in enumerate(run_rows, start=2):
        label = f"{RUN_MANIFEST_FILENAME} row {row_number}"
        run = row["run_accession"]
        if project and row["project_accession"] != project:
            errors.append(f"{label}: project accession mismatch")
        if not RUN_PATTERN.fullmatch(run):
            errors.append(f"{label}: invalid run accession")
        elif run in seen_runs:
            errors.append(f"{label}: duplicate run accession")
        seen_runs.add(run)
        if not SAMPLE_ID_PATTERN.fullmatch(row["sample_id"]):
            errors.append(f"{label}: invalid sample_id")
        selection = requested_by_accession.get(row["biosample_accession"])
        if selection is None or row["sample_id"] != row["biosample_accession"]:
            errors.append(f"{label}: run is outside the explicit BioSample selection")
        else:
            if row["sample_order"] != selection["selection_order"]:
                errors.append(f"{label}: sample_order differs from selection order")
            if row["group"] != selection["group"]:
                errors.append(f"{label}: group differs from explicit selection")
        if row["selection_file_sha256"] != selection_hash:
            errors.append(f"{label}: selection-file hash mismatch")
        for name, expected in (
            ("layout", "PAIRED"),
            ("strategy", "WGS"),
            ("source", "METAGENOMIC"),
            ("public_access", "true"),
            ("eligibility", "eligible"),
        ):
            if row[name] != expected:
                errors.append(f"{label}: {name} must be {expected}")
        if row["exclusion_reason"]:
            errors.append(f"{label}: eligible row has an exclusion reason")
        try:
            sample_order = int(row["sample_order"])
            run_order = int(row["run_order"])
            if sample_order <= 0 or run_order <= 0:
                raise ValueError
            expected_sort.append((sample_order, run_order))
            orders_by_sample[row["sample_id"]].append(run_order)
        except ValueError:
            errors.append(f"{label}: sample_order and run_order must be positive integers")
    if expected_sort != sorted(expected_sort):
        errors.append("run manifest is not in deterministic sample/run order")
    sample_order_values: dict[str, str] = {}
    for row in run_rows:
        previous = sample_order_values.setdefault(row["sample_id"], row["sample_order"])
        if previous != row["sample_order"]:
            errors.append(f"sample {row['sample_id']} has inconsistent sample_order")
    try:
        numeric_sample_orders = sorted(int(value) for value in sample_order_values.values())
        if numeric_sample_orders != list(range(1, len(numeric_sample_orders) + 1)):
            errors.append("sample_order values are not consecutive")
    except ValueError:
        pass
    for sample_id, orders in orders_by_sample.items():
        if orders != list(range(1, len(orders) + 1)):
            errors.append(f"run_order values are not consecutive for {sample_id}")

    rebuilt_samples = build_sample_rows(run_rows) if run_rows else []
    if sample_rows != rebuilt_samples:
        errors.append("sample manifest does not match the frozen run manifest")
    if sample_metadata != build_sample_metadata(sample_rows):
        errors.append("sample metadata does not match the frozen sample manifest")
    if {row["biosample_accession"] for row in sample_rows} != set(requested_by_accession):
        errors.append("not every requested BioSample has an eligible resolved sample")
    for row_number, row in enumerate(exclusion_rows, start=2):
        label = f"{EXCLUSIONS_FILENAME} row {row_number}"
        selection = requested_by_accession.get(row["biosample_accession"])
        if selection is None:
            errors.append(f"{label}: excluded run is outside the explicit selection")
        elif row["group"] != selection["group"]:
            errors.append(f"{label}: group differs from explicit selection")
        if row["selection_file_sha256"] != selection_hash:
            errors.append(f"{label}: selection-file hash mismatch")
        if row["eligibility"] != "excluded" or not row["exclusion_reason"]:
            errors.append(f"{label}: excluded row lacks an exclusion reason")
    expected_counts = {
        "record_count": len(run_rows) + len(exclusion_rows),
        "requested_biosample_count": len(requested_rows),
        "resolved_biosample_count": len(sample_rows),
        "eligible_run_count": len(run_rows),
        "excluded_run_count": len(exclusion_rows),
        "sample_count": len(sample_rows),
    }
    for key, value in expected_counts.items():
        if summary.get(key) != value:
            errors.append(f"summary {key} does not match report contents")
    return errors


def run_resolution(args: argparse.Namespace) -> int:
    requested: list[dict[str, str]] = []
    selection_raw = b""
    selection_hash = ""
    initial_errors: list[str] = []
    try:
        requested, selection_raw, selection_hash = load_selection(
            args.selection_file, args.group_column
        )
    except ResolverError as exc:
        initial_errors.append(str(exc))
        try:
            selection_raw = args.selection_file.read_bytes()
            selection_hash = _sha256(selection_raw)
        except OSError:
            pass
    try:
        project = normalize_project(args.project)
        platforms = parse_platforms(args.platforms)
    except ResolverError as exc:
        initial_errors.append(str(exc))
        project = (args.project or "").strip().upper()
        platforms = tuple()

    if initial_errors:
        resolution = Resolution(
            project=project,
            rows=[],
            errors=initial_errors,
            allowed_platforms=platforms,
            source="offline_runinfo" if args.runinfo_file else "ncbi_eutils",
            requested=requested,
            selection_file_sha256=selection_hash,
            selection_file_name=args.selection_file.name,
            selection_file_bytes=len(selection_raw),
            group_column=args.group_column,
        )
        raw = args.runinfo_file.read_bytes() if args.runinfo_file and args.runinfo_file.is_file() else b""
        write_reports(args.output_dir, render_reports(resolution, raw))
        for error in initial_errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 0 if args.write_invalid_and_succeed else 2

    source_details: dict[str, Any] = {}
    try:
        if args.runinfo_file:
            raw = args.runinfo_file.read_bytes()
            source = "offline_runinfo"
            source_details = {"fixture_name": args.runinfo_file.name}
        else:
            raw, source_details = fetch_runinfo(
                project,
                email=args.email,
                api_key=args.api_key,
                timeout=args.timeout,
                retries=args.retries,
            )
            source = "ncbi_eutils"
    except (OSError, ResolverError) as exc:
        raw = b""
        source = "offline_runinfo" if args.runinfo_file else "ncbi_eutils"
        resolution = Resolution(
            project=project,
            rows=[],
            errors=[str(exc)],
            allowed_platforms=platforms,
            source=source,
            requested=requested,
            selection_file_sha256=selection_hash,
            selection_file_name=args.selection_file.name,
            selection_file_bytes=len(selection_raw),
            group_column=args.group_column,
            source_details=source_details,
        )
        write_reports(args.output_dir, render_reports(resolution, raw))
        print(f"ERROR: {exc}", file=sys.stderr)
        return 0 if args.write_invalid_and_succeed else 2

    raw_rows, parse_errors = parse_runinfo(raw)
    rows, selection_errors = normalize_rows(
        project,
        raw_rows,
        platforms,
        requested,
        selection_hash,
    )
    resolution = Resolution(
        project=project,
        rows=rows,
        errors=[*parse_errors, *selection_errors],
        allowed_platforms=platforms,
        source=source,
        requested=requested,
        selection_file_sha256=selection_hash,
        selection_file_name=args.selection_file.name,
        selection_file_bytes=len(selection_raw),
        group_column=args.group_column,
        project_runinfo_record_count=len(raw_rows),
        source_details=source_details,
    )
    reports = render_reports(resolution, raw)
    write_reports(args.output_dir, reports)
    summary = json.loads(reports[SUMMARY_FILENAME])
    print(
        f"Resolved {summary['record_count']} RunInfo records into "
        f"{summary['sample_count']} samples and {summary['eligible_run_count']} eligible runs; "
        f"excluded {summary['excluded_run_count']} runs.",
        file=sys.stderr,
    )
    if resolution.valid:
        return 0
    for error in resolution.errors:
        print(f"ERROR: {error}", file=sys.stderr)
    print(
        f"ERROR: BioProject {project} does not satisfy the pipeline input contract; "
        f"see {args.output_dir / EXCLUSIONS_FILENAME} and {args.output_dir / SUMMARY_FILENAME}",
        file=sys.stderr,
    )
    return 0 if args.write_invalid_and_succeed else 2


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.validate_existing is not None:
        errors = validate_existing(
            args.validate_existing,
            args.selection_file,
            args.group_column,
        )
        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            return 0 if args.write_invalid_and_succeed else 2
        print(f"Validated frozen SRA reports in {args.validate_existing}", file=sys.stderr)
        return 0
    return run_resolution(args)


if __name__ == "__main__":
    raise SystemExit(main())
