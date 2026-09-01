#!/usr/bin/env python3
"""Create, validate, reconcile, seal results, and clean durable SRA checkpoints.

The checkpoint record is the commit marker: it is written atomically only after
both compressed FASTQ mates and the frozen-manifest association have been fully
validated.  Directory existence alone is never treated as completion.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import stat
import sys
import tempfile
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Iterator, Sequence


SCHEMA_VERSION = 1
CHECKPOINT_OWNER_NAME = "sra_checkpoint_owner.json"
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
RUN_FIELDS = (
    "project_accession",
    "sample_order",
    "sample_id",
    "identity_source",
    "biosample_accession",
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
CHECKPOINT_FIELDS = (
    "schema_version",
    "project_accession",
    "sample_order",
    "sample_id",
    "identity_source",
    "biosample_accession",
    "experiment_accessions",
    "run_count",
    "run_accessions",
    "run_manifest_sha256",
    "read_1",
    "read_1_bytes",
    "read_1_sha256",
    "read_2",
    "read_2_bytes",
    "read_2_sha256",
    "paired_fastq_records",
    "reports_json",
    "completed_at_utc",
    "status",
)
PENDING_FIELDS = ("sample_order", "sample_id", "reason")
SCIENTIFIC_RESULT_ROOTS = (
    "01_quality_control_and_filtering",
    "02_mag_construction",
    "03_taxonomic_classification_and_phylogenomics",
    "04_gene_prediction_and_functional_annotation",
    "05_mag_abundance_estimation",
    "06_global_processing_evaluation",
)
REQUIRED_SCIENTIFIC_ARTIFACTS = {
    "final_mags": (
        "02_mag_construction/final_catalog/final_catalog/*.fa",
    ),
    "final_catalog_provenance": (
        "02_mag_construction/final_catalog/final_catalog.provenance.tsv",
    ),
    "final_catalog_quality": (
        "02_mag_construction/final_catalog/final_catalog.quality.tsv",
    ),
    "final_checkm2": (
        "02_mag_construction/final_catalog/evaluation/checkm2/*.checkm2.quality_report.tsv",
    ),
    "final_gunc": (
        "02_mag_construction/final_catalog/evaluation/gunc/*.gunc.summary.tsv",
    ),
    "species_representatives": (
        "02_mag_construction/final_catalog/species_95/*.representatives/*.fa",
    ),
    "gtdbtk_summaries": (
        "03_taxonomic_classification_and_phylogenomics/gtdbtk/*.gtdbtk.bac120.summary.tsv",
        "03_taxonomic_classification_and_phylogenomics/gtdbtk/*.gtdbtk.ar53.summary.tsv",
    ),
    "phylogenomic_tree": (
        "03_taxonomic_classification_and_phylogenomics/phylophlan/*.phylophlan.tree.nwk",
    ),
    "functional_annotations": (
        "04_gene_prediction_and_functional_annotation/integrated/*/*.functional_annotations.tsv",
    ),
    "mag_abundance": (
        "05_mag_abundance_estimation/final_catalog.mag_abundance.long.tsv",
    ),
    "multiqc_report": (
        "06_global_processing_evaluation/global_processing_evaluation.multiqc.html",
    ),
    "software_versions": (
        "pipeline_info/software_versions.tsv",
    ),
}
BASELINE_OUTPUT_PATHS = {
    "multiqc_report": (
        "06_global_processing_evaluation/global_processing_evaluation.multiqc.html"
    ),
    "software_versions": "pipeline_info/software_versions.tsv",
    "mag_abundance": (
        "05_mag_abundance_estimation/final_catalog.mag_abundance.long.tsv"
    ),
}
ABUNDANCE_FIELDS = (
    "sample",
    "mag_id",
    "relative_abundance_percent",
    "mean_coverage",
    "covered_fraction",
    "genome_length",
)


class CheckpointError(RuntimeError):
    """A user-facing validation or checkpoint error."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write(path: Path, content: bytes) -> None:
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


def atomic_json(path: Path, value: Any) -> None:
    atomic_write(
        path,
        (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode(),
    )


def tsv_bytes(rows: Iterable[dict[str, Any]], fields: Sequence[str]) -> bytes:
    import io

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
    return stream.getvalue().encode()


def read_run_manifest(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(RUN_FIELDS):
                raise CheckpointError("frozen SRA run manifest has an unexpected header")
            rows = [dict(row) for row in reader]
    except (OSError, csv.Error) as exc:
        raise CheckpointError(f"cannot read frozen SRA run manifest {path}: {exc}") from exc
    if not rows:
        raise CheckpointError("frozen SRA run manifest contains no eligible runs")
    try:
        expected = sorted(
            rows,
            key=lambda row: (int(row["sample_order"]), int(row["run_order"])),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise CheckpointError(
            "frozen SRA run manifest contains an invalid sample_order or run_order"
        ) from exc
    if rows != expected:
        raise CheckpointError("frozen SRA run manifest is not deterministically ordered")
    for row in rows:
        if row["eligibility"] != "eligible" or row["exclusion_reason"]:
            raise CheckpointError(f"run {row['run_accession']} is not eligible")
    projects = {row["project_accession"] for row in rows}
    if len(projects) != 1:
        raise CheckpointError("frozen SRA run manifest contains multiple BioProjects")
    return rows


def expected_checkpoint_owner(
    rows: Sequence[dict[str, str]], manifest_hash: str
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "project_accession": rows[0]["project_accession"],
        "run_manifest_sha256": manifest_hash,
    }


def validate_checkpoint_owner(
    root: Path, project_accession: str, manifest_hash: str
) -> dict[str, Any]:
    owner_path = safe_child(root / CHECKPOINT_OWNER_NAME, root)
    if owner_path.is_symlink():
        raise CheckpointError("checkpoint ownership record must not be a symlink")
    try:
        with owner_path.open("r", encoding="utf-8") as handle:
            owner = json.load(handle)
    except FileNotFoundError as exc:
        raise CheckpointError(
            "checkpoint root has no ownership record; use a dedicated empty root"
        ) from exc
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckpointError(f"cannot validate checkpoint ownership record: {exc}") from exc
    if not isinstance(owner, dict) or owner.get("schema_version") != SCHEMA_VERSION:
        raise CheckpointError("checkpoint ownership record has an unsupported schema")
    if owner.get("project_accession") != project_accession:
        raise CheckpointError(
            "checkpoint root is owned by BioProject "
            f"{owner.get('project_accession')!r}, not {project_accession!r}"
        )
    if owner.get("run_manifest_sha256") != manifest_hash:
        raise CheckpointError(
            "checkpoint root is bound to a different frozen manifest; "
            "use a dedicated empty root for the new cohort"
        )
    return owner


def claim_checkpoint_root(
    root: Path, rows: Sequence[dict[str, str]], manifest_hash: str
) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=True)
    expected = expected_checkpoint_owner(rows, manifest_hash)
    owner_path = safe_child(root / CHECKPOINT_OWNER_NAME, root)
    if owner_path.exists() or owner_path.is_symlink():
        return validate_checkpoint_owner(
            root, expected["project_accession"], expected["run_manifest_sha256"]
        )
    try:
        existing_entries = list(root.iterdir())
    except OSError as exc:
        raise CheckpointError(f"cannot inspect checkpoint root {root}: {exc}") from exc
    if existing_entries:
        raise CheckpointError(
            "unclaimed checkpoint root is not empty; use a dedicated empty root "
            "to avoid overwriting unrelated data"
        )
    owner = {**expected, "created_at_utc": utc_now()}
    content = (json.dumps(owner, indent=2, sort_keys=True) + "\n").encode()
    try:
        descriptor = os.open(owner_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError:
        return validate_checkpoint_owner(
            root, expected["project_accession"], expected["run_manifest_sha256"]
        )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        try:
            owner_path.unlink()
        except FileNotFoundError:
            pass
        raise
    return owner


def require_active_checkpoint_lifecycle(root: Path) -> None:
    cleanup_record = safe_child(root / "sra_checkpoint_cleanup.json", root)
    if cleanup_record.exists() or cleanup_record.is_symlink():
        raise CheckpointError(
            "checkpoint cleanup has already started; use a dedicated empty root "
            "for a new processing lifecycle"
        )


def group_samples(rows: Sequence[dict[str, str]]) -> list[tuple[str, list[dict[str, str]]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["sample_id"], []).append(row)
    return sorted(
        grouped.items(),
        key=lambda item: (int(item[1][0]["sample_order"]), item[0]),
    )


def safe_child(path: Path, root: Path) -> Path:
    resolved = path.expanduser().resolve()
    resolved_root = root.expanduser().resolve()
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise CheckpointError(f"checkpoint path escapes configured root: {resolved}") from exc
    return resolved


def normalize_read_name(header: bytes, mate: int) -> bytes:
    if not header.startswith(b"@"):
        raise CheckpointError("FASTQ record header does not begin with '@'")
    fields = header[1:].strip().split()
    if not fields:
        raise CheckpointError("FASTQ record has an empty read name")
    name = fields[0]
    suffix = f"/{mate}".encode()
    other_suffix = f"/{3 - mate}".encode()
    if name.endswith(other_suffix):
        raise CheckpointError(f"mate {mate} FASTQ contains a mate-{3 - mate} read name")
    if name.endswith(suffix):
        name = name[: -len(suffix)]
    if len(fields) > 1 and fields[1][:2] in (b"1:", b"2:"):
        expected = f"{mate}:".encode()
        if fields[1][:2] != expected:
            raise CheckpointError(f"mate {mate} FASTQ has an inconsistent CASAVA mate token")
    return name


def fastq_records(path: Path, mate: int) -> Iterator[tuple[bytes, bytes, bytes, bytes]]:
    try:
        with gzip.open(path, "rb") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                if not sequence or not plus or not quality:
                    raise CheckpointError(f"truncated FASTQ record in {path}")
                if not plus.startswith(b"+"):
                    raise CheckpointError(f"FASTQ separator does not begin with '+' in {path}")
                if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
                    raise CheckpointError(f"sequence/quality length mismatch in {path}")
                normalize_read_name(header, mate)
                yield header, sequence, plus, quality
    except (OSError, EOFError) as exc:
        raise CheckpointError(f"invalid or truncated gzip FASTQ {path}: {exc}") from exc


def validate_pair(read_1: Path, read_2: Path) -> int:
    if not read_1.is_file() or not read_2.is_file():
        raise CheckpointError("both checkpoint FASTQ mates must exist")
    if read_1.stat().st_size <= 0 or read_2.stat().st_size <= 0:
        raise CheckpointError("both checkpoint FASTQ mates must be non-empty")
    iterator_1 = fastq_records(read_1, 1)
    iterator_2 = fastq_records(read_2, 2)
    count = 0
    while True:
        record_1 = next(iterator_1, None)
        record_2 = next(iterator_2, None)
        if record_1 is None or record_2 is None:
            if record_1 is not None or record_2 is not None:
                raise CheckpointError("checkpoint FASTQ mates contain different record counts")
            break
        if normalize_read_name(record_1[0], 1) != normalize_read_name(record_2[0], 2):
            raise CheckpointError(f"checkpoint FASTQ mate names differ at record {count + 1}")
        count += 1
    if count == 0:
        raise CheckpointError("checkpoint FASTQ pair contains no records")
    return count


def copy_atomic(source: Path, target: Path) -> Path:
    if not source.is_file() or source.stat().st_size <= 0:
        raise CheckpointError(f"cannot checkpoint missing or empty file: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.parent / f".{target.name}.partial.{os.getpid()}"
    try:
        with source.open("rb") as reader, temporary.open("xb") as writer:
            shutil.copyfileobj(reader, writer, length=1024 * 1024)
            writer.flush()
            os.fsync(writer.fileno())
        os.replace(temporary, target)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return target


def sample_record(
    sample_rows: Sequence[dict[str, str]],
    manifest_hash: str,
    read_1: Path,
    read_2: Path,
    record_count: int,
    reports: Sequence[Path],
) -> dict[str, Any]:
    first = sample_rows[0]
    experiments = sorted({row["experiment_accession"] for row in sample_rows if row["experiment_accession"]})
    return {
        "schema_version": SCHEMA_VERSION,
        "project_accession": first["project_accession"],
        "sample_order": int(first["sample_order"]),
        "sample_id": first["sample_id"],
        "identity_source": first["identity_source"],
        "biosample_accession": first["biosample_accession"],
        "experiment_accessions": experiments,
        "run_count": len(sample_rows),
        "run_accessions": [row["run_accession"] for row in sample_rows],
        "run_manifest_sha256": manifest_hash,
        "reads": {
            "read_1": {"path": str(read_1), "bytes": read_1.stat().st_size, "sha256": sha256_file(read_1)},
            "read_2": {"path": str(read_2), "bytes": read_2.stat().st_size, "sha256": sha256_file(read_2)},
            "paired_fastq_records": record_count,
        },
        "reports": [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in reports
        ],
        "completed_at_utc": utc_now(),
        "status": "complete",
    }


def persist(args: argparse.Namespace) -> None:
    rows = read_run_manifest(args.run_manifest)
    manifest_hash = sha256_file(args.run_manifest)
    selected = [row for row in rows if row["sample_id"] == args.sample_id]
    if not selected:
        raise CheckpointError(f"sample {args.sample_id!r} is absent from the frozen manifest")
    root = args.checkpoint_dir.expanduser().resolve()
    claim_checkpoint_root(root, rows, manifest_hash)
    require_active_checkpoint_lifecycle(root)
    reads_dir = root / "reads"
    reports_dir = root / "reports" / args.sample_id
    records_dir = root / "records"
    target_1 = safe_child(reads_dir / f"{args.sample_id}_host_removed_R1.fastq.gz", root)
    target_2 = safe_child(reads_dir / f"{args.sample_id}_host_removed_R2.fastq.gz", root)
    record_path = safe_child(records_dir / f"{args.sample_id}.checkpoint.json", root)

    # A valid completion record is immutable.  This makes a replay of the
    # uncached persistence process idempotent and prevents a retry from
    # replacing an already trusted pair with different upstream bytes.
    if record_path.is_symlink():
        raise CheckpointError("durable completion record must not be a symlink")
    if record_path.is_file():
        try:
            with record_path.open("r", encoding="utf-8") as handle:
                existing_record = json.load(handle)
            validate_record(existing_record, selected, manifest_hash, root)
        except (OSError, json.JSONDecodeError, CheckpointError, KeyError, TypeError, ValueError):
            pass
        else:
            if args.output_record:
                atomic_json(args.output_record, existing_record)
            print(
                f"Validated immutable checkpoint already exists for {args.sample_id}; "
                "no managed files were replaced."
            )
            return

    # Validate every source before the first managed target is changed.  If a
    # prior record is absent or invalid, any interrupted copy remains pending
    # until a later invocation successfully commits a new record.
    source_record_count = validate_pair(args.read_1, args.read_2)
    report_sources = sorted(
        (path.expanduser().resolve() for path in args.report), key=str
    )
    seen_names: set[str] = set()
    for report in report_sources:
        if not report.is_file():
            raise CheckpointError(f"report does not exist: {report}")
        if report.name in seen_names:
            raise CheckpointError(f"duplicate checkpoint report basename: {report.name}")
        seen_names.add(report.name)

    copy_atomic(args.read_1, target_1)
    copy_atomic(args.read_2, target_2)
    record_count = validate_pair(target_1, target_2)
    if record_count != source_record_count:
        raise CheckpointError("checkpoint FASTQ record count changed during persistence")

    persisted_reports: list[Path] = []
    for report in report_sources:
        target = safe_child(reports_dir / report.name, root)
        copy_atomic(report, target)
        persisted_reports.append(target)

    record = sample_record(
        selected,
        manifest_hash,
        target_1,
        target_2,
        record_count,
        persisted_reports,
    )
    atomic_json(record_path, record)
    if args.output_record:
        atomic_json(args.output_record, record)


def expected_record_values(sample_rows: Sequence[dict[str, str]], manifest_hash: str) -> dict[str, Any]:
    first = sample_rows[0]
    return {
        "project_accession": first["project_accession"],
        "sample_order": int(first["sample_order"]),
        "sample_id": first["sample_id"],
        "identity_source": first["identity_source"],
        "biosample_accession": first["biosample_accession"],
        "experiment_accessions": sorted({row["experiment_accession"] for row in sample_rows if row["experiment_accession"]}),
        "run_count": len(sample_rows),
        "run_accessions": [row["run_accession"] for row in sample_rows],
        "run_manifest_sha256": manifest_hash,
        "status": "complete",
    }


def validate_record(
    record: dict[str, Any],
    sample_rows: Sequence[dict[str, str]],
    manifest_hash: str,
    checkpoint_root: Path,
) -> dict[str, Any]:
    if not isinstance(record, dict):
        raise CheckpointError("checkpoint record must be a JSON object")
    if record.get("schema_version") != SCHEMA_VERSION:
        raise CheckpointError("unsupported checkpoint schema version")
    for field, expected in expected_record_values(sample_rows, manifest_hash).items():
        if record.get(field) != expected:
            raise CheckpointError(f"checkpoint field {field!r} does not match the frozen manifest")
    try:
        reads = record["reads"]
        if not isinstance(reads, dict):
            raise TypeError
        read_metadata_1 = reads["read_1"]
        read_metadata_2 = reads["read_2"]
        if not isinstance(read_metadata_1, dict) or not isinstance(read_metadata_2, dict):
            raise TypeError
        read_1 = safe_child(Path(read_metadata_1["path"]), checkpoint_root)
        read_2 = safe_child(Path(read_metadata_2["path"]), checkpoint_root)
    except (KeyError, TypeError) as exc:
        raise CheckpointError("checkpoint record has invalid read metadata") from exc
    for label, path in (("read_1", read_1), ("read_2", read_2)):
        metadata = reads[label]
        if not path.is_file() or path.stat().st_size <= 0:
            raise CheckpointError(f"checkpoint {label} is missing or empty")
        if path.stat().st_size != metadata.get("bytes"):
            raise CheckpointError(f"checkpoint {label} size differs from its completion record")
        if sha256_file(path) != metadata.get("sha256"):
            raise CheckpointError(f"checkpoint {label} SHA-256 differs from its completion record")
    paired_records = validate_pair(read_1, read_2)
    if paired_records != reads.get("paired_fastq_records"):
        raise CheckpointError("checkpoint paired record count differs from its completion record")
    reports: list[str] = []
    raw_reports = record.get("reports", [])
    if not isinstance(raw_reports, list):
        raise CheckpointError("checkpoint reports metadata must be a list")
    seen_reports: set[Path] = set()
    for raw_report in raw_reports:
        if not isinstance(raw_report, dict):
            raise CheckpointError("checkpoint report metadata must be an object")
        raw_path = raw_report.get("path")
        if not isinstance(raw_path, str):
            raise CheckpointError("checkpoint report path must be a string")
        report = safe_child(Path(raw_path), checkpoint_root)
        if report in seen_reports:
            raise CheckpointError(f"checkpoint report metadata is duplicated: {report}")
        seen_reports.add(report)
        if not report.is_file() or report.stat().st_size <= 0:
            raise CheckpointError(f"checkpoint report is missing or empty: {report}")
        try:
            expected_report_size = int(raw_report.get("bytes"))
        except (TypeError, ValueError) as exc:
            raise CheckpointError(f"checkpoint report has an invalid byte count: {report}") from exc
        if expected_report_size <= 0 or report.stat().st_size != expected_report_size:
            raise CheckpointError(f"checkpoint report size differs from its completion record: {report}")
        if sha256_file(report) != raw_report.get("sha256"):
            raise CheckpointError(f"checkpoint report SHA-256 differs from its completion record: {report}")
        reports.append(str(report))
    return {
        "schema_version": record["schema_version"],
        "project_accession": record["project_accession"],
        "sample_order": record["sample_order"],
        "sample_id": record["sample_id"],
        "identity_source": record["identity_source"],
        "biosample_accession": record["biosample_accession"],
        "experiment_accessions": ";".join(record["experiment_accessions"]),
        "run_count": record["run_count"],
        "run_accessions": ";".join(record["run_accessions"]),
        "run_manifest_sha256": record["run_manifest_sha256"],
        "read_1": str(read_1),
        "read_1_bytes": record["reads"]["read_1"]["bytes"],
        "read_1_sha256": record["reads"]["read_1"]["sha256"],
        "read_2": str(read_2),
        "read_2_bytes": record["reads"]["read_2"]["bytes"],
        "read_2_sha256": record["reads"]["read_2"]["sha256"],
        "paired_fastq_records": paired_records,
        "reports_json": json.dumps(reports, separators=(",", ":")),
        "completed_at_utc": record["completed_at_utc"],
        "status": "complete",
    }


def reconcile(args: argparse.Namespace) -> int:
    rows = read_run_manifest(args.run_manifest)
    manifest_hash = sha256_file(args.run_manifest)
    root = args.checkpoint_dir.expanduser().resolve()
    claim_checkpoint_root(root, rows, manifest_hash)
    require_active_checkpoint_lifecycle(root)
    complete: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    for sample_id, sample_rows in group_samples(rows):
        record_path = root / "records" / f"{sample_id}.checkpoint.json"
        if not record_path.is_file():
            pending.append({"sample_order": sample_rows[0]["sample_order"], "sample_id": sample_id, "reason": "completion_record_missing"})
            continue
        try:
            with record_path.open("r", encoding="utf-8") as handle:
                record = json.load(handle)
            complete.append(validate_record(record, sample_rows, manifest_hash, root))
        except (OSError, json.JSONDecodeError, CheckpointError, KeyError, TypeError, ValueError) as exc:
            pending.append({"sample_order": sample_rows[0]["sample_order"], "sample_id": sample_id, "reason": str(exc).replace("\t", " ").replace("\n", " ")})
    complete.sort(key=lambda row: (int(row["sample_order"]), row["sample_id"]))
    pending.sort(key=lambda row: (int(row["sample_order"]), row["sample_id"]))
    atomic_write(args.output_manifest, tsv_bytes(complete, CHECKPOINT_FIELDS))
    atomic_write(args.pending_output, tsv_bytes(pending, PENDING_FIELDS))
    status = {
        "schema_version": SCHEMA_VERSION,
        "project_accession": rows[0]["project_accession"],
        "expected_samples": len(group_samples(rows)),
        "complete_samples": len(complete),
        "pending_samples": len(pending),
        "complete": not pending,
        "checked_at_utc": utc_now(),
    }
    if args.status_output:
        atomic_json(args.status_output, status)
    print(json.dumps(status, sort_keys=True))
    return 2 if args.require_complete and pending else 0


def validate_sample(args: argparse.Namespace) -> None:
    """Revalidate one durable sample before its disposable work is removed."""
    rows = read_run_manifest(args.run_manifest)
    sample_rows = [row for row in rows if row["sample_id"] == args.sample_id]
    if not sample_rows:
        raise CheckpointError(
            f"sample {args.sample_id!r} is absent from the frozen manifest"
        )
    root = args.checkpoint_dir.expanduser().resolve()
    manifest_hash = sha256_file(args.run_manifest)
    claim_checkpoint_root(root, rows, manifest_hash)
    require_active_checkpoint_lifecycle(root)
    record_path = safe_child(
        root / "records" / f"{args.sample_id}.checkpoint.json",
        root,
    )
    if not record_path.is_file():
        raise CheckpointError(
            f"durable completion record is missing for sample {args.sample_id}"
        )
    try:
        with record_path.open("r", encoding="utf-8") as handle:
            record = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckpointError(
            f"cannot read durable completion record for sample {args.sample_id}: {exc}"
        ) from exc
    validated = validate_record(record, sample_rows, manifest_hash, root)
    print(json.dumps(validated, sort_keys=True))


def load_global_success_marker(path: Path) -> dict[str, Any]:
    """Load a global-success marker without following a marker symlink."""
    marker = Path(os.path.abspath(path.expanduser()))
    if marker.is_symlink():
        raise CheckpointError("global-success marker must not be a symlink")
    try:
        with marker.open("r", encoding="utf-8") as handle:
            success = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckpointError(f"cannot validate global-success marker: {exc}") from exc
    if not isinstance(success, dict):
        raise CheckpointError("global-success marker must be a JSON object")
    return success


def validate_global_success_marker(
    success: dict[str, Any],
    checkpoint_manifest: Path,
    *,
    project_accession: str | None = None,
    results_dir: Path | None = None,
) -> dict[str, Any]:
    """Validate the finalizer marker and its three baseline durable outputs."""
    if success.get("status") != "complete":
        raise CheckpointError("global-success marker does not record complete status")
    if success.get("schema_version") != SCHEMA_VERSION:
        raise CheckpointError("global-success marker has an unsupported schema version")
    marker_project = success.get("project_accession")
    if not isinstance(marker_project, str) or not marker_project:
        raise CheckpointError("global-success marker has no project accession")
    if project_accession is not None and marker_project != project_accession:
        raise CheckpointError(
            "global-success marker project differs from the requested BioProject"
        )

    manifest_description = success.get("checkpoint_manifest")
    if not isinstance(manifest_description, dict):
        raise CheckpointError("global-success marker is not bound to a checkpoint manifest")
    manifest = Path(os.path.abspath(checkpoint_manifest.expanduser()))
    if manifest.is_symlink():
        raise CheckpointError(
            "checkpoint manifest must be a durable regular copy, not a symlink"
        )
    if not manifest.is_file():
        raise CheckpointError("checkpoint manifest is missing")
    described_manifest = Path(
        str(manifest_description.get("path", ""))
    ).expanduser()
    if (
        not described_manifest.is_absolute()
        or described_manifest.resolve() != manifest.resolve()
    ):
        raise CheckpointError(
            "checkpoint manifest path differs from the global-success marker"
        )
    if manifest.stat().st_size != manifest_description.get("bytes"):
        raise CheckpointError(
            "checkpoint manifest size differs from the global-success marker"
        )
    if sha256_file(manifest) != manifest_description.get("sha256"):
        raise CheckpointError(
            "checkpoint manifest SHA-256 differs from the global-success marker"
        )

    results_root: Path | None = None
    if results_dir is not None:
        unresolved_root = Path(os.path.abspath(results_dir.expanduser()))
        if unresolved_root.is_symlink():
            raise CheckpointError("scientific results root must not be a symlink")
        if not unresolved_root.is_dir():
            raise CheckpointError("scientific results root is missing")
        results_root = unresolved_root.resolve()

    outputs = success.get("outputs")
    if not isinstance(outputs, dict):
        raise CheckpointError("global-success marker has no durable output descriptions")
    for label in ("multiqc_report", "software_versions", "mag_abundance"):
        description = outputs.get(label)
        if not isinstance(description, dict):
            raise CheckpointError(f"global-success marker has no {label} description")
        output = Path(str(description.get("path", ""))).expanduser()
        if not output.is_absolute():
            raise CheckpointError(f"durable global output is missing: {label}")
        if output.is_symlink():
            raise CheckpointError(
                f"durable global output must be a regular copy, not a symlink: {label}"
            )
        if not output.is_file():
            raise CheckpointError(f"durable global output is missing: {label}")
        if results_root is not None:
            expected = (results_root / BASELINE_OUTPUT_PATHS[label]).resolve()
            if output.resolve() != expected:
                raise CheckpointError(
                    f"durable global output path differs from the results root: {label}"
                )
        if output.stat().st_size != description.get("bytes"):
            raise CheckpointError(f"durable global output size differs: {label}")
        if sha256_file(output) != description.get("sha256"):
            raise CheckpointError(f"durable global output SHA-256 differs: {label}")
    return manifest_description


def _inventory_regular_files(directory: Path, results_root: Path) -> list[dict[str, Any]]:
    """Hash one result subtree deterministically, rejecting every special entry."""
    descriptions: list[dict[str, Any]] = []
    try:
        with os.scandir(directory) as iterator:
            entries = sorted(iterator, key=lambda entry: entry.name)
    except OSError as exc:
        raise CheckpointError(f"cannot inspect scientific result directory {directory}: {exc}") from exc
    for entry in entries:
        path = Path(entry.path)
        relative = path.relative_to(results_root).as_posix()
        if entry.is_symlink():
            raise CheckpointError(
                f"scientific result inventory refuses symbolic link: {relative}"
            )
        try:
            metadata_before = path.lstat()
        except OSError as exc:
            raise CheckpointError(f"cannot inspect scientific artifact {relative}: {exc}") from exc
        if stat.S_ISDIR(metadata_before.st_mode):
            descriptions.extend(_inventory_regular_files(path, results_root))
            metadata_after = path.lstat()
            if (
                not stat.S_ISDIR(metadata_after.st_mode)
                or metadata_before.st_dev != metadata_after.st_dev
                or metadata_before.st_ino != metadata_after.st_ino
                or metadata_before.st_mtime_ns != metadata_after.st_mtime_ns
            ):
                raise CheckpointError(
                    "scientific result directory changed while its inventory was "
                    f"created: {relative}"
                )
            continue
        if not stat.S_ISREG(metadata_before.st_mode):
            raise CheckpointError(
                f"scientific result inventory refuses non-regular entry: {relative}"
            )
        digest = sha256_file(path)
        try:
            metadata_after = path.lstat()
        except OSError as exc:
            raise CheckpointError(f"cannot revalidate scientific artifact {relative}: {exc}") from exc
        if (
            metadata_before.st_size != metadata_after.st_size
            or metadata_before.st_mtime_ns != metadata_after.st_mtime_ns
            or metadata_before.st_dev != metadata_after.st_dev
            or metadata_before.st_ino != metadata_after.st_ino
            or not stat.S_ISREG(metadata_after.st_mode)
        ):
            raise CheckpointError(
                f"scientific artifact changed while its inventory was created: {relative}"
            )
        descriptions.append(
            {
                "relative_path": relative,
                "bytes": metadata_after.st_size,
                "sha256": digest,
            }
        )
    return descriptions


def build_scientific_inventory(results_dir: Path) -> dict[str, Any]:
    """Describe every durable scientific file and enforce required deliverables."""
    unresolved_root = Path(os.path.abspath(results_dir.expanduser()))
    if unresolved_root.is_symlink():
        raise CheckpointError("scientific results root must not be a symlink")
    if not unresolved_root.is_dir():
        raise CheckpointError("scientific results root is missing")
    results_root = unresolved_root.resolve()
    descriptions: list[dict[str, Any]] = []
    for root_name in SCIENTIFIC_RESULT_ROOTS:
        scientific_root = results_root / root_name
        if scientific_root.is_symlink():
            raise CheckpointError(
                f"scientific result root must not be a symlink: {root_name}"
            )
        if not scientific_root.is_dir():
            raise CheckpointError(f"mandatory scientific result root is missing: {root_name}")
        root_files = _inventory_regular_files(scientific_root, results_root)
        if not root_files:
            raise CheckpointError(
                f"mandatory scientific result root contains no regular files: {root_name}"
            )
        descriptions.extend(root_files)

    pipeline_info = results_root / "pipeline_info"
    if pipeline_info.is_symlink():
        raise CheckpointError("pipeline_info must not be a symlink")
    if not pipeline_info.is_dir():
        raise CheckpointError("pipeline_info directory is missing")
    software_versions = pipeline_info / "software_versions.tsv"
    if software_versions.is_symlink():
        raise CheckpointError(
            "scientific result inventory refuses symbolic link: pipeline_info/software_versions.tsv"
        )
    if not software_versions.is_file():
        raise CheckpointError("mandatory scientific artifact is missing: software_versions")
    software_metadata = software_versions.lstat()
    software_digest = sha256_file(software_versions)
    software_after = software_versions.lstat()
    if (
        not stat.S_ISREG(software_after.st_mode)
        or software_metadata.st_dev != software_after.st_dev
        or software_metadata.st_ino != software_after.st_ino
        or software_metadata.st_size != software_after.st_size
        or software_metadata.st_mtime_ns != software_after.st_mtime_ns
    ):
        raise CheckpointError(
            "scientific artifact changed while its inventory was created: "
            "pipeline_info/software_versions.tsv"
        )
    descriptions.append(
        {
            "relative_path": "pipeline_info/software_versions.tsv",
            "bytes": software_after.st_size,
            "sha256": software_digest,
        }
    )
    descriptions.sort(key=lambda item: item["relative_path"])

    required: dict[str, list[str]] = {}
    by_relative = {item["relative_path"]: item for item in descriptions}
    for label, patterns in REQUIRED_SCIENTIFIC_ARTIFACTS.items():
        matches = sorted(
            relative
            for relative in by_relative
            if any(PurePosixPath(relative).match(pattern) for pattern in patterns)
        )
        if not matches:
            raise CheckpointError(f"mandatory scientific artifact is missing: {label}")
        if any(int(by_relative[relative]["bytes"]) <= 0 for relative in matches):
            raise CheckpointError(f"mandatory scientific artifact is empty: {label}")
        required[label] = matches

    final_mag_ids = {
        PurePosixPath(relative).stem for relative in required["final_mags"]
    }
    annotation_suffix = ".functional_annotations.tsv"
    annotation_paths = required["functional_annotations"]
    annotation_ids = {
        PurePosixPath(relative).name[: -len(annotation_suffix)]
        for relative in annotation_paths
    }
    if annotation_ids != final_mag_ids or len(annotation_paths) != len(final_mag_ids):
        raise CheckpointError(
            "functional annotation MAG identifiers do not exactly match the final catalog"
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "results_root": str(results_root),
        "file_count": len(descriptions),
        "total_bytes": sum(int(item["bytes"]) for item in descriptions),
        "files": descriptions,
        "required_artifacts": required,
    }


def _validate_abundance_number(
    raw: str | None, metric: str, sample_id: str, mag_id: str
) -> None:
    value = (raw or "").strip()
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise CheckpointError(
            f"final MAG abundance contains non-numeric {metric} for "
            f"sample {sample_id!r}, MAG {mag_id!r}"
        ) from exc
    if not number.is_finite():
        raise CheckpointError(
            f"final MAG abundance contains non-finite {metric} for "
            f"sample {sample_id!r}, MAG {mag_id!r}"
        )
    if metric == "genome_length":
        valid = number > 0
    elif metric == "covered_fraction":
        valid = Decimal(0) <= number <= Decimal(1)
    elif metric == "relative_abundance_percent":
        valid = Decimal(0) <= number <= Decimal(100)
    else:
        valid = number >= 0
    if not valid:
        raise CheckpointError(
            f"final MAG abundance contains out-of-range {metric} for "
            f"sample {sample_id!r}, MAG {mag_id!r}"
        )


def validate_final_mag_abundance(
    results_root: Path,
    checkpoint_manifest: Path,
    project_accession: str,
) -> None:
    """Require one valid abundance row for every checkpoint-sample/MAG pair."""
    try:
        with checkpoint_manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(CHECKPOINT_FIELDS):
                raise CheckpointError("checkpoint manifest has an unexpected header")
            checkpoint_rows = list(reader)
    except (OSError, csv.Error) as exc:
        raise CheckpointError(f"cannot read checkpoint manifest: {exc}") from exc
    if not checkpoint_rows:
        raise CheckpointError("checkpoint manifest contains no completed samples")
    if any(
        row.get("project_accession") != project_accession
        or row.get("status") != "complete"
        for row in checkpoint_rows
    ):
        raise CheckpointError(
            "checkpoint manifest project/status does not match final MAG abundance"
        )
    sample_ids = [row.get("sample_id", "") for row in checkpoint_rows]
    if len(set(sample_ids)) != len(sample_ids) or any(
        not SAMPLE_ID_PATTERN.fullmatch(sample_id) for sample_id in sample_ids
    ):
        raise CheckpointError(
            "checkpoint manifest contains duplicate or unsafe sample identifiers"
        )
    expected_samples = set(sample_ids)

    final_mag_dir = (
        results_root / "02_mag_construction" / "final_catalog" / "final_catalog"
    )
    final_mag_paths = sorted(final_mag_dir.glob("*.fa"), key=lambda path: path.name)
    final_mag_ids = [path.stem for path in final_mag_paths]
    if not final_mag_ids or len(set(final_mag_ids)) != len(final_mag_ids):
        raise CheckpointError(
            "final MAG catalog has no unique identifiers for abundance validation"
        )
    expected_mags = set(final_mag_ids)

    abundance_path = (
        results_root
        / "05_mag_abundance_estimation"
        / "final_catalog.mag_abundance.long.tsv"
    )
    try:
        with abundance_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(ABUNDANCE_FIELDS):
                raise CheckpointError(
                    "final MAG abundance table has an unexpected header"
                )
            seen_pairs: set[tuple[str, str]] = set()
            row_count = 0
            for row_number, row in enumerate(reader, start=2):
                row_count += 1
                if None in row:
                    raise CheckpointError(
                        f"final MAG abundance row {row_number} contains unexpected fields"
                    )
                sample_id = (row.get("sample") or "").strip()
                mag_id = (row.get("mag_id") or "").strip()
                if sample_id not in expected_samples or mag_id not in expected_mags:
                    raise CheckpointError(
                        "final MAG abundance contains a row outside the checkpoint-sample "
                        f"and final-MAG cohorts: sample={sample_id!r}, mag_id={mag_id!r}"
                    )
                pair = (sample_id, mag_id)
                if pair in seen_pairs:
                    raise CheckpointError(
                        "final MAG abundance contains a duplicate sample/MAG pair: "
                        f"sample={sample_id!r}, mag_id={mag_id!r}"
                    )
                seen_pairs.add(pair)
                for metric in ABUNDANCE_FIELDS[2:]:
                    _validate_abundance_number(row.get(metric), metric, sample_id, mag_id)
    except (OSError, csv.Error) as exc:
        raise CheckpointError(f"cannot read final MAG abundance table: {exc}") from exc
    if row_count == 0:
        raise CheckpointError("final MAG abundance table contains no data rows")

    expected_count = len(expected_samples) * len(expected_mags)
    if len(seen_pairs) != expected_count:
        raise CheckpointError(
            "final MAG abundance does not contain the exact checkpoint-sample x "
            f"final-MAG matrix: expected {expected_count} rows, found {len(seen_pairs)}"
        )


def validate_scientific_inventory(
    success: dict[str, Any], expected_results_dir: Path | None = None
) -> Path:
    """Rebuild and exactly compare the sealed scientific output inventory."""
    sealed = success.get("scientific_outputs")
    if not isinstance(sealed, dict):
        raise CheckpointError(
            "global-success marker has no sealed scientific output inventory"
        )
    expected_keys = {
        "schema_version",
        "results_root",
        "file_count",
        "total_bytes",
        "files",
        "required_artifacts",
        "sealed_at_utc",
    }
    if set(sealed) != expected_keys or sealed.get("schema_version") != SCHEMA_VERSION:
        raise CheckpointError("sealed scientific output inventory has an invalid schema")
    results_root_text = sealed.get("results_root")
    if not isinstance(results_root_text, str):
        raise CheckpointError("sealed scientific output inventory has no results root")
    results_root = Path(results_root_text).expanduser()
    if not results_root.is_absolute():
        raise CheckpointError("sealed scientific output inventory results root is not absolute")
    if expected_results_dir is not None:
        expected_root = Path(os.path.abspath(expected_results_dir.expanduser()))
        if results_root.resolve() != expected_root.resolve():
            raise CheckpointError(
                "sealed scientific output inventory is bound to a different results root"
            )
    sealed_at = sealed.get("sealed_at_utc")
    if not isinstance(sealed_at, str) or not re.fullmatch(
        r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", sealed_at
    ):
        raise CheckpointError("sealed scientific output inventory has an invalid timestamp")

    current = build_scientific_inventory(results_root)
    expected = {**current, "sealed_at_utc": sealed_at}
    if sealed != expected:
        raise CheckpointError(
            "sealed scientific output inventory differs from current results "
            "(missing, tampered, or extra artifact)"
        )

    for label, relative in BASELINE_OUTPUT_PATHS.items():
        description = success["outputs"][label]
        described = Path(str(description.get("path", ""))).expanduser()
        expected_output = results_root / relative
        if described.resolve() != expected_output.resolve():
            raise CheckpointError(
                f"durable global output path differs from sealed results: {label}"
            )
    manifest_description = success.get("checkpoint_manifest")
    if not isinstance(manifest_description, dict):
        raise CheckpointError("global-success marker is not bound to a checkpoint manifest")
    checkpoint_manifest = Path(
        str(manifest_description.get("path", ""))
    ).expanduser()
    validate_final_mag_abundance(
        results_root.resolve(),
        checkpoint_manifest,
        str(success.get("project_accession", "")),
    )
    return results_root.resolve()


def seal_global(args: argparse.Namespace) -> None:
    """Extend the finalizer marker with an immutable complete-results seal."""
    success = load_global_success_marker(args.success_marker)
    validate_global_success_marker(
        success,
        args.checkpoint_manifest,
        project_accession=args.project_accession,
        results_dir=args.results_dir,
    )
    if "scientific_outputs" in success:
        validate_scientific_inventory(success, args.results_dir)
        print("Scientific output inventory was already sealed and validated.")
        return

    inventory = build_scientific_inventory(args.results_dir)
    validate_final_mag_abundance(
        Path(str(inventory["results_root"])),
        args.checkpoint_manifest,
        args.project_accession,
    )
    # Recheck the finalizer's baseline outputs after the full hashing pass so
    # a concurrent publication cannot be committed into a contradictory seal.
    validate_global_success_marker(
        success,
        args.checkpoint_manifest,
        project_accession=args.project_accession,
        results_dir=args.results_dir,
    )
    success["scientific_outputs"] = {
        **inventory,
        "sealed_at_utc": utc_now(),
    }
    atomic_json(args.success_marker, success)
    print(
        json.dumps(
            {
                "file_count": inventory["file_count"],
                "results_root": inventory["results_root"],
                "total_bytes": inventory["total_bytes"],
            },
            sort_keys=True,
        )
    )


def cleanup(args: argparse.Namespace) -> None:
    success = load_global_success_marker(args.success_marker)
    manifest_description = validate_global_success_marker(
        success, args.checkpoint_manifest
    )
    validate_scientific_inventory(success)
    try:
        with args.checkpoint_manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(CHECKPOINT_FIELDS):
                raise CheckpointError("checkpoint manifest has an unexpected header")
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        raise CheckpointError(f"cannot read checkpoint manifest: {exc}") from exc
    if not rows:
        raise CheckpointError("checkpoint manifest contains no completed samples")
    projects = {row["project_accession"] for row in rows}
    if projects != {success.get("project_accession")}:
        raise CheckpointError("checkpoint manifest project does not match global-success marker")
    if any(row["status"] != "complete" for row in rows):
        raise CheckpointError("checkpoint manifest contains a non-complete row")
    sample_ids = [row["sample_id"] for row in rows]
    if len(set(sample_ids)) != len(sample_ids) or any(
        not SAMPLE_ID_PATTERN.fullmatch(sample_id) for sample_id in sample_ids
    ):
        raise CheckpointError("checkpoint manifest contains duplicate or unsafe sample identifiers")
    root = args.checkpoint_dir.expanduser().resolve()
    manifest_hashes = {row["run_manifest_sha256"] for row in rows}
    if len(manifest_hashes) != 1:
        raise CheckpointError(
            "checkpoint manifest contains multiple frozen-manifest hashes"
        )
    validate_checkpoint_owner(
        root, str(success.get("project_accession")), next(iter(manifest_hashes))
    )
    cleanup_record_path = root / "sra_checkpoint_cleanup.json"
    previous_cleanup: dict[str, Any] | None = None
    if cleanup_record_path.is_file():
        try:
            loaded_cleanup = json.loads(cleanup_record_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise CheckpointError(f"cannot validate prior checkpoint cleanup record: {exc}") from exc
        if not isinstance(loaded_cleanup, dict):
            raise CheckpointError("prior checkpoint cleanup record must be a JSON object")
        if (
            loaded_cleanup.get("schema_version") != SCHEMA_VERSION
            or loaded_cleanup.get("project_accession") != success.get("project_accession")
        ):
            raise CheckpointError("prior checkpoint cleanup record does not match this project")
        previous_cleanup = loaded_cleanup

    targets: list[tuple[dict[str, str], str, Path, int, str]] = []
    for row in rows:
        for field in ("read_1", "read_2"):
            path = safe_child(Path(row[field]), root)
            mate = 1 if field == "read_1" else 2
            if path.name != f"{row['sample_id']}_host_removed_R{mate}.fastq.gz":
                raise CheckpointError(f"cleanup safety check rejected {path}")
            try:
                expected_size = int(row[f"{field}_bytes"])
            except (KeyError, TypeError, ValueError) as exc:
                raise CheckpointError(
                    f"cleanup found an invalid expected byte count for {path}"
                ) from exc
            if expected_size <= 0:
                raise CheckpointError(
                    f"cleanup found a non-positive expected byte count for {path}"
                )
            expected_hash = row[f"{field}_sha256"]
            if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
                raise CheckpointError(f"cleanup found an invalid expected SHA-256 for {path}")
            targets.append((row, field, path, expected_size, expected_hash))

    planned_files = [
        {
            "sample_id": row["sample_id"],
            "path": str(path),
            "bytes": expected_size,
            "sha256": expected_hash,
        }
        for row, _field, path, expected_size, expected_hash in targets
    ]
    planned_by_path = {item["path"]: item for item in planned_files}

    def validate_intact(path: Path, expected_size: int, expected_hash: str) -> None:
        if not path.is_file():
            raise CheckpointError(f"cleanup refused because a checkpoint read is missing: {path}")
        if path.stat().st_size != expected_size:
            raise CheckpointError(f"cleanup refused size-mismatched checkpoint: {path}")
        if sha256_file(path) != expected_hash:
            raise CheckpointError(f"cleanup refused hash-mismatched checkpoint: {path}")

    removed_by_path: dict[str, dict[str, Any]] = {}
    cleanup_record: dict[str, Any]
    if previous_cleanup is None:
        for _row, _field, path, expected_size, expected_hash in targets:
            validate_intact(path, expected_size, expected_hash)
        if args.keep:
            print("Checkpoint retention requested; no checkpoint reads were removed.")
            return
        now = utc_now()
        cleanup_record = {
            "schema_version": SCHEMA_VERSION,
            "project_accession": success.get("project_accession", ""),
            "status": "in_progress",
            "started_at_utc": now,
            "updated_at_utc": now,
            "checkpoint_manifest": manifest_description,
            "planned_files": planned_files,
            "removed_files": [],
            "removed_bytes": 0,
            "reports_and_provenance_retained": True,
        }
        # Persist the complete, validated deletion intent before the first
        # unlink.  If the process dies between unlink and the next journal
        # update, a later invocation can safely infer and finish that planned
        # removal instead of becoming permanently wedged.
        atomic_json(cleanup_record_path, cleanup_record)
    else:
        raw_removed = previous_cleanup.get("removed_files")
        if not isinstance(raw_removed, list) or not all(
            isinstance(item, dict) for item in raw_removed
        ):
            raise CheckpointError("prior checkpoint cleanup record has invalid removed_files")
        for item in raw_removed:
            raw_path = item.get("path")
            if not isinstance(raw_path, str) or raw_path not in planned_by_path:
                raise CheckpointError("prior cleanup record contains invalid path metadata")
            if raw_path in removed_by_path:
                raise CheckpointError("prior cleanup record contains duplicate paths")
            expected = planned_by_path[raw_path]
            try:
                parsed_bytes = int(item.get("bytes"))
            except (TypeError, ValueError) as exc:
                raise CheckpointError("prior cleanup record contains invalid byte metadata") from exc
            if (
                item.get("sample_id") != expected["sample_id"]
                or parsed_bytes != expected["bytes"]
                or (item.get("sha256") not in (None, expected["sha256"]))
            ):
                raise CheckpointError("prior cleanup record does not match the checkpoint manifest")
            removed_by_path[raw_path] = dict(expected)

        cleanup_status = previous_cleanup.get("status")
        if cleanup_status in (None, "complete"):
            # ``status``/``planned_files`` were added with the recoverable
            # journal.  Continue to accept a fully matching legacy completion
            # record produced by earlier repository revisions.
            if (
                "checkpoint_manifest" in previous_cleanup
                and previous_cleanup.get("checkpoint_manifest") != manifest_description
            ):
                raise CheckpointError(
                    "completed cleanup record is bound to a different manifest"
                )
            if (
                "planned_files" in previous_cleanup
                and previous_cleanup.get("planned_files") != planned_files
            ):
                raise CheckpointError(
                    "completed cleanup plan does not match the checkpoint manifest"
                )
            if "planned_files" in previous_cleanup and raw_removed != planned_files:
                raise CheckpointError(
                    "completed cleanup removals do not match its durable plan"
                )
            if set(removed_by_path) != set(planned_by_path):
                raise CheckpointError("prior cleanup record does not match the checkpoint manifest")
            if previous_cleanup.get("removed_bytes") != sum(
                int(item["bytes"]) for item in planned_files
            ):
                raise CheckpointError("prior cleanup byte total does not match the checkpoint manifest")
            if any(path.exists() for _row, _field, path, _size, _hash in targets):
                raise CheckpointError("checkpoint read exists despite a completed cleanup record")
            print("Checkpoint cleanup was already completed and validated; no files were removed.")
            return
        if cleanup_status != "in_progress":
            raise CheckpointError("prior checkpoint cleanup record has an unsupported status")
        if previous_cleanup.get("checkpoint_manifest") != manifest_description:
            raise CheckpointError("in-progress cleanup journal is bound to a different manifest")
        if previous_cleanup.get("planned_files") != planned_files:
            raise CheckpointError("in-progress cleanup journal plan does not match the manifest")
        if args.keep:
            raise CheckpointError(
                "checkpoint cleanup already started and cannot be changed to retention"
            )
        cleanup_record = previous_cleanup

    for row, _field, path, expected_size, expected_hash in targets:
        path_text = str(path)
        if path_text in removed_by_path:
            if path.exists():
                raise CheckpointError(
                    f"checkpoint read exists despite its recorded removal: {path}"
                )
            continue
        if path.exists():
            validate_intact(path, expected_size, expected_hash)
            path.unlink()
        # A missing, not-yet-recorded target is accepted only after a durable
        # in-progress journal exists: this is the crash window between unlink
        # and the journal update, and the path/hash were already validated in
        # the persisted deletion plan.
        removed_by_path[path_text] = dict(planned_by_path[path_text])
        removed = [
            removed_by_path[item["path"]]
            for item in planned_files
            if item["path"] in removed_by_path
        ]
        cleanup_record["removed_files"] = removed
        cleanup_record["removed_bytes"] = sum(int(item["bytes"]) for item in removed)
        cleanup_record["updated_at_utc"] = utc_now()
        atomic_json(cleanup_record_path, cleanup_record)

    cleanup_record["status"] = "complete"
    cleanup_record["cleaned_at_utc"] = utc_now()
    cleanup_record["updated_at_utc"] = cleanup_record["cleaned_at_utc"]
    cleanup_record["removed_files"] = planned_files
    cleanup_record["removed_bytes"] = sum(int(item["bytes"]) for item in planned_files)
    atomic_json(cleanup_record_path, cleanup_record)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    create = commands.add_parser("persist", help="atomically persist one validated sample pair")
    create.add_argument("--run-manifest", required=True, type=Path)
    create.add_argument("--sample-id", required=True)
    create.add_argument("--read-1", required=True, type=Path)
    create.add_argument("--read-2", required=True, type=Path)
    create.add_argument("--checkpoint-dir", required=True, type=Path)
    create.add_argument("--report", action="append", default=[], type=Path)
    create.add_argument("--output-record", type=Path)

    check = commands.add_parser("reconcile", help="validate all records against a frozen cohort")
    check.add_argument("--run-manifest", required=True, type=Path)
    check.add_argument("--checkpoint-dir", required=True, type=Path)
    check.add_argument("--output-manifest", required=True, type=Path)
    check.add_argument("--pending-output", required=True, type=Path)
    check.add_argument("--status-output", type=Path)
    check.add_argument("--require-complete", action="store_true")

    sample_check = commands.add_parser(
        "validate-sample",
        help="revalidate one durable sample before deleting its disposable work",
    )
    sample_check.add_argument("--run-manifest", required=True, type=Path)
    sample_check.add_argument("--checkpoint-dir", required=True, type=Path)
    sample_check.add_argument("--sample-id", required=True)

    seal = commands.add_parser(
        "seal-global",
        help="seal every mandatory durable scientific result before read cleanup",
    )
    seal.add_argument("--success-marker", required=True, type=Path)
    seal.add_argument("--checkpoint-manifest", required=True, type=Path)
    seal.add_argument("--results-dir", required=True, type=Path)
    seal.add_argument("--project-accession", required=True)

    remove = commands.add_parser("cleanup", help="remove read pairs only after global success")
    remove.add_argument("--checkpoint-manifest", required=True, type=Path)
    remove.add_argument("--checkpoint-dir", required=True, type=Path)
    remove.add_argument("--success-marker", required=True, type=Path)
    remove.add_argument("--keep", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "sample_id", None) and not SAMPLE_ID_PATTERN.fullmatch(args.sample_id):
        print("ERROR: invalid sample identifier", file=sys.stderr)
        return 2
    try:
        if args.command == "persist":
            persist(args)
            return 0
        if args.command == "reconcile":
            return reconcile(args)
        if args.command == "validate-sample":
            validate_sample(args)
            return 0
        if args.command == "seal-global":
            seal_global(args)
            return 0
        if args.command == "cleanup":
            cleanup(args)
            return 0
        raise CheckpointError(f"unsupported command: {args.command}")
    except CheckpointError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except (OSError, csv.Error, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        print(f"ERROR: checkpoint lifecycle I/O or data failure: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
