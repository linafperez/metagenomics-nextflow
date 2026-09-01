#!/usr/bin/env python3
"""Merge Nextflow traces and create deterministic resource-accounting reports.

The input registry makes every Nextflow invocation explicit.  Cached trace rows
remain visible for provenance but contribute no resource consumption, which is
essential when a run is resumed or when the pipeline is split into several
Nextflow invocations.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import re
import statistics
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


MISSING = {"", "-", "na", "n/a", "nan", "null", "none"}
NUMBER_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
DURATION_PART_RE = re.compile(
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*"
    r"(milliseconds?|msecs?|ms|seconds?|secs?|s|minutes?|mins?|m|hours?|hrs?|h|days?|d)",
    re.IGNORECASE,
)
MEMORY_RE = re.compile(
    r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*"
    r"(bytes?|b|kib|kb|mib|mb|gib|gb|tib|tb|pib|pb)?$",
    re.IGNORECASE,
)


TASK_FIELDS = (
    "invocation_id",
    "session_id",
    "invocation_stage",
    "trace_row",
    "task_id",
    "hash",
    "native_id",
    "process",
    "module",
    "parent_subworkflow",
    "subworkflow",
    "tag",
    "name",
    "status",
    "exit_code",
    "attempt",
    "submit_time",
    "start_time",
    "complete_time",
    "duration_seconds",
    "realtime_seconds",
    "requested_cpus",
    "requested_memory_bytes",
    "requested_time_seconds",
    "cpu_percent",
    "memory_percent",
    "memory_efficiency_percent",
    "allocated_cpu_hours",
    "observed_cpu_hours",
    "peak_rss_bytes",
    "peak_vmem_bytes",
    "rss_bytes",
    "vmem_bytes",
    "rchar_bytes",
    "wchar_bytes",
    "read_bytes",
    "write_bytes",
    "voluntary_context_switches",
    "involuntary_context_switches",
    "task_peak_work_bytes",
    "accelerator_requested",
    "accelerator_type",
    "gpu_hours",
    "gpu_models",
    "gpu_utilization_mean_percent",
    "gpu_utilization_max_percent",
    "peak_gpu_memory_bytes",
    "gpu_metric_samples",
    "accounted",
    "workdir",
    "queue",
    "hostname",
    "container",
    "error_action",
)


AGGREGATE_FIELDS = (
    "scope",
    "subworkflow",
    "parent_subworkflow",
    "process",
    "module",
    "process_count",
    "task_runs",
    "executed_tasks",
    "cached_tasks",
    "completed_tasks",
    "failed_tasks",
    "aborted_tasks",
    "unique_task_hashes",
    "cumulative_realtime_seconds",
    "mean_realtime_seconds",
    "median_realtime_seconds",
    "max_realtime_seconds",
    "activity_wall_time_seconds",
    "requested_cpu_measurements",
    "max_requested_cpus",
    "median_requested_cpus",
    "mean_requested_cpus",
    "requested_memory_measurements",
    "max_requested_memory_bytes",
    "median_requested_memory_bytes",
    "mean_requested_memory_bytes",
    "requested_time_measurements",
    "max_requested_time_seconds",
    "median_requested_time_seconds",
    "mean_requested_time_seconds",
    "allocated_cpu_hours",
    "observed_cpu_hours",
    "cpu_efficiency_percent",
    "cpu_percent_measurements",
    "peak_rss_measurements",
    "max_peak_rss_bytes",
    "median_peak_rss_bytes",
    "mean_peak_rss_bytes",
    "memory_efficiency_measurements",
    "max_memory_efficiency_percent",
    "median_memory_efficiency_percent",
    "mean_memory_efficiency_percent",
    "max_peak_vmem_bytes",
    "task_disk_measurements",
    "max_task_peak_work_bytes",
    "median_task_peak_work_bytes",
    "mean_task_peak_work_bytes",
    "sampled_concurrent_work_measurements",
    "sampled_peak_concurrent_work_bytes",
    "sampled_stage_storage_measurements",
    "sampled_stage_peak_dynamic_storage_bytes",
    "rchar_bytes",
    "wchar_bytes",
    "read_bytes",
    "write_bytes",
    "requested_accelerator_task_total",
    "max_requested_accelerators",
    "accelerator_hours",
    "gpu_models",
    "gpu_utilization_mean_percent",
    "gpu_utilization_max_percent",
    "peak_gpu_memory_bytes",
    "gpu_metric_samples",
)


OUTER_SCOPE_TOKENS = {
    "QUALITY_CONTROL_AND_FILTERING": "quality_control_and_filtering",
    "MAG_CONSTRUCTION": "mag_construction",
    "TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS": (
        "taxonomic_classification_and_phylogenomics"
    ),
    "GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION": (
        "gene_prediction_and_functional_annotation"
    ),
    "MAG_ABUNDANCE_ESTIMATION": "mag_abundance_estimation",
    "GLOBAL_PROCESSING_EVALUATION": "global_processing_evaluation",
    "SRA_ACQUISITION_AND_PREPROCESSING": "sra_acquisition_and_preprocessing",
    "SRA_ACQUISITION": "sra_acquisition_and_preprocessing",
    "SRA_SAMPLE_PREPROCESSING": "sra_acquisition_and_preprocessing",
    "SRA_PREPROCESSING": "sra_acquisition_and_preprocessing",
}

INPUT_DIRECT_TASKS = {
    "CHECK_SAMPLESHEET",
    "RESOLVE_SRA_PROJECT",
    "RESOLVE_SRA_INPUT",
    "VALIDATE_SRA_MANIFEST",
    "VALIDATE_SRA_PROJECT",
    "CHECK_SRA_CHECKPOINTS",
}
REPORTING_DIRECT_TASKS = {
    "COLLECT_VERSIONS",
    "SUMMARIZE_RESOURCES",
    "MONITOR_STORAGE",
    "FINALIZE_SRA_GLOBAL_RUN",
}


@dataclass(frozen=True)
class Invocation:
    invocation_id: str
    trace_path: Path
    session_id: str
    stage: str
    started_at: str
    finished_at: str
    trace_raw: bool
    launch_dir: Path
    status: str = ""
    exit_code: int | None = None


@dataclass
class TaskRecord:
    values: dict[str, Any]
    raw: dict[str, str]


@dataclass
class BuildContext:
    limitations: set[str] = field(default_factory=set)
    warnings: set[str] = field(default_factory=set)
    raw_columns: set[str] = field(default_factory=set)
    skipped_duplicate_rows: int = 0


def clean_text(value: object) -> str:
    return "" if value is None else str(value).strip()


def is_missing(value: object) -> bool:
    return clean_text(value).lower() in MISSING


def parse_float(value: object) -> float | None:
    text = clean_text(value).replace(",", "")
    if text.lower() in MISSING:
        return None
    if not NUMBER_RE.fullmatch(text):
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def parse_int(value: object) -> int | None:
    number = parse_float(value)
    if number is None:
        return None
    return int(number)


def parse_bool(value: object, *, default: bool) -> bool:
    text = clean_text(value).lower()
    if not text:
        return default
    if text in {"true", "yes", "y", "1"}:
        return True
    if text in {"false", "no", "n", "0"}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def parse_duration(value: object, *, raw_numeric: bool = True) -> float | None:
    text = clean_text(value).lower().replace(",", "")
    if text in MISSING:
        return None
    if NUMBER_RE.fullmatch(text):
        number = float(text)
        return number / 1000.0 if raw_numeric else number

    if ":" in text and re.fullmatch(r"\d+(?::\d+){1,2}(?:\.\d+)?", text):
        parts = [float(part) for part in text.split(":")]
        if len(parts) == 3:
            return parts[0] * 3600 + parts[1] * 60 + parts[2]
        return parts[0] * 60 + parts[1]

    units = {
        "ms": 0.001,
        "msec": 0.001,
        "msecs": 0.001,
        "millisecond": 0.001,
        "milliseconds": 0.001,
        "s": 1.0,
        "sec": 1.0,
        "secs": 1.0,
        "second": 1.0,
        "seconds": 1.0,
        "m": 60.0,
        "min": 60.0,
        "mins": 60.0,
        "minute": 60.0,
        "minutes": 60.0,
        "h": 3600.0,
        "hr": 3600.0,
        "hrs": 3600.0,
        "hour": 3600.0,
        "hours": 3600.0,
        "d": 86400.0,
        "day": 86400.0,
        "days": 86400.0,
    }
    total = 0.0
    consumed: list[tuple[int, int]] = []
    for match in DURATION_PART_RE.finditer(text):
        total += float(match.group(1)) * units[match.group(2).lower()]
        consumed.append(match.span())
    if not consumed:
        return None
    residue = list(text)
    for start, end in consumed:
        residue[start:end] = " " * (end - start)
    if "".join(residue).strip():
        return None
    return total


def parse_bytes(value: object) -> int | None:
    text = clean_text(value).replace(",", "")
    if text.lower() in MISSING:
        return None
    match = MEMORY_RE.fullmatch(text)
    if not match:
        return None
    number = float(match.group(1))
    unit = (match.group(2) or "b").lower()
    powers = {
        "b": 0,
        "byte": 0,
        "bytes": 0,
        "kb": 1,
        "kib": 1,
        "mb": 2,
        "mib": 2,
        "gb": 3,
        "gib": 3,
        "tb": 4,
        "tib": 4,
        "pb": 5,
        "pib": 5,
    }
    return int(round(number * (1024 ** powers[unit])))


def parse_percent(value: object) -> float | None:
    text = clean_text(value)
    if text.endswith("%"):
        text = text[:-1].strip()
    return parse_float(text)


def parse_timestamp(value: object) -> datetime | None:
    text = clean_text(value)
    if text.lower() in MISSING:
        return None
    if NUMBER_RE.fullmatch(text):
        number = float(text)
        # Nextflow raw trace dates are epoch milliseconds.
        if abs(number) > 100_000_000_000:
            number /= 1000.0
        try:
            return datetime.fromtimestamp(number, tz=timezone.utc)
        except (OverflowError, OSError, ValueError):
            return None
    normalized = text[:-1] + "+00:00" if text.endswith("Z") else text
    for parser in (
        lambda item: datetime.fromisoformat(item),
        lambda item: datetime.strptime(item, "%Y-%m-%d %H:%M:%S.%f"),
        lambda item: datetime.strptime(item, "%Y-%m-%d %H:%M:%S"),
    ):
        try:
            parsed = parser(normalized)
        except ValueError:
            continue
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc)
    return None


def timestamp_text(value: object) -> str:
    parsed = parse_timestamp(value)
    if parsed is None:
        return clean_text(value) if not is_missing(value) else ""
    return parsed.isoformat(timespec="milliseconds").replace("+00:00", "Z")


def parse_accelerator_count(value: object) -> float:
    number = parse_float(value)
    if number is not None:
        return max(0.0, number)
    text = clean_text(value)
    if not text:
        return 0.0
    match = re.search(r"(?:request|count|limit)\s*[:=]\s*(\d+(?:\.\d+)?)", text)
    return float(match.group(1)) if match else 0.0


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return format(value, ".12g")
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return str(value)


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    from io import StringIO

    output = StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: fmt(row.get(field)) for field in fields})
    atomic_write_text(path, output.getvalue())


def read_registry(path: Path) -> list[Invocation]:
    registry_dir = path.parent.absolute()
    invocations: dict[str, Invocation] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"invocation_id", "trace_path"}.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                "trace registry must contain invocation_id and trace_path columns"
            )
        for line_number, row in enumerate(reader, start=2):
            invocation_id = clean_text(row.get("invocation_id"))
            trace_value = clean_text(row.get("trace_path"))
            if not invocation_id or not trace_value:
                raise ValueError(f"incomplete trace registry row {line_number}")
            trace_path = Path(trace_value)
            if not trace_path.is_absolute():
                trace_path = registry_dir / trace_path
            launch_value = clean_text(row.get("launch_dir"))
            launch_dir = Path(launch_value) if launch_value else registry_dir
            if not launch_dir.is_absolute():
                launch_dir = registry_dir / launch_dir
            invocation = Invocation(
                invocation_id=invocation_id,
                trace_path=Path(os.path.abspath(trace_path)),
                session_id=clean_text(row.get("session_id")),
                stage=clean_text(row.get("stage")),
                started_at=clean_text(row.get("started_at")),
                finished_at=clean_text(row.get("finished_at")),
                trace_raw=parse_bool(row.get("trace_raw"), default=True),
                launch_dir=Path(os.path.abspath(launch_dir)),
                status=clean_text(row.get("status")).lower(),
                exit_code=parse_int(row.get("exit_code")),
            )
            previous = invocations.get(invocation_id)
            if previous is None:
                invocations[invocation_id] = invocation
            elif previous != invocation:
                raise ValueError(
                    f"conflicting registry entries for invocation_id {invocation_id!r}"
                )
    return [invocations[key] for key in sorted(invocations)]


def outer_scope(process: str, stage: str = "") -> str:
    tokens = [part.strip() for part in process.split(":") if part.strip()]
    upper_tokens = [part.upper() for part in tokens]
    stage_upper = stage.upper().replace("-", "_").replace(" ", "_")
    if stage_upper == "SRA_PREPROCESS":
        # One staged invocation owns the full acquisition -> QC -> trimming ->
        # host-removal -> checkpoint lifecycle for a biological sample.
        return "sra_acquisition_and_preprocessing"
    for token in upper_tokens:
        if token in OUTER_SCOPE_TOKENS:
            return OUTER_SCOPE_TOKENS[token]
    terminal = upper_tokens[-1] if upper_tokens else ""
    if terminal in INPUT_DIRECT_TASKS:
        return "input_validation"
    if terminal in REPORTING_DIRECT_TASKS:
        return "pipeline_reporting"
    if "SRA" in stage_upper:
        return "sra_acquisition_and_preprocessing"
    return "unmapped"


def scope_parts(process: str, stage: str) -> tuple[str, str, str]:
    tokens = [part.strip() for part in process.split(":") if part.strip()]
    module = tokens[-1] if tokens else process
    parent = tokens[-2] if len(tokens) > 1 else ""
    return module, parent, outer_scope(process, stage)


def normalise_workdir(value: object, invocation: Invocation) -> str:
    text = clean_text(value)
    if not text:
        return ""
    path = Path(text)
    if not path.is_absolute():
        path = invocation.launch_dir / path
    return os.path.abspath(path)


def raw_value(row: Mapping[str, str], *names: str) -> str:
    for name in names:
        if name in row:
            return clean_text(row.get(name))
    return ""


def build_task(
    row: dict[str, str], row_number: int, invocation: Invocation
) -> TaskRecord:
    process = raw_value(row, "process", "name")
    module, parent, subworkflow = scope_parts(process, invocation.stage)
    status = raw_value(row, "status").upper()
    accounted = status != "CACHED"
    realtime = parse_duration(raw_value(row, "realtime"), raw_numeric=invocation.trace_raw)
    duration = parse_duration(raw_value(row, "duration"), raw_numeric=invocation.trace_raw)
    cpus = parse_float(raw_value(row, "cpus"))
    cpu_percent = parse_percent(raw_value(row, "%cpu", "cpu"))
    accelerator = parse_accelerator_count(raw_value(row, "accelerator"))
    allocated_cpu_hours = (
        cpus * realtime / 3600.0
        if accounted and cpus is not None and realtime is not None
        else 0.0
    )
    observed_cpu_hours = (
        (cpu_percent / 100.0) * realtime / 3600.0
        if accounted and cpu_percent is not None and realtime is not None
        else 0.0
    )
    gpu_hours = (
        accelerator * realtime / 3600.0
        if accounted and realtime is not None
        else 0.0
    )
    requested_memory = parse_bytes(raw_value(row, "memory"))
    peak_rss = parse_bytes(raw_value(row, "peak_rss"))
    memory_efficiency = (
        peak_rss / requested_memory * 100.0
        if accounted
        and peak_rss is not None
        and requested_memory is not None
        and requested_memory > 0
        else None
    )
    values: dict[str, Any] = {
        "invocation_id": invocation.invocation_id,
        "session_id": invocation.session_id,
        "invocation_stage": invocation.stage,
        "trace_row": row_number,
        "task_id": raw_value(row, "task_id"),
        "hash": raw_value(row, "hash"),
        "native_id": raw_value(row, "native_id"),
        "process": process,
        "module": module,
        "parent_subworkflow": parent,
        "subworkflow": subworkflow,
        "tag": raw_value(row, "tag"),
        "name": raw_value(row, "name"),
        "status": status,
        "exit_code": parse_int(raw_value(row, "exit")),
        "attempt": parse_int(raw_value(row, "attempt")),
        "submit_time": timestamp_text(raw_value(row, "submit")),
        "start_time": timestamp_text(raw_value(row, "start")),
        "complete_time": timestamp_text(raw_value(row, "complete")),
        "duration_seconds": duration,
        "realtime_seconds": realtime,
        "requested_cpus": cpus,
        "requested_memory_bytes": requested_memory,
        "requested_time_seconds": parse_duration(
            raw_value(row, "time"), raw_numeric=invocation.trace_raw
        ),
        "cpu_percent": cpu_percent,
        "memory_percent": parse_percent(raw_value(row, "%mem")),
        "memory_efficiency_percent": memory_efficiency,
        "allocated_cpu_hours": allocated_cpu_hours,
        "observed_cpu_hours": observed_cpu_hours,
        "peak_rss_bytes": peak_rss,
        "peak_vmem_bytes": parse_bytes(raw_value(row, "peak_vmem")),
        "rss_bytes": parse_bytes(raw_value(row, "rss")),
        "vmem_bytes": parse_bytes(raw_value(row, "vmem")),
        "rchar_bytes": parse_bytes(raw_value(row, "rchar")),
        "wchar_bytes": parse_bytes(raw_value(row, "wchar")),
        "read_bytes": parse_bytes(raw_value(row, "read_bytes")),
        "write_bytes": parse_bytes(raw_value(row, "write_bytes")),
        "voluntary_context_switches": parse_int(raw_value(row, "vol_ctxt")),
        "involuntary_context_switches": parse_int(raw_value(row, "inv_ctxt")),
        "task_peak_work_bytes": None,
        "accelerator_requested": accelerator,
        "accelerator_type": raw_value(row, "accelerator_type"),
        "gpu_hours": gpu_hours,
        "gpu_models": "",
        "gpu_utilization_mean_percent": None,
        "gpu_utilization_max_percent": None,
        "peak_gpu_memory_bytes": None,
        "gpu_metric_samples": 0,
        "accounted": accounted,
        "workdir": normalise_workdir(raw_value(row, "workdir"), invocation),
        "queue": raw_value(row, "queue"),
        "hostname": raw_value(row, "hostname"),
        "container": raw_value(row, "container"),
        "error_action": raw_value(row, "error_action"),
    }
    return TaskRecord(values, dict(row))


def dedupe_key(task: TaskRecord) -> tuple[Any, ...]:
    values = task.values
    if values["session_id"]:
        return (
            "session",
            values["session_id"],
            values["task_id"],
            values["attempt"],
            values["hash"],
            values["process"],
            values["start_time"],
            values["complete_time"],
            values["status"],
        )
    return ("invocation", values["invocation_id"], values["trace_row"])


def load_tasks(invocations: Sequence[Invocation], context: BuildContext) -> list[TaskRecord]:
    tasks: list[TaskRecord] = []
    seen: set[tuple[Any, ...]] = set()
    for invocation in invocations:
        if not invocation.trace_path.exists():
            raise FileNotFoundError(f"trace not found: {invocation.trace_path}")
        with invocation.trace_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                context.warnings.add(f"empty trace: {invocation.trace_path}")
                continue
            context.raw_columns.update(reader.fieldnames)
            for row_number, row in enumerate(reader, start=2):
                task = build_task(
                    {key: value or "" for key, value in row.items()},
                    row_number,
                    invocation,
                )
                key = dedupe_key(task)
                if key in seen:
                    context.skipped_duplicate_rows += 1
                    continue
                seen.add(key)
                tasks.append(task)

    tasks.sort(
        key=lambda task: (
            task.values["invocation_id"],
            task.values["process"],
            _sort_number(task.values["task_id"]),
            task.values["attempt"] or 0,
            task.values["trace_row"],
        )
    )
    return tasks


def _sort_number(value: object) -> tuple[int, float | str]:
    number = parse_float(value)
    return (0, number) if number is not None else (1, clean_text(value))


def load_task_peaks(path: Path | None) -> dict[str, int]:
    if path is None or not path.exists():
        return {}
    peaks: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"workdir", "peak_work_bytes"}.issubset(
            reader.fieldnames
        ):
            raise ValueError("task peak file requires workdir and peak_work_bytes")
        for row in reader:
            workdir = clean_text(row.get("workdir"))
            size = parse_int(row.get("peak_work_bytes"))
            if workdir and size is not None:
                key = os.path.abspath(workdir)
                peaks[key] = max(peaks.get(key, 0), size)
    return peaks


def attach_task_peaks(tasks: Sequence[TaskRecord], peaks: Mapping[str, int]) -> None:
    by_suffix: dict[str, list[int]] = defaultdict(list)
    for path, size in peaks.items():
        parts = Path(path).parts
        if len(parts) >= 2:
            by_suffix[os.path.join(parts[-2], parts[-1])].append(size)
    for task in tasks:
        workdir = task.values["workdir"]
        if workdir in peaks:
            task.values["task_peak_work_bytes"] = peaks[workdir]
            continue
        if workdir:
            parts = Path(workdir).parts
            suffix = os.path.join(parts[-2], parts[-1]) if len(parts) >= 2 else ""
            candidates = by_suffix.get(suffix, [])
            if len(candidates) == 1:
                task.values["task_peak_work_bytes"] = candidates[0]


def load_gpu_metrics(
    path: Path | None, context: BuildContext
) -> dict[tuple[str, str, str, int], dict[str, Any]]:
    grouped: dict[tuple[str, str, str, int], dict[str, Any]] = {}
    if path is None or not path.exists():
        return grouped
    files = sorted(path.glob("*.gpu_metrics.tsv")) if path.is_dir() else [path]
    for metric_file in files:
        try:
            with metric_file.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {
                    "process",
                    "sample_id",
                    "session_id",
                    "attempt",
                    "gpu_name",
                    "utilization_gpu_percent",
                    "memory_used_mib",
                }
                if not reader.fieldnames or not required.issubset(reader.fieldnames):
                    context.warnings.add(f"GPU metric file has an incompatible header: {metric_file}")
                    continue
                for row in reader:
                    process = clean_text(row.get("process")).upper()
                    sample_id = clean_text(row.get("sample_id"))
                    session_id = clean_text(row.get("session_id"))
                    attempt = parse_int(row.get("attempt"))
                    utilization = parse_float(row.get("utilization_gpu_percent"))
                    memory_mib = parse_float(row.get("memory_used_mib"))
                    if (
                        not process
                        or not sample_id
                        or not session_id
                        or attempt is None
                        or utilization is None
                        or memory_mib is None
                    ):
                        continue
                    item = grouped.setdefault(
                        (session_id, process, sample_id, attempt),
                        {"utilizations": [], "peak_memory_bytes": 0, "models": set()},
                    )
                    item["utilizations"].append(utilization)
                    item["peak_memory_bytes"] = max(item["peak_memory_bytes"], int(memory_mib * 1024 * 1024))
                    model = clean_text(row.get("gpu_name"))
                    if model:
                        item["models"].add(model)
        except OSError as exc:
            context.warnings.add(f"Cannot read GPU metric file {metric_file}: {exc}")
    return grouped


def attach_gpu_metrics(
    tasks: Sequence[TaskRecord],
    metrics: Mapping[tuple[str, str, str, int], Mapping[str, Any]],
) -> None:
    for task in tasks:
        session_id = clean_text(task.values.get("session_id"))
        process = clean_text(task.values.get("module")).upper()
        tag = clean_text(task.values.get("tag"))
        attempt = task.values.get("attempt")
        if not session_id or not process or not tag or attempt is None:
            continue
        value = metrics.get((session_id, process, tag, int(attempt)))
        if value is None:
            continue
        utilizations = value["utilizations"]
        task.values["gpu_models"] = ";".join(sorted(value["models"]))
        task.values["gpu_utilization_mean_percent"] = statistics.fmean(utilizations) if utilizations else None
        task.values["gpu_utilization_max_percent"] = max(utilizations) if utilizations else None
        task.values["peak_gpu_memory_bytes"] = value["peak_memory_bytes"]
        task.values["gpu_metric_samples"] = len(utilizations)


def non_null(values: Iterable[Any]) -> list[Any]:
    return [value for value in values if value is not None]


def total_metric(tasks: Iterable[TaskRecord], field: str) -> float | int:
    return sum(
        task.values[field]
        for task in tasks
        if task.values["accounted"] and task.values[field] is not None
    )


def elapsed_between(tasks: Sequence[TaskRecord]) -> float | None:
    starts = non_null(parse_timestamp(task.values["start_time"]) for task in tasks)
    finishes = non_null(parse_timestamp(task.values["complete_time"]) for task in tasks)
    if not starts or not finishes:
        return None
    elapsed = (max(finishes) - min(starts)).total_seconds()
    return max(0.0, elapsed)


def aggregate_row(
    scope: str, process: str, subworkflow: str, grouped: Sequence[TaskRecord]
) -> dict[str, Any]:
    executed = [task for task in grouped if task.values["accounted"]]
    statuses = [task.values["status"] for task in grouped]
    realtime = non_null(task.values["realtime_seconds"] for task in executed)
    rss = non_null(task.values["peak_rss_bytes"] for task in executed)
    requested_cpus = non_null(task.values["requested_cpus"] for task in executed)
    requested_memory = non_null(
        task.values["requested_memory_bytes"] for task in executed
    )
    requested_time = non_null(
        task.values["requested_time_seconds"] for task in executed
    )
    memory_efficiency = non_null(
        task.values["memory_efficiency_percent"] for task in executed
    )
    vmem = non_null(task.values["peak_vmem_bytes"] for task in executed)
    disk = non_null(task.values["task_peak_work_bytes"] for task in executed)
    allocated = float(sum(task.values["allocated_cpu_hours"] for task in executed))
    observed = float(sum(task.values["observed_cpu_hours"] for task in executed))
    cpu_covered_allocated = sum(
        task.values["allocated_cpu_hours"]
        for task in executed
        if task.values["cpu_percent"] is not None
    )
    module_names = sorted({task.values["module"] for task in grouped})
    parent_names = sorted(
        {task.values["parent_subworkflow"] for task in grouped if task.values["parent_subworkflow"]}
    )
    gpu_sample_count = sum(task.values["gpu_metric_samples"] for task in executed)
    gpu_utilization_weighted = sum(
        task.values["gpu_utilization_mean_percent"] * task.values["gpu_metric_samples"]
        for task in executed
        if task.values["gpu_utilization_mean_percent"] is not None
    )
    gpu_models = sorted(
        {
            model
            for task in executed
            for model in clean_text(task.values["gpu_models"]).split(";")
            if model
        }
    )
    gpu_peak_memory = non_null(task.values["peak_gpu_memory_bytes"] for task in executed)
    gpu_max_utilization = non_null(task.values["gpu_utilization_max_percent"] for task in executed)
    return {
        "scope": scope,
        "subworkflow": subworkflow,
        "parent_subworkflow": parent_names[0] if len(parent_names) == 1 else ("MULTIPLE" if parent_names else ""),
        "process": process,
        "module": module_names[0] if len(module_names) == 1 else "MULTIPLE",
        "process_count": len({task.values["process"] for task in grouped}),
        "task_runs": len(grouped),
        "executed_tasks": len(executed),
        "cached_tasks": statuses.count("CACHED"),
        "completed_tasks": statuses.count("COMPLETED"),
        "failed_tasks": statuses.count("FAILED"),
        "aborted_tasks": statuses.count("ABORTED"),
        "unique_task_hashes": len(
            {task.values["hash"] for task in grouped if task.values["hash"]}
        ),
        "cumulative_realtime_seconds": sum(realtime),
        "mean_realtime_seconds": statistics.fmean(realtime) if realtime else None,
        "median_realtime_seconds": statistics.median(realtime) if realtime else None,
        "max_realtime_seconds": max(realtime) if realtime else None,
        "activity_wall_time_seconds": elapsed_between(executed),
        "requested_cpu_measurements": len(requested_cpus),
        "max_requested_cpus": max(requested_cpus) if requested_cpus else None,
        "median_requested_cpus": statistics.median(requested_cpus) if requested_cpus else None,
        "mean_requested_cpus": statistics.fmean(requested_cpus) if requested_cpus else None,
        "requested_memory_measurements": len(requested_memory),
        "max_requested_memory_bytes": max(requested_memory) if requested_memory else None,
        "median_requested_memory_bytes": statistics.median(requested_memory) if requested_memory else None,
        "mean_requested_memory_bytes": statistics.fmean(requested_memory) if requested_memory else None,
        "requested_time_measurements": len(requested_time),
        "max_requested_time_seconds": max(requested_time) if requested_time else None,
        "median_requested_time_seconds": statistics.median(requested_time) if requested_time else None,
        "mean_requested_time_seconds": statistics.fmean(requested_time) if requested_time else None,
        "allocated_cpu_hours": allocated,
        "observed_cpu_hours": observed,
        "cpu_efficiency_percent": (
            observed / cpu_covered_allocated * 100.0
            if cpu_covered_allocated > 0
            else None
        ),
        "cpu_percent_measurements": sum(
            task.values["cpu_percent"] is not None for task in executed
        ),
        "peak_rss_measurements": len(rss),
        "max_peak_rss_bytes": max(rss) if rss else None,
        "median_peak_rss_bytes": statistics.median(rss) if rss else None,
        "mean_peak_rss_bytes": statistics.fmean(rss) if rss else None,
        "memory_efficiency_measurements": len(memory_efficiency),
        "max_memory_efficiency_percent": max(memory_efficiency) if memory_efficiency else None,
        "median_memory_efficiency_percent": statistics.median(memory_efficiency) if memory_efficiency else None,
        "mean_memory_efficiency_percent": statistics.fmean(memory_efficiency) if memory_efficiency else None,
        "max_peak_vmem_bytes": max(vmem) if vmem else None,
        "task_disk_measurements": len(disk),
        "max_task_peak_work_bytes": max(disk) if disk else None,
        "median_task_peak_work_bytes": statistics.median(disk) if disk else None,
        "mean_task_peak_work_bytes": statistics.fmean(disk) if disk else None,
        "sampled_concurrent_work_measurements": 0,
        "sampled_peak_concurrent_work_bytes": None,
        "sampled_stage_storage_measurements": 0,
        "sampled_stage_peak_dynamic_storage_bytes": None,
        "rchar_bytes": total_metric(executed, "rchar_bytes"),
        "wchar_bytes": total_metric(executed, "wchar_bytes"),
        "read_bytes": total_metric(executed, "read_bytes"),
        "write_bytes": total_metric(executed, "write_bytes"),
        "requested_accelerator_task_total": sum(
            task.values["accelerator_requested"] for task in executed
        ),
        "max_requested_accelerators": max(
            (task.values["accelerator_requested"] for task in executed),
            default=0,
        ),
        "accelerator_hours": sum(task.values["gpu_hours"] for task in executed),
        "gpu_models": ";".join(gpu_models),
        "gpu_utilization_mean_percent": (
            gpu_utilization_weighted / gpu_sample_count if gpu_sample_count else None
        ),
        "gpu_utilization_max_percent": max(gpu_max_utilization) if gpu_max_utilization else None,
        "peak_gpu_memory_bytes": max(gpu_peak_memory) if gpu_peak_memory else None,
        "gpu_metric_samples": gpu_sample_count,
    }


def aggregate_by_process(tasks: Sequence[TaskRecord]) -> list[dict[str, Any]]:
    groups: dict[str, list[TaskRecord]] = defaultdict(list)
    for task in tasks:
        groups[task.values["process"]].append(task)
    rows = []
    for process in sorted(groups):
        grouped = groups[process]
        subworkflows = sorted({task.values["subworkflow"] for task in grouped})
        subworkflow = subworkflows[0] if len(subworkflows) == 1 else "multiple"
        rows.append(aggregate_row("process", process, subworkflow, grouped))
    return rows


def aggregate_by_subworkflow(tasks: Sequence[TaskRecord]) -> list[dict[str, Any]]:
    groups: dict[str, list[TaskRecord]] = defaultdict(list)
    for task in tasks:
        groups[task.values["subworkflow"]].append(task)
    return [
        aggregate_row("subworkflow", "", subworkflow, groups[subworkflow])
        for subworkflow in sorted(groups)
    ]


def attach_sampled_concurrent_work(
    path: Path | None,
    tasks: Sequence[TaskRecord],
    process_rows: Sequence[dict[str, Any]],
    subworkflow_rows: Sequence[dict[str, Any]],
    context: BuildContext,
) -> None:
    if path is None or not path.exists():
        context.limitations.add(
            "Task work-directory timeseries was not provided; sampled concurrent process/subworkflow disk peaks are unavailable."
        )
        return

    by_workdir: dict[str, tuple[str, str]] = {}
    by_suffix: dict[str, set[tuple[str, str]]] = defaultdict(set)
    ambiguous: set[str] = set()
    for task in tasks:
        workdir = clean_text(task.values.get("workdir"))
        if not workdir:
            continue
        key = os.path.normcase(os.path.abspath(workdir))
        scope = (task.values["process"], task.values["subworkflow"])
        parts = Path(key).parts
        if len(parts) >= 2:
            by_suffix[os.path.normcase(os.path.join(parts[-2], parts[-1]))].add(
                scope
            )
        previous = by_workdir.get(key)
        if previous is not None and previous != scope:
            ambiguous.add(key)
        else:
            by_workdir[key] = scope
    for key in ambiguous:
        by_workdir.pop(key, None)
    if ambiguous:
        context.warnings.add(
            f"Ignored {len(ambiguous)} task work directories with ambiguous scope attribution."
        )

    snapshots: dict[str, dict[str, int]] = defaultdict(dict)
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"timestamp", "workdir", "work_bytes"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(
                    "task work timeseries requires timestamp, workdir, and work_bytes"
                )
            for row in reader:
                timestamp = clean_text(row.get("timestamp"))
                workdir = clean_text(row.get("workdir"))
                size = parse_int(row.get("work_bytes"))
                if not timestamp or not workdir or size is None or size < 0:
                    continue
                key = os.path.normcase(os.path.abspath(workdir))
                snapshots[timestamp][key] = max(
                    snapshots[timestamp].get(key, 0), size
                )
    except (OSError, csv.Error, ValueError) as exc:
        context.warnings.add(f"Cannot read optional task work timeseries: {exc}")
        return

    process_samples: dict[str, list[int]] = defaultdict(list)
    subworkflow_samples: dict[str, list[int]] = defaultdict(list)
    for timestamp in sorted(snapshots):
        process_totals: dict[str, int] = defaultdict(int)
        subworkflow_totals: dict[str, int] = defaultdict(int)
        for workdir, size in snapshots[timestamp].items():
            scope = by_workdir.get(workdir)
            if scope is None:
                parts = Path(workdir).parts
                suffix = (
                    os.path.normcase(os.path.join(parts[-2], parts[-1]))
                    if len(parts) >= 2
                    else ""
                )
                candidates = by_suffix.get(suffix, set())
                scope = next(iter(candidates)) if len(candidates) == 1 else None
            if scope is None:
                continue
            process, subworkflow = scope
            process_totals[process] += size
            subworkflow_totals[subworkflow] += size
        for process, total in process_totals.items():
            process_samples[process].append(total)
        for subworkflow, total in subworkflow_totals.items():
            subworkflow_samples[subworkflow].append(total)

    for row in process_rows:
        samples = process_samples.get(str(row["process"]), [])
        row["sampled_concurrent_work_measurements"] = len(samples)
        row["sampled_peak_concurrent_work_bytes"] = max(samples) if samples else None
    for row in subworkflow_rows:
        samples = subworkflow_samples.get(str(row["subworkflow"]), [])
        row["sampled_concurrent_work_measurements"] = len(samples)
        row["sampled_peak_concurrent_work_bytes"] = max(samples) if samples else None

    context.limitations.add(
        "Sampled concurrent process/subworkflow work peaks sum observed task work directories at each monitor timestamp; they exclude external SRA scratch/cache/temp roots and can miss short-lived maxima between samples."
    )


def load_storage(path: Path | None, context: BuildContext) -> dict[str, Any]:
    empty = {
        "samples": 0,
        "peak_total_dynamic_bytes": None,
        "peak_total_measured_bytes": None,
        "peaks": {},
        "last": {},
        "by_invocation": {},
        "by_invocation_stage": {},
        "measurement_complete": None,
    }
    if path is None or not path.exists():
        context.limitations.add(
            "Global storage timeseries was not provided; peak and final global storage are unavailable."
        )
        return empty
    rows: list[dict[str, str]] = []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows.extend(
                {key: value or "" for key, value in row.items()} for row in reader
            )
    except (OSError, csv.Error) as exc:
        context.warnings.add(f"Cannot read optional global storage timeseries: {exc}")
        return empty
    if not rows:
        context.limitations.add("Global storage timeseries is empty.")
        return empty
    totals = non_null(parse_int(row.get("total_measured_bytes")) for row in rows)
    dynamic_totals = non_null(parse_int(row.get("total_dynamic_bytes")) for row in rows)
    if not dynamic_totals:
        # Compatibility with early telemetry fixtures: derive dynamic storage
        # only when both measured total and database bytes are present.
        dynamic_totals = non_null(
            total - database
            if (total := parse_int(row.get("total_measured_bytes"))) is not None
            and (database := parse_int(row.get("database_bytes"))) is not None
            else None
            for row in rows
        )
    complete = all(clean_text(row.get("measurement_complete")).lower() == "true" for row in rows)
    if not complete:
        context.limitations.add(
            "At least one storage sample was incomplete or used overlapping configured roots."
        )
    last = rows[-1]
    category_fields = (
        "work_bytes",
        "checkpoint_bytes",
        "sra_cache_bytes",
        "sra_scratch_bytes",
        "sra_temp_bytes",
        "results_bytes",
        "database_bytes",
    )
    def grouped_storage(grouped_rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
        grouped_dynamic = non_null(
            parse_int(row.get("total_dynamic_bytes")) for row in grouped_rows
        )
        grouped_measured = non_null(
            parse_int(row.get("total_measured_bytes")) for row in grouped_rows
        )
        return {
            "samples": len(grouped_rows),
            "peak_total_dynamic_bytes": (
                max(grouped_dynamic) if grouped_dynamic else None
            ),
            "peak_total_measured_bytes": (
                max(grouped_measured) if grouped_measured else None
            ),
            "peaks": {
                field: max(values) if values else None
                for field in category_fields
                if (
                    values := non_null(
                        parse_int(row.get(field)) for row in grouped_rows
                    )
                )
            },
            "measurement_complete": all(
                clean_text(row.get("measurement_complete")).lower() == "true"
                for row in grouped_rows
            ),
        }

    rows_by_invocation: dict[str, list[dict[str, str]]] = defaultdict(list)
    rows_by_stage: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        invocation_id = clean_text(row.get("invocation_id"))
        invocation_stage = clean_text(row.get("invocation_stage"))
        if invocation_id:
            rows_by_invocation[invocation_id].append(row)
        if invocation_stage:
            rows_by_stage[invocation_stage].append(row)
    return {
        "samples": len(rows),
        "peak_total_dynamic_bytes": max(dynamic_totals) if dynamic_totals else None,
        "peak_total_measured_bytes": max(totals) if totals else None,
        "peaks": {
            field: max(values) if values else None
            for field in category_fields
            if (values := non_null(parse_int(row.get(field)) for row in rows))
        },
        "last": {
            field: parse_int(last.get(field))
            for field in (*category_fields, "total_dynamic_bytes", "total_measured_bytes")
        },
        "by_invocation": {
            key: grouped_storage(rows_by_invocation[key])
            for key in sorted(rows_by_invocation)
        },
        "by_invocation_stage": {
            key: grouped_storage(rows_by_stage[key]) for key in sorted(rows_by_stage)
        },
        "measurement_complete": complete,
    }


def attach_stage_storage(
    subworkflow_rows: Sequence[dict[str, Any]], storage: Mapping[str, Any]
) -> None:
    stage_reports = storage.get("by_invocation_stage", {})
    if not isinstance(stage_reports, Mapping):
        return
    preprocessing = stage_reports.get("sra-preprocess")
    if not isinstance(preprocessing, Mapping):
        return
    for row in subworkflow_rows:
        if row.get("subworkflow") == "sra_acquisition_and_preprocessing":
            row["sampled_stage_storage_measurements"] = parse_int(
                preprocessing.get("samples")
            ) or 0
            row["sampled_stage_peak_dynamic_storage_bytes"] = parse_int(
                preprocessing.get("peak_total_dynamic_bytes")
            )


def _manifest_reader(handle: Any, path: Path) -> csv.DictReader:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    return csv.DictReader(handle, delimiter=delimiter)


def load_manifest_counts(
    path: Path | None, context: BuildContext, manifest_mode: str = "auto"
) -> dict[str, int | None]:
    if path is None or not path.exists():
        context.limitations.add(
            "An input manifest was not provided; biological sample and SRA run counts are unavailable."
        )
        return {"biological_samples": None, "sra_runs": None}
    sample_candidates = (
        "internal_sample_id",
        "sample_id",
        "sample",
        "biosample_accession",
        "biosample",
    )
    run_candidates = ("run_accession", "sra_run_accession", "run")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = _manifest_reader(handle, path)
        field_map = {field.lower(): field for field in (reader.fieldnames or [])}
        sample_field = next((field_map[name] for name in sample_candidates if name in field_map), None)
        run_field = next((field_map[name] for name in run_candidates if name in field_map), None)
        resolved_mode = manifest_mode
        if resolved_mode == "auto":
            resolved_mode = "sra" if run_field is not None else "local"
        samples: set[str] = set()
        runs: set[str] = set()
        for row in reader:
            if sample_field and clean_text(row.get(sample_field)):
                samples.add(clean_text(row.get(sample_field)))
            if run_field and clean_text(row.get(run_field)):
                runs.add(clean_text(row.get(run_field)))
    if sample_field is None:
        context.limitations.add("Input manifest has no recognized biological-sample column.")
    if resolved_mode == "sra" and run_field is None:
        context.limitations.add("SRA manifest has no recognized run-accession column.")
    return {
        "biological_samples": len(samples) if sample_field else None,
        "sra_runs": len(runs) if run_field else (0 if resolved_mode == "local" else None),
    }


def workflow_wall_time(
    invocations: Sequence[Invocation], tasks: Sequence[TaskRecord], context: BuildContext
) -> tuple[float | None, str, str, str]:
    registry_starts = non_null(parse_timestamp(item.started_at) for item in invocations)
    registry_finishes = non_null(parse_timestamp(item.finished_at) for item in invocations)
    source = "registry"
    if not registry_starts or not registry_finishes:
        starts = non_null(parse_timestamp(task.values["start_time"]) for task in tasks)
        finishes = non_null(parse_timestamp(task.values["complete_time"]) for task in tasks)
        registry_starts = registry_starts or starts
        registry_finishes = registry_finishes or finishes
        source = "trace bounds"
        context.limitations.add(
            "Workflow wall time was inferred from trace task bounds because complete invocation timestamps were unavailable."
        )
    if not registry_starts or not registry_finishes:
        return None, "", "", "unavailable"
    start = min(registry_starts)
    finish = max(registry_finishes)
    elapsed = max(0.0, (finish - start).total_seconds())
    return (
        elapsed,
        start.isoformat(timespec="milliseconds").replace("+00:00", "Z"),
        finish.isoformat(timespec="milliseconds").replace("+00:00", "Z"),
        source,
    )


def largest_row(rows: Sequence[Mapping[str, Any]], metric: str, label: str) -> str | None:
    candidates = [row for row in rows if row.get(metric) is not None]
    if not candidates:
        return None
    winner = sorted(candidates, key=lambda row: (-float(row[metric]), str(row[label])))[0]
    return str(winner[label])


def build_summary(
    invocations: Sequence[Invocation],
    tasks: Sequence[TaskRecord],
    process_rows: Sequence[Mapping[str, Any]],
    subworkflow_rows: Sequence[Mapping[str, Any]],
    storage: Mapping[str, Any],
    manifest_counts: Mapping[str, int | None],
    project_status: Mapping[str, Any],
    slurm_accounting: Mapping[str, Any],
    context: BuildContext,
) -> dict[str, Any]:
    executed = [task for task in tasks if task.values["accounted"]]
    cached = [task for task in tasks if not task.values["accounted"]]
    rss = non_null(task.values["peak_rss_bytes"] for task in executed)
    task_disk = non_null(task.values["task_peak_work_bytes"] for task in executed)
    realtime = non_null(task.values["realtime_seconds"] for task in executed)
    requested_cpus = non_null(task.values["requested_cpus"] for task in executed)
    requested_memory = non_null(
        task.values["requested_memory_bytes"] for task in executed
    )
    requested_time = non_null(
        task.values["requested_time_seconds"] for task in executed
    )
    memory_efficiency = non_null(
        task.values["memory_efficiency_percent"] for task in executed
    )
    cpu_measured = [task for task in executed if task.values["cpu_percent"] is not None]
    allocated = float(sum(task.values["allocated_cpu_hours"] for task in executed))
    observed = float(sum(task.values["observed_cpu_hours"] for task in executed))
    covered_allocated = sum(task.values["allocated_cpu_hours"] for task in cpu_measured)
    wall, start, finish, wall_source = workflow_wall_time(invocations, tasks, context)
    accelerator_hours = sum(task.values["gpu_hours"] for task in executed)
    accelerator_tasks = sum(task.values["accelerator_requested"] > 0 for task in executed)
    gpu_metric_samples = sum(task.values["gpu_metric_samples"] for task in executed)
    gpu_memory = non_null(task.values["peak_gpu_memory_bytes"] for task in executed)
    gpu_utilization = [
        (task.values["gpu_utilization_mean_percent"], task.values["gpu_metric_samples"])
        for task in executed
        if task.values["gpu_utilization_mean_percent"] is not None
    ]
    gpu_models = sorted({model for task in executed for model in clean_text(task.values["gpu_models"]).split(";") if model})

    if len(cpu_measured) < len(executed):
        context.limitations.add(
            "Observed CPU hours exclude tasks without %cpu; Nextflow utilization is sampled and can be inaccurate for short tasks."
        )
    if len(rss) < len(executed):
        context.limitations.add("RSS statistics exclude tasks without peak_rss measurements.")
    if len(task_disk) < len(executed):
        context.limitations.add(
            "Task disk statistics exclude work directories not observed by the storage monitor."
        )
    if any(task.values["subworkflow"] == "unmapped" for task in tasks):
        context.limitations.add(
            "Some processes did not match a known outer workflow scope and are grouped as unmapped."
        )
    context.limitations.add(
        "RSS is reported as maximum, median, and mean per task; values are never summed and do not estimate concurrent pipeline RAM."
    )
    context.limitations.add(
        "rchar/wchar are process system-call bytes and may include cache or pipes; read_bytes/write_bytes are block-I/O counters."
    )
    context.limitations.add(
        "Per-task peak work size and concurrent group work size are sampled; max_task_peak_work_bytes remains the largest individual task, while sampled_peak_concurrent_work_bytes sums same-snapshot task directories."
    )
    context.limitations.add(
        "Published output bytes cannot be attributed reliably to individual processes or subworkflows from portable Nextflow traces; final results size is measured only at project scope."
    )
    context.limitations.add(
        "The final results-size sample precedes generation of the resource summary files themselves."
    )
    if accelerator_tasks == 0:
        context.limitations.add(
            "No task requested a GPU/accelerator; GPU hours are zero and GPU memory/utilization are not measured."
        )
    elif not gpu_metric_samples:
        context.limitations.add(
            "GPU hours are requested accelerator count multiplied by task realtime; nvidia-smi samples were unavailable."
        )
    else:
        context.limitations.add(
            "GPU utilization and memory are device-level nvidia-smi samples for the single visible allocated device, not PID-isolated measurements."
        )
        tasks_with_gpu_metrics = sum(
            task.values["gpu_metric_samples"] > 0 for task in executed
        )
        if tasks_with_gpu_metrics < accelerator_tasks:
            context.limitations.add(
                "GPU utilization coverage excludes accelerator task attempts without a session/attempt-matched metric file."
            )

    return {
        "schema_version": "1.0",
        "project": {
            "status": clean_text(project_status.get("status")) or "unknown",
            "exit_code": parse_int(project_status.get("exit_code")),
            "status_updated_at_utc": clean_text(project_status.get("updated_at_utc")),
            "successful": clean_text(project_status.get("status")).lower() == "complete",
        },
        "counts": {
            "nextflow_invocations": len(invocations),
            "successful_nextflow_invocations": sum(item.status == "completed" for item in invocations),
            "failed_nextflow_invocations": sum(item.status == "failed" for item in invocations),
            "trace_rows": len(tasks),
            "executed_task_runs": len(executed),
            "cached_task_rows_excluded": len(cached),
            "duplicate_trace_rows_excluded": context.skipped_duplicate_rows,
            "completed_task_rows": sum(task.values["status"] == "COMPLETED" for task in tasks),
            "failed_task_rows": sum(task.values["status"] == "FAILED" for task in tasks),
            "aborted_task_rows": sum(task.values["status"] == "ABORTED" for task in tasks),
            "process_scopes": len(process_rows),
            "outer_workflow_scopes": len(subworkflow_rows),
            **manifest_counts,
        },
        "time": {
            "workflow_start": start,
            "workflow_finish": finish,
            "workflow_wall_time_seconds": wall,
            "wall_time_source": wall_source,
            "cumulative_task_realtime_seconds": sum(realtime),
            "cumulative_task_duration_seconds": total_metric(executed, "duration_seconds"),
        },
        "cpu": {
            "requested_cpu_task_measurements": len(requested_cpus),
            "max_requested_cpus_per_task": max(requested_cpus) if requested_cpus else None,
            "median_requested_cpus_per_task": statistics.median(requested_cpus) if requested_cpus else None,
            "mean_requested_cpus_per_task": statistics.fmean(requested_cpus) if requested_cpus else None,
            "allocated_cpu_hours": allocated,
            "observed_cpu_hours": observed,
            "observed_cpu_efficiency_percent": (
                observed / covered_allocated * 100.0 if covered_allocated > 0 else None
            ),
            "cpu_percent_task_coverage": len(cpu_measured),
            "executed_task_count": len(executed),
            "allocated_formula": "requested_cpus * realtime_hours",
            "observed_formula": "(cpu_percent / 100) * realtime_hours",
        },
        "memory": {
            "requested_memory_task_measurements": len(requested_memory),
            "max_requested_memory_bytes": max(requested_memory) if requested_memory else None,
            "median_requested_memory_bytes": statistics.median(requested_memory) if requested_memory else None,
            "mean_requested_memory_bytes": statistics.fmean(requested_memory) if requested_memory else None,
            "requested_time_task_measurements": len(requested_time),
            "max_requested_time_seconds": max(requested_time) if requested_time else None,
            "peak_rss_task_coverage": len(rss),
            "max_peak_rss_bytes": max(rss) if rss else None,
            "median_peak_rss_bytes": statistics.median(rss) if rss else None,
            "mean_peak_rss_bytes": statistics.fmean(rss) if rss else None,
            "memory_efficiency_task_coverage": len(memory_efficiency),
            "max_memory_efficiency_percent": max(memory_efficiency) if memory_efficiency else None,
            "median_memory_efficiency_percent": statistics.median(memory_efficiency) if memory_efficiency else None,
            "mean_memory_efficiency_percent": statistics.fmean(memory_efficiency) if memory_efficiency else None,
            "memory_efficiency_formula": "peak_rss_bytes / requested_memory_bytes * 100",
            "aggregation_rule": "max/median/mean per-task peak RSS; never sum",
        },
        "io": {
            "rchar_bytes": total_metric(executed, "rchar_bytes"),
            "wchar_bytes": total_metric(executed, "wchar_bytes"),
            "read_bytes": total_metric(executed, "read_bytes"),
            "write_bytes": total_metric(executed, "write_bytes"),
        },
        "storage": {
            **storage,
            "task_peak_coverage": len(task_disk),
            "max_task_peak_work_bytes": max(task_disk) if task_disk else None,
            "median_task_peak_work_bytes": (
                statistics.median(task_disk) if task_disk else None
            ),
            "mean_task_peak_work_bytes": (
                statistics.fmean(task_disk) if task_disk else None
            ),
        },
        "gpu": {
            "requested": accelerator_tasks > 0,
            "telemetry_observed": gpu_metric_samples > 0,
            "used": accelerator_tasks > 0,
            "tasks_requesting_accelerator": accelerator_tasks,
            "requested_accelerator_task_total": sum(
                task.values["accelerator_requested"] for task in executed
            ),
            "allocated_gpu_hours": accelerator_hours,
            "accelerator_hours": accelerator_hours,
            "allocated_formula": "requested_accelerators * task_realtime_hours",
            "models": gpu_models,
            "metric_samples": gpu_metric_samples,
            "observed_utilization_mean_percent": (
                sum(mean * count for mean, count in gpu_utilization) / sum(count for _mean, count in gpu_utilization)
                if gpu_utilization else None
            ),
            "peak_gpu_memory_bytes": max(gpu_memory) if gpu_memory else None,
        },
        "slurm_accounting": dict(slurm_accounting),
        "largest_consumers": {
            "process_by_allocated_cpu_hours": largest_row(
                process_rows, "allocated_cpu_hours", "process"
            ),
            "process_by_peak_rss": largest_row(
                process_rows, "max_peak_rss_bytes", "process"
            ),
            "process_by_task_peak_disk": largest_row(
                process_rows, "max_task_peak_work_bytes", "process"
            ),
            "process_by_sampled_concurrent_work": largest_row(
                process_rows, "sampled_peak_concurrent_work_bytes", "process"
            ),
            "process_by_cumulative_realtime": largest_row(
                process_rows, "cumulative_realtime_seconds", "process"
            ),
            "process_by_longest_task": largest_row(
                process_rows, "max_realtime_seconds", "process"
            ),
            "subworkflow_by_allocated_cpu_hours": largest_row(
                subworkflow_rows, "allocated_cpu_hours", "subworkflow"
            ),
            "subworkflow_by_task_peak_disk": largest_row(
                subworkflow_rows, "max_task_peak_work_bytes", "subworkflow"
            ),
            "subworkflow_by_sampled_concurrent_work": largest_row(
                subworkflow_rows,
                "sampled_peak_concurrent_work_bytes",
                "subworkflow",
            ),
            "subworkflow_by_sampled_stage_dynamic_storage": largest_row(
                subworkflow_rows,
                "sampled_stage_peak_dynamic_storage_bytes",
                "subworkflow",
            ),
            "subworkflow_by_cumulative_realtime": largest_row(
                subworkflow_rows, "cumulative_realtime_seconds", "subworkflow"
            ),
        },
        "warnings": sorted(context.warnings),
        "limitations": sorted(context.limitations),
    }


def flatten_summary(summary: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def visit(prefix: str, value: Any) -> None:
        if isinstance(value, Mapping):
            for key in sorted(value):
                visit(f"{prefix}.{key}" if prefix else str(key), value[key])
        elif isinstance(value, list):
            for index, item in enumerate(value, start=1):
                visit(f"{prefix}.{index}", item)
        else:
            rows.append({"metric": prefix, "value": value})

    visit("", summary)
    return rows


def render_html(
    summary: Mapping[str, Any],
    process_rows: Sequence[Mapping[str, Any]],
    subworkflow_rows: Sequence[Mapping[str, Any]],
) -> str:
    def cell(value: Any) -> str:
        return html.escape(fmt(value))

    def table(fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> str:
        header = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
        body = "".join(
            "<tr>" + "".join(f"<td>{cell(row.get(field))}</td>" for field in fields) + "</tr>"
            for row in rows
        )
        return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"

    overview = [
        {"metric": row["metric"], "value": row["value"]}
        for row in flatten_summary(summary)
        if not row["metric"].startswith(("limitations.", "warnings."))
    ]
    limitations = "".join(
        f"<li>{html.escape(str(item))}</li>" for item in summary["limitations"]
    )
    process_fields = (
        "process",
        "subworkflow",
        "executed_tasks",
        "cached_tasks",
        "cumulative_realtime_seconds",
        "allocated_cpu_hours",
        "observed_cpu_hours",
        "max_peak_rss_bytes",
        "max_task_peak_work_bytes",
        "sampled_peak_concurrent_work_bytes",
    )
    subworkflow_fields = (
        "subworkflow",
        "executed_tasks",
        "cached_tasks",
        "cumulative_realtime_seconds",
        "activity_wall_time_seconds",
        "allocated_cpu_hours",
        "observed_cpu_hours",
        "max_peak_rss_bytes",
        "max_task_peak_work_bytes",
        "sampled_peak_concurrent_work_bytes",
        "sampled_stage_peak_dynamic_storage_bytes",
    )
    return """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Metagenomics resource usage</title>
<style>
body{font-family:system-ui,sans-serif;margin:2rem;color:#1f2937}table{border-collapse:collapse;width:100%;margin-bottom:2rem;font-size:.9rem}th,td{border:1px solid #d1d5db;padding:.35rem;text-align:right}th:first-child,td:first-child{text-align:left}thead{background:#f3f4f6}code{white-space:nowrap}h1,h2{color:#111827}
</style></head><body>
<h1>Metagenomics resource usage</h1>
<p>Cached rows remain in the merged trace but are excluded from consumption totals.</p>
<h2>Summary</h2>
""" + table(("metric", "value"), overview) + """
<h2>Outer workflow scopes</h2>
""" + table(subworkflow_fields, subworkflow_rows) + """
<h2>Process scopes</h2>
""" + table(process_fields, process_rows) + f"""
<h2>Limitations</h2><ul>{limitations}</ul>
</body></html>
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge registered Nextflow traces and summarize resource usage."
    )
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--task-peaks", type=Path)
    parser.add_argument("--task-work-timeseries", type=Path)
    parser.add_argument("--storage-timeseries", type=Path)
    parser.add_argument("--input-manifest", "--sra-manifest", dest="input_manifest", type=Path)
    parser.add_argument(
        "--manifest-mode", choices=("auto", "local", "sra"), default="auto"
    )
    parser.add_argument("--project-status", type=Path)
    parser.add_argument("--gpu-metrics-dir", type=Path)
    parser.add_argument("--slurm-accounting", type=Path)
    return parser


def load_project_status(path: Path | None, context: BuildContext) -> dict[str, Any]:
    if path is None or not path.exists():
        context.limitations.add("Project outcome file was not provided; overall success is unknown.")
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        context.warnings.add(f"Cannot read project outcome file: {exc}")
        return {}
    return value if isinstance(value, dict) else {}


def load_slurm_accounting(path: Path | None, context: BuildContext) -> dict[str, Any]:
    if path is None or not path.exists():
        context.limitations.add("Portable Nextflow telemetry was used without optional sacct enrichment.")
        return {"available": False, "records": 0, "source_path": ""}
    status_path = path.with_suffix(".status.json")
    collection_status: dict[str, Any] = {}
    if status_path.exists():
        try:
            loaded_status = json.loads(status_path.read_text(encoding="utf-8"))
            if isinstance(loaded_status, dict):
                collection_status = loaded_status
        except (OSError, json.JSONDecodeError) as exc:
            context.warnings.add(f"Cannot read SLURM accounting status: {exc}")
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, csv.Error) as exc:
        context.warnings.add(f"Cannot read optional SLURM accounting: {exc}")
        return {"available": False, "records": 0, "source_path": os.path.abspath(path)}
    if collection_status and not bool(collection_status.get("available")):
        context.limitations.add(
            "Optional sacct enrichment did not complete; portable Nextflow telemetry remains authoritative."
        )
        return {
            "available": False,
            "records": 0,
            "state": clean_text(collection_status.get("state")),
            "error": clean_text(collection_status.get("error")),
            "source_path": os.path.abspath(path),
            "status_path": os.path.abspath(status_path),
        }
    if not rows and not collection_status:
        context.limitations.add(
            "The optional sacct file contained no records and no successful collection status."
        )
        return {"available": False, "records": 0, "source_path": os.path.abspath(path)}
    states: dict[str, int] = defaultdict(int)
    gpu_records = 0
    for row in rows:
        states[clean_text(row.get("state")) or "unknown"] += 1
        tres = clean_text(row.get("allocated_tres")).lower()
        gpu_records += int("gpu" in tres)
    return {
        "available": True,
        "records": len(rows),
        "states": dict(sorted(states.items())),
        "records_with_gpu_tres": gpu_records,
        "source_path": os.path.abspath(path),
        "status_path": os.path.abspath(status_path) if status_path.exists() else "",
        "aggregation_rule": "preserved independently; values do not replace Nextflow observations",
    }


def run(args: argparse.Namespace) -> int:
    context = BuildContext()
    invocations = read_registry(args.registry)
    if not invocations:
        raise ValueError("trace registry contains no invocations")
    tasks = load_tasks(invocations, context)
    if not tasks:
        context.limitations.add("No task rows were present in the registered traces.")

    try:
        peaks = load_task_peaks(args.task_peaks)
    except (OSError, csv.Error, ValueError) as exc:
        context.warnings.add(f"Cannot read optional task work-directory peaks: {exc}")
        peaks = {}
    attach_task_peaks(tasks, peaks)
    attach_gpu_metrics(tasks, load_gpu_metrics(args.gpu_metrics_dir, context))
    if args.task_peaks is None or not args.task_peaks.exists():
        context.limitations.add(
            "Task work-directory peaks were not provided; per-task disk metrics are unavailable."
        )

    process_rows = aggregate_by_process(tasks)
    subworkflow_rows = aggregate_by_subworkflow(tasks)
    attach_sampled_concurrent_work(
        args.task_work_timeseries,
        tasks,
        process_rows,
        subworkflow_rows,
        context,
    )
    storage = load_storage(args.storage_timeseries, context)
    attach_stage_storage(subworkflow_rows, storage)
    manifest_counts = load_manifest_counts(
        args.input_manifest, context, args.manifest_mode
    )
    project_status = load_project_status(args.project_status, context)
    slurm_accounting = load_slurm_accounting(args.slurm_accounting, context)
    summary = build_summary(
        invocations,
        tasks,
        process_rows,
        subworkflow_rows,
        storage,
        manifest_counts,
        project_status,
        slurm_accounting,
        context,
    )

    output_dir = args.output_dir
    resource_dir = output_dir / "resources"
    execution_fields = (
        "invocation_id",
        "session_id",
        "invocation_stage",
        "trace_row",
        *[field for field in TASK_FIELDS if field not in {"invocation_id", "session_id", "invocation_stage", "trace_row"}],
        *[f"raw_{field}" for field in sorted(context.raw_columns)],
    )
    execution_rows = []
    for task in tasks:
        row = dict(task.values)
        row.update({f"raw_{key}": task.raw.get(key, "") for key in context.raw_columns})
        execution_rows.append(row)

    write_tsv(output_dir / "execution_trace.tsv", execution_fields, execution_rows)
    write_tsv(resource_dir / "resource_usage_by_task.tsv", TASK_FIELDS, (task.values for task in tasks))
    write_tsv(resource_dir / "resource_usage_by_process.tsv", AGGREGATE_FIELDS, process_rows)
    write_tsv(resource_dir / "resource_usage_by_subworkflow.tsv", AGGREGATE_FIELDS, subworkflow_rows)
    write_tsv(
        resource_dir / "resource_usage_summary.tsv",
        ("metric", "value"),
        flatten_summary(summary),
    )
    atomic_write_text(
        resource_dir / "resource_usage_summary.json",
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
    )
    atomic_write_text(
        resource_dir / "resource_usage_summary.html",
        render_html(summary, process_rows, subworkflow_rows),
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, csv.Error) as exc:
        print(f"resource summarizer error: {exc}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
