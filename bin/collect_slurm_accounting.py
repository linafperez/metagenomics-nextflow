#!/usr/bin/env python3
"""Best-effort preservation of sacct records for registered Nextflow tasks."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


JOB_PATTERN = re.compile(r"^[0-9]+(?:_[0-9]+)?(?:\.[A-Za-z0-9_-]+)?$")
FIELDS = (
    "invocation_id",
    "native_id",
    "job_id_raw",
    "state",
    "elapsed_raw_seconds",
    "total_cpu",
    "max_rss",
    "allocated_cpus",
    "allocated_tres",
    "requested_tres",
)


def atomic_tsv(path: Path, rows: Iterable[dict[str, str]]) -> None:
    from io import StringIO

    stream = StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(stream.getvalue())
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def atomic_json(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def registered_jobs(registry: Path) -> dict[str, set[str]]:
    jobs: dict[str, set[str]] = {}
    with registry.open("r", encoding="utf-8", newline="") as handle:
        for invocation in csv.DictReader(handle, delimiter="\t"):
            invocation_id = (invocation.get("invocation_id") or "").strip()
            trace_path = Path((invocation.get("trace_path") or "").strip())
            if not trace_path.is_absolute():
                trace_path = registry.parent / trace_path
            if not trace_path.exists():
                continue
            with trace_path.open("r", encoding="utf-8", newline="") as trace_handle:
                for row in csv.DictReader(trace_handle, delimiter="\t"):
                    native_id = (row.get("native_id") or "").strip()
                    if JOB_PATTERN.fullmatch(native_id):
                        jobs.setdefault(native_id, set()).add(invocation_id)
    return jobs


def chunks(values: Sequence[str], size: int = 100) -> Iterable[Sequence[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def natural_job_key(value: str) -> tuple[tuple[int, object], ...]:
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part)
        for part in re.split(r"([0-9]+)", value)
    )


def collect(
    sacct: str, jobs: dict[str, set[str]], timeout_seconds: float
) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    requested = sorted(jobs, key=natural_job_key)
    deadline = time.monotonic() + timeout_seconds
    for batch in chunks(requested):
        command = [
            sacct,
            "--noheader",
            "--parsable2",
            "--jobs",
            ",".join(batch),
            "--format=JobIDRaw,State,ElapsedRaw,TotalCPU,MaxRSS,AllocCPUS,AllocTRES,ReqTRES",
        ]
        remaining_seconds = deadline - time.monotonic()
        if remaining_seconds <= 0:
            raise subprocess.TimeoutExpired(sacct, timeout_seconds)
        result = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            timeout=remaining_seconds,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or f"sacct exited {result.returncode}")
        for line in result.stdout.splitlines():
            values = line.split("|")
            if len(values) < 8:
                continue
            job_id = values[0].strip()
            base_id = job_id.split(".", 1)[0]
            invocation_ids = sorted(jobs.get(job_id, jobs.get(base_id, {""})))
            for invocation_id in invocation_ids:
                records.append(
                    dict(
                        zip(
                            FIELDS,
                            (
                                invocation_id,
                                base_id,
                                job_id,
                                values[1].strip(),
                                values[2].strip(),
                                values[3].strip(),
                                values[4].strip(),
                                values[5].strip(),
                                values[6].strip(),
                                values[7].strip(),
                            ),
                        )
                    )
                )
    return sorted(records, key=lambda row: (row["invocation_id"], row["native_id"], row["job_id_raw"]))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sacct", default="sacct")
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--status-output", type=Path)
    args = parser.parse_args(argv)
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be greater than zero")
    status_output = args.status_output or args.output.with_suffix(".status.json")
    error = ""
    try:
        jobs = registered_jobs(args.registry)
        rows = collect(args.sacct, jobs, args.timeout_seconds) if jobs else []
        available = bool(jobs)
        state = "collected" if jobs else "no_registered_slurm_jobs"
    except (OSError, csv.Error, RuntimeError, subprocess.TimeoutExpired) as exc:
        print(f"SLURM accounting warning: {exc}", file=sys.stderr)
        rows = []
        available = False
        state = "error"
        error = f"{exc.__class__.__name__}: {exc}"
    atomic_tsv(args.output, rows)
    atomic_json(
        status_output,
        {
            "schema_version": 1,
            "available": available,
            "state": state,
            "records": len(rows),
            "error": error,
            "collected_at_utc": datetime.now(timezone.utc)
            .isoformat(timespec="seconds")
            .replace("+00:00", "Z"),
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
