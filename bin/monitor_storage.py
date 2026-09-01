#!/usr/bin/env python3
"""Best-effort storage telemetry for the metagenomics workflow.

The monitor deliberately uses only the Python standard library.  It measures
allocated filesystem bytes (the same quantity that ``du`` approximates), does
not follow symbolic links, and continuously checkpoints the largest observed
size of every Nextflow task work directory.

Runtime measurement errors are recorded in the timeseries and reported on
stderr, but result in a successful exit so that optional telemetry cannot turn
a scientific workflow into a failed workflow.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import signal
import stat
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


TIMESERIES_FIELDS = (
    "timestamp",
    "invocation_id",
    "invocation_stage",
    "work_bytes",
    "checkpoint_bytes",
    "sra_cache_bytes",
    "sra_scratch_bytes",
    "sra_temp_bytes",
    "results_bytes",
    "database_bytes",
    "total_dynamic_bytes",
    "total_measured_bytes",
    "measurement_complete",
    "errors",
)

TASK_TIMESERIES_FIELDS = (
    "timestamp",
    "workdir",
    "work_bytes",
)

TASK_PEAK_FIELDS = (
    "workdir",
    "peak_work_bytes",
    "first_observed",
    "last_observed",
    "samples",
)

WORK_PREFIX_RE = re.compile(r"^[0-9a-fA-F]{2}$")


@dataclass
class ScanResult:
    allocated_bytes: int | None
    errors: list[str] = field(default_factory=list)
    task_bytes: dict[str, int] = field(default_factory=dict)


@dataclass
class PeakRecord:
    peak_work_bytes: int
    first_observed: str
    last_observed: str
    samples: int


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z"
    )


def allocated_bytes(info: os.stat_result) -> int:
    """Return allocated bytes, falling back to logical size where unavailable."""

    blocks = getattr(info, "st_blocks", None)
    # CPython exposes ``st_blocks == 0`` for ordinary non-empty files on some
    # Windows filesystems even though allocated-block accounting is unavailable.
    # Treat that platform sentinel as unavailable instead of reporting every
    # task and project tree as zero bytes.
    if blocks is not None and not (os.name == "nt" and blocks == 0 and info.st_size > 0):
        return int(blocks) * 512
    return int(info.st_size)


def _display_path(path: Path) -> str:
    return os.path.abspath(os.fspath(path))


def file_identity(info: os.stat_result, path: Path) -> tuple[object, ...]:
    """Return a hard-link identity, falling back to the path when unavailable."""

    inode = int(getattr(info, "st_ino", 0) or 0)
    device = int(getattr(info, "st_dev", 0) or 0)
    if inode:
        return ("inode", device, inode)
    return ("path", os.path.normcase(_display_path(path)))


def scan_tree(root: Path, *, identify_tasks: bool = False) -> ScanResult:
    """Measure *root* without following symlinks.

    When ``identify_tasks`` is true, two-level directories below a hexadecimal
    Nextflow work prefix are treated as task work directories.  Global totals
    deduplicate hard-linked inodes.  Per-task totals deduplicate hard links
    independently within each task, which avoids arbitrary attribution between
    otherwise independent task directories.
    """

    root = Path(os.path.abspath(os.fspath(root)))
    errors: list[str] = []
    task_totals: dict[str, int] = {}
    seen_global: set[tuple[object, ...]] = set()
    seen_tasks: dict[str, set[tuple[object, ...]]] = {}

    try:
        root_info = os.lstat(root)
    except FileNotFoundError:
        return ScanResult(None, [f"missing:{root}"])
    except OSError as exc:
        return ScanResult(None, [f"stat:{root}:{exc.__class__.__name__}"])

    root_key = file_identity(root_info, root)
    total = allocated_bytes(root_info)
    seen_global.add(root_key)

    if stat.S_ISLNK(root_info.st_mode):
        return ScanResult(total, [f"root_is_symlink:{root}"])
    if not stat.S_ISDIR(root_info.st_mode):
        return ScanResult(total)

    # (directory, depth relative to root, current task directory)
    stack: list[tuple[Path, int, str | None]] = [(root, 0, None)]
    while stack:
        directory, depth, current_task = stack.pop()
        try:
            with os.scandir(directory) as iterator:
                entries = sorted(iterator, key=lambda item: item.name, reverse=True)
        except OSError as exc:
            errors.append(f"scan:{directory}:{exc.__class__.__name__}")
            continue

        for entry in entries:
            entry_path = Path(entry.path)
            try:
                info = entry.stat(follow_symlinks=False)
            except OSError as exc:
                errors.append(f"stat:{entry_path}:{exc.__class__.__name__}")
                continue

            mode = info.st_mode
            inode_key = file_identity(info, entry_path)
            entry_task = current_task
            if (
                identify_tasks
                and current_task is None
                and WORK_PREFIX_RE.fullmatch(directory.name)
                and stat.S_ISDIR(mode)
                and not stat.S_ISLNK(mode)
            ):
                entry_task = _display_path(entry_path)
                task_totals.setdefault(entry_task, 0)
                seen_tasks.setdefault(entry_task, set())

            size = allocated_bytes(info)
            if inode_key not in seen_global:
                total += size
                seen_global.add(inode_key)

            if entry_task is not None:
                task_seen = seen_tasks.setdefault(entry_task, set())
                if inode_key not in task_seen:
                    task_totals[entry_task] = task_totals.get(entry_task, 0) + size
                    task_seen.add(inode_key)

            if stat.S_ISDIR(mode) and not stat.S_ISLNK(mode):
                stack.append((entry_path, depth + 1, entry_task))

    return ScanResult(total, errors, task_totals)


def load_task_peaks(path: Path) -> dict[str, PeakRecord]:
    if not path.exists():
        return {}
    peaks: dict[str, PeakRecord] = {}
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                workdir = row.get("workdir", "").strip()
                if not workdir:
                    continue
                try:
                    record = PeakRecord(
                        peak_work_bytes=int(row.get("peak_work_bytes", "0")),
                        first_observed=row.get("first_observed", ""),
                        last_observed=row.get("last_observed", ""),
                        samples=int(row.get("samples", "0")),
                    )
                except (TypeError, ValueError):
                    continue
                previous = peaks.get(workdir)
                if previous is None or record.peak_work_bytes > previous.peak_work_bytes:
                    peaks[workdir] = record
    except OSError:
        return {}
    return peaks


def atomic_write_tsv(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def persist_task_peaks(path: Path, peaks: Mapping[str, PeakRecord]) -> None:
    rows = (
        {
            "workdir": workdir,
            "peak_work_bytes": record.peak_work_bytes,
            "first_observed": record.first_observed,
            "last_observed": record.last_observed,
            "samples": record.samples,
        }
        for workdir, record in sorted(peaks.items())
    )
    atomic_write_tsv(path, TASK_PEAK_FIELDS, rows)


def append_timeseries(path: Path, row: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    needs_header = not path.exists() or path.stat().st_size == 0
    if not needs_header:
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                existing_header = handle.readline().rstrip("\r\n").split("\t")
            if existing_header != list(TIMESERIES_FIELDS):
                raise ValueError("existing timeseries has an incompatible header")
        except OSError as exc:
            raise RuntimeError(f"cannot inspect {path}: {exc}") from exc
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TIMESERIES_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        if needs_header:
            writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def append_task_timeseries(
    path: Path, timestamp: str, task_bytes: Mapping[str, int]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    needs_header = not path.exists() or path.stat().st_size == 0
    if not needs_header:
        with path.open("r", encoding="utf-8", newline="") as handle:
            existing_header = handle.readline().rstrip("\r\n").split("\t")
        if existing_header != list(TASK_TIMESERIES_FIELDS):
            raise ValueError("existing task timeseries has an incompatible header")
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=TASK_TIMESERIES_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        if needs_header:
            writer.writeheader()
        writer.writerows(
            {
                "timestamp": timestamp,
                "workdir": workdir,
                "work_bytes": size,
            }
            for workdir, size in sorted(task_bytes.items())
        )
        handle.flush()
        os.fsync(handle.fileno())


def overlapping_path_warnings(paths: Mapping[str, Path | None]) -> list[str]:
    configured = [
        (name, Path(os.path.abspath(os.fspath(path))))
        for name, path in paths.items()
        if path is not None
    ]
    warnings: list[str] = []
    for index, (left_name, left_path) in enumerate(configured):
        for right_name, right_path in configured[index + 1 :]:
            try:
                common = Path(os.path.commonpath((left_path, right_path)))
            except ValueError:
                continue
            if common == left_path or common == right_path:
                warnings.append(f"overlapping_roots:{left_name}:{right_name}")
    return warnings


class StorageMonitor:
    def __init__(
        self,
        *,
        paths: Mapping[str, Path | None],
        timeseries_path: Path,
        task_peaks_path: Path,
        task_timeseries_path: Path | None = None,
        stage_file: Path | None = None,
    ) -> None:
        self.paths = dict(paths)
        self.timeseries_path = timeseries_path
        self.task_peaks_path = task_peaks_path
        self.task_timeseries_path = (
            task_timeseries_path
            if task_timeseries_path is not None
            else timeseries_path.parent / "task_workdir_timeseries.tsv"
        )
        self.stage_file = stage_file
        self.peaks = load_task_peaks(task_peaks_path)
        self.static_database = ScanResult(None)
        self._database_thread: threading.Thread | None = None
        self.configuration_warnings = overlapping_path_warnings(self.paths)

    def measure_database(self) -> None:
        path = self.paths.get("database")
        if path is not None:
            self.static_database = scan_tree(path)

    def start_database_measurement(self) -> None:
        if self.paths.get("database") is None or self._database_thread is not None:
            return
        self._database_thread = threading.Thread(
            target=self.measure_database,
            name="database-storage-measurement",
            daemon=True,
        )
        self._database_thread.start()

    def sample(self) -> None:
        timestamp = utc_now()
        invocation_id = ""
        invocation_stage = ""
        stage_warning = ""
        if self.stage_file is not None:
            try:
                stage_fields = self.stage_file.read_text(encoding="utf-8").rstrip(
                    "\r\n"
                ).split("\t", 1)
                invocation_id = stage_fields[0].strip()
                invocation_stage = (
                    stage_fields[1].strip() if len(stage_fields) > 1 else ""
                )
            except FileNotFoundError:
                pass
            except OSError as exc:
                stage_warning = f"stage_context:{exc.__class__.__name__}"
        scans: dict[str, ScanResult] = {}
        for category in (
            "work",
            "checkpoint",
            "sra_cache",
            "sra_scratch",
            "sra_temp",
            "results",
        ):
            path = self.paths.get(category)
            if path is None:
                scans[category] = ScanResult(None)
            else:
                scans[category] = scan_tree(path, identify_tasks=category == "work")
        scans["database"] = self.static_database

        for workdir, size in scans["work"].task_bytes.items():
            previous = self.peaks.get(workdir)
            if previous is None:
                self.peaks[workdir] = PeakRecord(size, timestamp, timestamp, 1)
            else:
                previous.peak_work_bytes = max(previous.peak_work_bytes, size)
                previous.last_observed = timestamp
                previous.samples += 1
        persist_task_peaks(self.task_peaks_path, self.peaks)
        append_task_timeseries(
            self.task_timeseries_path, timestamp, scans["work"].task_bytes
        )

        errors = list(self.configuration_warnings)
        if stage_warning:
            errors.append(stage_warning)
        for category in (
            "work",
            "checkpoint",
            "sra_cache",
            "sra_scratch",
            "sra_temp",
            "results",
            "database",
        ):
            errors.extend(f"{category}:{item}" for item in scans[category].errors)

        dynamic_categories = (
            "work",
            "checkpoint",
            "sra_cache",
            "sra_scratch",
            "sra_temp",
            "results",
        )
        dynamic_measured = [
            scans[category].allocated_bytes
            for category in dynamic_categories
            if scans[category].allocated_bytes is not None
        ]
        total_dynamic = sum(dynamic_measured)
        database_bytes = scans["database"].allocated_bytes
        row: dict[str, object] = {
            "timestamp": timestamp,
            "invocation_id": invocation_id,
            "invocation_stage": invocation_stage,
            "work_bytes": _optional_number(scans["work"].allocated_bytes),
            "checkpoint_bytes": _optional_number(scans["checkpoint"].allocated_bytes),
            "sra_cache_bytes": _optional_number(scans["sra_cache"].allocated_bytes),
            "sra_scratch_bytes": _optional_number(
                scans["sra_scratch"].allocated_bytes
            ),
            "sra_temp_bytes": _optional_number(scans["sra_temp"].allocated_bytes),
            "results_bytes": _optional_number(scans["results"].allocated_bytes),
            "database_bytes": _optional_number(database_bytes),
            "total_dynamic_bytes": total_dynamic,
            "total_measured_bytes": total_dynamic + (database_bytes or 0),
            "measurement_complete": "true" if not errors else "false",
            "errors": ";".join(sorted(set(errors))),
        }
        append_timeseries(self.timeseries_path, row)
        if errors:
            print(
                "storage monitor warning: " + "; ".join(sorted(set(errors))),
                file=sys.stderr,
                flush=True,
            )


def _optional_number(value: int | None) -> int | str:
    return "" if value is None else value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Sample workflow storage without following symlinks. Runtime "
            "measurement errors are best-effort and do not produce failure exits."
        )
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--sra-cache-dir", type=Path)
    parser.add_argument("--sra-scratch-dir", type=Path)
    parser.add_argument("--sra-temp-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--database-dir", type=Path)
    parser.add_argument("--timeseries", type=Path)
    parser.add_argument("--task-peaks", type=Path)
    parser.add_argument("--task-timeseries", type=Path)
    parser.add_argument("--stage-file", type=Path)
    parser.add_argument("--stop-file", type=Path)
    parser.add_argument(
        "--sample-request-file",
        type=Path,
        help="Request/acknowledgement file used to force a sample before cleanup.",
    )
    parser.add_argument("--interval-seconds", type=float, default=60.0)
    parser.add_argument(
        "--once", action="store_true", help="Take exactly one sample and exit."
    )
    return parser


def run(args: argparse.Namespace) -> int:
    if args.interval_seconds <= 0:
        raise ValueError("--interval-seconds must be greater than zero")

    output_dir = Path(args.output_dir)
    timeseries = args.timeseries or output_dir / "storage_usage_timeseries.tsv"
    task_peaks = args.task_peaks or output_dir / "task_workdir_peaks.tsv"
    task_timeseries = (
        args.task_timeseries or output_dir / "task_workdir_timeseries.tsv"
    )
    paths = {
        "work": args.work_dir,
        "checkpoint": args.checkpoint_dir,
        "sra_cache": args.sra_cache_dir,
        "sra_scratch": args.sra_scratch_dir,
        "sra_temp": args.sra_temp_dir,
        "results": args.results_dir,
        "database": args.database_dir,
    }
    monitor = StorageMonitor(
        paths=paths,
        timeseries_path=timeseries,
        task_peaks_path=task_peaks,
        task_timeseries_path=task_timeseries,
        stage_file=args.stage_file,
    )

    if args.once:
        monitor.measure_database()
        monitor.sample()
        return 0

    stop_event = threading.Event()

    def request_stop(_signum: int, _frame: object) -> None:
        stop_event.set()

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)

    def request_token() -> str | None:
        if args.sample_request_file is None:
            return None
        try:
            return args.sample_request_file.read_text(encoding="utf-8")
        except FileNotFoundError:
            return None
        except OSError as exc:
            print(f"storage monitor warning: cannot read sample request: {exc}", file=sys.stderr)
            return None

    def acknowledge_request(token: str | None) -> None:
        if token is None or args.sample_request_file is None:
            return
        try:
            if args.sample_request_file.read_text(encoding="utf-8") == token:
                args.sample_request_file.unlink()
        except FileNotFoundError:
            pass
        except OSError as exc:
            print(f"storage monitor warning: cannot acknowledge sample request: {exc}", file=sys.stderr)

    # Persist a dynamic snapshot before the potentially expensive one-time
    # database scan starts.  The database measurement runs in the background.
    initial_token = request_token()
    monitor.sample()
    acknowledge_request(initial_token)
    monitor.start_database_measurement()

    while not stop_event.is_set():
        if args.stop_file is not None and args.stop_file.exists():
            break
        deadline = time.monotonic() + args.interval_seconds
        while not stop_event.is_set() and time.monotonic() < deadline:
            if args.stop_file is not None and args.stop_file.exists():
                stop_event.set()
                break
            token = request_token()
            if token is not None:
                monitor.sample()
                acknowledge_request(token)
                deadline = time.monotonic() + args.interval_seconds
            stop_event.wait(min(1.0, max(0.0, deadline - time.monotonic())))
        if not stop_event.is_set() and time.monotonic() >= deadline:
            monitor.sample()
    # Capture the post-cleanup/failure state instead of ending with a stale
    # pre-stop sample.  Telemetry exceptions remain non-fatal in main().
    final_token = request_token()
    monitor.sample()
    acknowledge_request(final_token)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as exc:  # telemetry must never determine scientific success
        print(
            f"storage monitor warning: telemetry stopped after {exc.__class__.__name__}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
