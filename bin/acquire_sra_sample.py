#!/usr/bin/env python3
"""Acquire exactly one frozen SRA biological sample with bounded disk usage.

Runs are handled sequentially.  Each ``fasterq-dump`` pair is compressed before
the next run starts, and the compressed run pairs are stream-recompressed into
one gzip stream per mate.  All large temporary directories are explicit; this
program never falls back to the system temporary directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable, Sequence


PROJECT_PATTERN = re.compile(r"^PRJ(?:NA|EB|DB)[0-9]+$")
RUN_PATTERN = re.compile(r"^[SED]RR[0-9]+$")
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


class AcquisitionError(RuntimeError):
    """A user-facing acquisition or manifest error."""


@dataclass(frozen=True)
class Toolchain:
    prefix: tuple[str, ...]
    prefetch: tuple[str, ...]
    validate: tuple[str, ...]
    fasterq: tuple[str, ...]
    pigz: tuple[str, ...]

    def command(self, tool: str, *arguments: str | Path) -> list[str]:
        executable = {
            "prefetch": self.prefetch,
            "vdb-validate": self.validate,
            "fasterq-dump": self.fasterq,
            "pigz": self.pigz,
        }[tool]
        return [*self.prefix, *executable, *(str(value) for value in arguments)]


class CommandRunner:
    def __init__(self, *, dry_run: bool) -> None:
        self.dry_run = dry_run

    @staticmethod
    def display(command: Sequence[str]) -> str:
        return shlex.join(command)

    def run(self, command: Sequence[str]) -> None:
        print(f"+ {self.display(command)}", file=sys.stderr)
        if self.dry_run:
            return
        try:
            subprocess.run(list(command), check=True)
        except FileNotFoundError as exc:
            raise AcquisitionError(f"executable not found: {command[0]}") from exc
        except subprocess.CalledProcessError as exc:
            raise AcquisitionError(
                f"command failed with exit status {exc.returncode}: {self.display(command)}"
            ) from exc


def _command_words(value: str, option_name: str) -> tuple[str, ...]:
    # Offline tests and site wrappers may inject an absolute script path.  On
    # Windows an extensionless shebang script is not directly executable, so
    # run it with the current interpreter.  This branch is never used by the
    # pinned Linux production executables.
    direct_path = Path(value)
    if direct_path.is_file():
        if os.name == "nt":
            try:
                if direct_path.read_bytes()[:2] == b"#!":
                    return (sys.executable, str(direct_path))
            except OSError as exc:
                raise AcquisitionError(f"cannot inspect {option_name}: {exc}") from exc
        return (str(direct_path),)
    try:
        words = tuple(shlex.split(value, posix=os.name != "nt"))
    except ValueError as exc:
        raise AcquisitionError(f"invalid {option_name}: {exc}") from exc
    if not words:
        raise AcquisitionError(f"{option_name} cannot be empty")
    return words


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Acquire and merge exactly one sample from a frozen SRA run manifest."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--scratch-dir",
        required=True,
        type=Path,
        help="Large task-local scratch root (never defaults to /tmp)",
    )
    parser.add_argument(
        "--prefetch-dir",
        type=Path,
        help="Optional prefetch scratch root; a sample-owned subdirectory is created",
    )
    parser.add_argument(
        "--temp-dir",
        type=Path,
        help="Optional fasterq temporary root; a sample-owned subdirectory is created",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--max-size",
        default="u",
        help="Value for prefetch --max-size (default: u, unlimited per selected run)",
    )
    parser.add_argument(
        "--tool-prefix",
        default="",
        help="Command prefix applied separately to every tool, parsed with shell quoting",
    )
    parser.add_argument("--prefetch-executable", default="prefetch")
    parser.add_argument("--vdb-validate-executable", default="vdb-validate")
    parser.add_argument("--fasterq-dump-executable", default="fasterq-dump")
    parser.add_argument("--pigz-executable", default="pigz")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Atomically replace an existing final sample pair",
    )
    args = parser.parse_args(argv)
    if args.threads <= 0:
        parser.error("--threads must be positive")
    if not args.max_size.strip():
        parser.error("--max-size cannot be empty")
    if not SAMPLE_ID_PATTERN.fullmatch(args.sample_id):
        parser.error(
            "--sample-id must begin with an alphanumeric character and contain only "
            "letters, digits, '.', '_', or '-'"
        )
    return args


def load_sample_rows(manifest: Path, sample_id: str) -> list[dict[str, str]]:
    try:
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(RUN_FIELDS):
                raise AcquisitionError(
                    "manifest header does not match the frozen SRA run-manifest schema"
                )
            all_rows = [dict(row) for row in reader]
    except OSError as exc:
        raise AcquisitionError(f"cannot read manifest {manifest}: {exc}") from exc
    except csv.Error as exc:
        raise AcquisitionError(f"cannot parse manifest {manifest}: {exc}") from exc

    rows = [row for row in all_rows if row["sample_id"] == sample_id]
    if not rows:
        available = sorted({row["sample_id"] for row in all_rows})
        suffix = f"; available samples: {', '.join(available[:10])}" if available else ""
        raise AcquisitionError(f"sample {sample_id!r} is absent from the frozen manifest{suffix}")

    seen_runs: set[str] = set()
    sample_orders: set[str] = set()
    parsed: list[tuple[int, dict[str, str]]] = []
    for row in rows:
        run = row["run_accession"]
        if not RUN_PATTERN.fullmatch(run):
            raise AcquisitionError(f"invalid run accession for {sample_id}: {run!r}")
        if run in seen_runs:
            raise AcquisitionError(f"duplicate run accession for {sample_id}: {run}")
        seen_runs.add(run)
        if not PROJECT_PATTERN.fullmatch(row["project_accession"]):
            raise AcquisitionError(f"invalid project accession for run {run}")
        if row["eligibility"] != "eligible" or row["exclusion_reason"]:
            raise AcquisitionError(f"run {run} is not marked eligible")
        expected_values = {
            "layout": "PAIRED",
            "strategy": "WGS",
            "source": "METAGENOMIC",
            "public_access": "true",
        }
        for field_name, expected in expected_values.items():
            if row[field_name] != expected:
                raise AcquisitionError(
                    f"run {run} violates frozen input contract: {field_name}={row[field_name]!r}"
                )
        if not row["platform"]:
            raise AcquisitionError(f"run {run} has no platform")
        if row["spots"] and row["spots_with_mates"]:
            try:
                if int(row["spots"]) != int(row["spots_with_mates"]):
                    raise AcquisitionError(f"run {run} contains spots without mates")
            except ValueError as exc:
                raise AcquisitionError(f"run {run} has invalid spot counts") from exc
        try:
            order = int(row["run_order"])
            if order <= 0:
                raise ValueError
        except ValueError as exc:
            raise AcquisitionError(f"run {run} has invalid run_order") from exc
        parsed.append((order, row))
        sample_orders.add(row["sample_order"])

    if len(sample_orders) != 1:
        raise AcquisitionError(f"sample {sample_id} has inconsistent sample_order values")
    parsed.sort(key=lambda item: item[0])
    if [order for order, _ in parsed] != list(range(1, len(parsed) + 1)):
        raise AcquisitionError(f"sample {sample_id} has non-consecutive run_order values")
    if [row["run_accession"] for _, row in parsed] != sorted(seen_runs):
        raise AcquisitionError(f"sample {sample_id} runs are not in deterministic accession order")
    return [row for _, row in parsed]


def _owned_path(base: Path, sample_id: str, role: str) -> Path:
    return base.expanduser().resolve() / f"{sample_id}.{role}"


def _prepare_owned_directory(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.mkdir()
    except FileExistsError as exc:
        raise AcquisitionError(
            f"task-local directory already exists; refusing to reuse or delete it: {path}"
        ) from exc


def _remove_owned_directory(path: Path, sample_id: str, role: str) -> None:
    if path.name != f"{sample_id}.{role}":
        raise AcquisitionError(f"internal safety check refused cleanup of {path}")
    if path.exists():
        shutil.rmtree(path)


def _locate_prefetched_accession(prefetch_root: Path, run: str) -> Path:
    """Return the accession-as-directory produced by ``prefetch``.

    NCBI explicitly recommends passing this directory to ``fasterq-dump``
    because a run can contain more than the top-level ``.sra`` file.  The
    recursive fallback tolerates older/fake ``prefetch`` layouts while still
    returning the directory whose basename is the accession.
    """

    run_dir = prefetch_root / run
    expected = run_dir / f"{run}.sra"
    if expected.is_file() and run_dir.name == run:
        return run_dir
    matches = sorted(run_dir.rglob(f"{run}.sra")) if run_dir.is_dir() else []
    if len(matches) == 1 and matches[0].is_file() and matches[0].parent.name == run:
        return matches[0].parent
    raise AcquisitionError(
        f"prefetch completed but did not create exactly one accession directory for {run} below {run_dir}"
    )


def _assert_fasterq_pair(output_dir: Path, run: str) -> tuple[Path, Path]:
    expected = (output_dir / f"{run}_1.fastq", output_dir / f"{run}_2.fastq")
    entries = sorted(path for path in output_dir.iterdir())
    unexpected = [path.name for path in entries if path not in expected]
    missing = [path.name for path in expected if not path.is_file()]
    empty = [path.name for path in expected if path.is_file() and path.stat().st_size == 0]
    if unexpected or missing or empty:
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if empty:
            details.append("empty=" + ",".join(empty))
        if unexpected:
            details.append("unexpected/singleton=" + ",".join(unexpected))
        raise AcquisitionError(
            f"fasterq-dump output for {run} is not exactly one non-empty paired FASTQ set "
            f"({' ; '.join(details)})"
        )
    return expected


def _compress_mate(
    source: Path,
    *,
    toolchain: Toolchain,
    runner: CommandRunner,
    threads: int,
) -> Path:
    compressed = Path(f"{source}.gz")
    runner.run(toolchain.command("pigz", "-p", str(threads), source))
    if runner.dry_run:
        return compressed
    if source.exists():
        raise AcquisitionError(f"pigz did not remove uncompressed FASTQ: {source}")
    if not compressed.is_file() or compressed.stat().st_size == 0:
        raise AcquisitionError(f"pigz did not create a non-empty gzip file: {compressed}")
    return compressed


def _verify_nonempty_gzip(path: Path) -> None:
    try:
        with gzip.open(path, "rb") as handle:
            first = handle.read(1024 * 1024)
            if not first:
                raise AcquisitionError(f"merged FASTQ is empty: {path}")
            while handle.read(1024 * 1024):
                pass
    except (OSError, EOFError) as exc:
        raise AcquisitionError(f"merged FASTQ is not a readable gzip stream: {path}") from exc


def _paired_read_name(header: bytes, mate: int) -> bytes:
    if not header.startswith(b"@"):
        raise AcquisitionError(f"mate {mate} FASTQ header does not begin with '@'")
    fields = header[1:].strip().split()
    if not fields:
        raise AcquisitionError(f"mate {mate} FASTQ has an empty read name")
    name = fields[0]
    expected_suffix = f"/{mate}".encode()
    other_suffix = f"/{3 - mate}".encode()
    if name.endswith(other_suffix):
        raise AcquisitionError(f"mate {mate} FASTQ contains a mate-{3 - mate} read name")
    if name.endswith(expected_suffix):
        name = name[: -len(expected_suffix)]
    if len(fields) > 1 and fields[1][:2] in (b"1:", b"2:"):
        expected_token = f"{mate}:".encode()
        if fields[1][:2] != expected_token:
            raise AcquisitionError(f"mate {mate} FASTQ has an inconsistent CASAVA mate token")
    return name


def _next_fastq_record(handle: BinaryIO, path: Path, mate: int) -> tuple[bytes, bytes, bytes, bytes] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise AcquisitionError(f"truncated FASTQ record in {path}")
    if not plus.startswith(b"+"):
        raise AcquisitionError(f"FASTQ separator does not begin with '+' in {path}")
    if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
        raise AcquisitionError(f"sequence/quality length mismatch in {path}")
    _paired_read_name(header, mate)
    return header, sequence, plus, quality


def _validate_fastq_pair(read_1: Path, read_2: Path) -> int:
    """Fully validate gzip integrity, FASTQ structure, pairing, and record count."""

    count = 0
    try:
        with gzip.open(read_1, "rb") as first, gzip.open(read_2, "rb") as second:
            while True:
                record_1 = _next_fastq_record(first, read_1, 1)
                record_2 = _next_fastq_record(second, read_2, 2)
                if record_1 is None or record_2 is None:
                    if record_1 is not None or record_2 is not None:
                        raise AcquisitionError("merged FASTQ mates contain different record counts")
                    break
                if _paired_read_name(record_1[0], 1) != _paired_read_name(record_2[0], 2):
                    raise AcquisitionError(f"merged FASTQ mate names differ at record {count + 1}")
                count += 1
    except (OSError, EOFError) as exc:
        raise AcquisitionError(f"merged FASTQ pair is not valid gzip data: {exc}") from exc
    if count == 0:
        raise AcquisitionError("merged FASTQ pair contains no records")
    return count


def _stream_recompress(
    inputs: Sequence[Path],
    output: Path,
    *,
    toolchain: Toolchain,
    runner: CommandRunner,
    threads: int,
) -> None:
    decompress = toolchain.command("pigz", "-d", "-c", *inputs)
    compress = toolchain.command("pigz", "-p", str(threads), "-c")
    print(
        f"+ {runner.display(decompress)} | {runner.display(compress)} > {shlex.quote(str(output))}",
        file=sys.stderr,
    )
    if runner.dry_run:
        return
    try:
        with output.open("xb") as output_handle:
            left = subprocess.Popen(decompress, stdout=subprocess.PIPE)
            assert left.stdout is not None
            try:
                right = subprocess.Popen(compress, stdin=left.stdout, stdout=output_handle)
            except BaseException:
                left.kill()
                left.wait()
                raise
            left.stdout.close()
            right_status = right.wait()
            left_status = left.wait()
    except FileNotFoundError as exc:
        raise AcquisitionError(f"pigz executable not found: {exc.filename}") from exc
    except OSError as exc:
        raise AcquisitionError(f"cannot create merged FASTQ {output}: {exc}") from exc
    if left_status != 0 or right_status != 0:
        raise AcquisitionError(
            "stream recompression failed: "
            f"decompressor exit={left_status}, compressor exit={right_status}"
        )
    if output.stat().st_size == 0:
        raise AcquisitionError(f"stream recompression created an empty file: {output}")
    _verify_nonempty_gzip(output)


def _unlink_uncompressed_fastq(paths: Iterable[Path]) -> None:
    for directory in paths:
        if directory.is_dir():
            for fastq in directory.rglob("*.fastq"):
                try:
                    fastq.unlink()
                except FileNotFoundError:
                    pass


def build_toolchain(args: argparse.Namespace) -> Toolchain:
    prefix = _command_words(args.tool_prefix, "--tool-prefix") if args.tool_prefix else ()
    return Toolchain(
        prefix=prefix,
        prefetch=_command_words(args.prefetch_executable, "--prefetch-executable"),
        validate=_command_words(args.vdb_validate_executable, "--vdb-validate-executable"),
        fasterq=_command_words(args.fasterq_dump_executable, "--fasterq-dump-executable"),
        pigz=_command_words(args.pigz_executable, "--pigz-executable"),
    )


def dry_run_plan(
    rows: Sequence[dict[str, str]],
    *,
    args: argparse.Namespace,
    toolchain: Toolchain,
) -> None:
    runner = CommandRunner(dry_run=True)
    scratch = _owned_path(args.scratch_dir, args.sample_id, "acquisition")
    prefetch = _owned_path(args.prefetch_dir or args.scratch_dir, args.sample_id, "prefetch")
    temporary = _owned_path(args.temp_dir or args.scratch_dir, args.sample_id, "fasterq-temp")
    pairs: list[tuple[Path, Path]] = []
    for row in rows:
        run = row["run_accession"]
        fastq_dir = scratch / "runs" / run
        accession_dir = prefetch / run
        runner.run(
            toolchain.command(
                "prefetch", "--output-directory", prefetch, "--max-size", args.max_size, run
            )
        )
        runner.run(toolchain.command("vdb-validate", accession_dir))
        runner.run(
            toolchain.command(
                "fasterq-dump",
                "--split-files",
                "--skip-technical",
                "--threads",
                str(args.threads),
                "--temp",
                temporary / run,
                "--outdir",
                fastq_dir,
                accession_dir,
            )
        )
        pair = (
            _compress_mate(
                fastq_dir / f"{run}_1.fastq",
                toolchain=toolchain,
                runner=runner,
                threads=args.threads,
            ),
            _compress_mate(
                fastq_dir / f"{run}_2.fastq",
                toolchain=toolchain,
                runner=runner,
                threads=args.threads,
            ),
        )
        pairs.append(pair)
    output_1 = args.output_dir / f"{args.sample_id}_R1.fastq.gz"
    output_2 = args.output_dir / f"{args.sample_id}_R2.fastq.gz"
    _stream_recompress(
        [pair[0] for pair in pairs],
        output_1,
        toolchain=toolchain,
        runner=runner,
        threads=args.threads,
    )
    _stream_recompress(
        [pair[1] for pair in pairs],
        output_2,
        toolchain=toolchain,
        runner=runner,
        threads=args.threads,
    )


def acquire(
    rows: Sequence[dict[str, str]],
    *,
    args: argparse.Namespace,
    toolchain: Toolchain,
) -> tuple[Path, Path]:
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    final_1 = output_dir / f"{args.sample_id}_R1.fastq.gz"
    final_2 = output_dir / f"{args.sample_id}_R2.fastq.gz"
    existing = [path for path in (final_1, final_2) if path.exists()]
    if existing and not args.force:
        raise AcquisitionError(
            "refusing to replace existing output(s): " + ", ".join(str(path) for path in existing)
        )

    scratch = _owned_path(args.scratch_dir, args.sample_id, "acquisition")
    prefetch = _owned_path(args.prefetch_dir or args.scratch_dir, args.sample_id, "prefetch")
    temporary = _owned_path(args.temp_dir or args.scratch_dir, args.sample_id, "fasterq-temp")
    owned = [(scratch, "acquisition"), (prefetch, "prefetch"), (temporary, "fasterq-temp")]
    if len({path for path, _ in owned}) != len(owned):
        raise AcquisitionError("scratch, prefetch, and fasterq temporary paths must be distinct")

    partial_1 = output_dir / f".{args.sample_id}_R1.fastq.gz.partial.{os.getpid()}"
    partial_2 = output_dir / f".{args.sample_id}_R2.fastq.gz.partial.{os.getpid()}"
    created: list[tuple[Path, str]] = []
    run_pairs: list[tuple[Path, Path]] = []
    runner = CommandRunner(dry_run=False)
    try:
        for path, role in owned:
            _prepare_owned_directory(path)
            created.append((path, role))
        runs_root = scratch / "runs"
        runs_root.mkdir()

        for row in rows:
            run = row["run_accession"]
            print(f"Acquiring {run} for biological sample {args.sample_id}", file=sys.stderr)
            run_fastq_dir = runs_root / run
            run_fastq_dir.mkdir()
            run_temp_dir = temporary / run
            run_temp_dir.mkdir()
            runner.run(
                toolchain.command(
                    "prefetch",
                    "--output-directory",
                    prefetch,
                    "--max-size",
                    args.max_size,
                    run,
                )
            )
            accession_dir = _locate_prefetched_accession(prefetch, run)
            runner.run(toolchain.command("vdb-validate", accession_dir))
            runner.run(
                toolchain.command(
                    "fasterq-dump",
                    "--split-files",
                    "--skip-technical",
                    "--threads",
                    str(args.threads),
                    "--temp",
                    run_temp_dir,
                    "--outdir",
                    run_fastq_dir,
                    accession_dir,
                )
            )
            mate_1, mate_2 = _assert_fasterq_pair(run_fastq_dir, run)
            compressed_1 = _compress_mate(
                mate_1, toolchain=toolchain, runner=runner, threads=args.threads
            )
            compressed_2 = _compress_mate(
                mate_2, toolchain=toolchain, runner=runner, threads=args.threads
            )
            run_pairs.append((compressed_1, compressed_2))

            # Once both run mates are compressed, the SRA object and fasterq
            # temporary data are no longer needed.  Removing them here bounds
            # acquisition storage to one SRA/uncompressed run at a time.
            run_prefetch_dir = prefetch / run
            run_temp_dir = temporary / run
            if run_prefetch_dir.is_dir():
                shutil.rmtree(run_prefetch_dir)
            if run_temp_dir.is_dir():
                shutil.rmtree(run_temp_dir)

        _stream_recompress(
            [pair[0] for pair in run_pairs],
            partial_1,
            toolchain=toolchain,
            runner=runner,
            threads=args.threads,
        )
        _stream_recompress(
            [pair[1] for pair in run_pairs],
            partial_2,
            toolchain=toolchain,
            runner=runner,
            threads=args.threads,
        )
        _validate_fastq_pair(partial_1, partial_2)
        os.replace(partial_1, final_1)
        os.replace(partial_2, final_2)
        return final_1, final_2
    finally:
        _unlink_uncompressed_fastq(path for path, _ in created)
        for partial in (partial_1, partial_2):
            try:
                partial.unlink()
            except FileNotFoundError:
                pass
        cleanup_errors: list[str] = []
        for path, role in reversed(created):
            try:
                _remove_owned_directory(path, args.sample_id, role)
            except (OSError, AcquisitionError) as exc:
                cleanup_errors.append(str(exc))
        if cleanup_errors:
            print("WARNING: task-local cleanup incomplete: " + "; ".join(cleanup_errors), file=sys.stderr)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        rows = load_sample_rows(args.manifest, args.sample_id)
        toolchain = build_toolchain(args)
        if args.dry_run:
            dry_run_plan(rows, args=args, toolchain=toolchain)
            return 0
        output_1, output_2 = acquire(rows, args=args, toolchain=toolchain)
        print(f"Created {output_1} and {output_2}", file=sys.stderr)
        return 0
    except AcquisitionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
