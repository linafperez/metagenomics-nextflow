#!/usr/bin/env python3
"""Launch the tests-only assembler and binner benchmark matrix."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ASSEMBLERS = ("megahit", "spades", "both")
BINNERS = ("comebin", "metabat2", "semibin2", "vamb", "all")
PROJECT_DIR = Path(__file__).resolve().parents[2]


def variant_names() -> list[str]:
    return [f"{assembler}_{binner}" for assembler in ASSEMBLERS for binner in BINNERS]


def split_variant(name: str) -> tuple[str, str]:
    assembler, binner = name.split("_", maxsplit=1)
    if assembler not in ASSEMBLERS or binner not in BINNERS:
        raise ValueError(f"Unsupported benchmark variant: {name}")
    return assembler, binner


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one or all of the 15 metagenomics benchmark variants."
    )
    environment = parser.add_mutually_exclusive_group(required=True)
    environment.add_argument("--local", action="store_const", const="local", dest="environment")
    environment.add_argument("--hpc", action="store_const", const="hpc", dest="environment")

    runtime = parser.add_mutually_exclusive_group()
    runtime.add_argument("--docker", action="store_const", const="docker", dest="runtime")
    runtime.add_argument("--conda", action="store_const", const="conda", dest="runtime")
    runtime.add_argument("--apptainer", action="store_const", const="apptainer", dest="runtime")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--stub", action="store_const", const="stub", dest="mode")
    mode.add_argument("--run", action="store_const", const="run", dest="mode")

    parser.add_argument(
        "--variant",
        action="append",
        choices=variant_names(),
        help="Run only this variant; repeat to select multiple variants. Default: all 15.",
    )
    parser.add_argument("--input", type=Path)
    parser.add_argument("--checkm2-db", type=Path)
    parser.add_argument("--gunc-db", type=Path)
    parser.add_argument("--results-root", type=Path, default=PROJECT_DIR / "tests" / "results")
    parser.add_argument("--work-root", type=Path, default=PROJECT_DIR / "tests" / "work" / "variants")
    parser.add_argument("--nextflow", default="nextflow")
    parser.add_argument("--slurm-account")
    parser.add_argument("--slurm-queue")
    parser.add_argument("--slurm-qos")
    parser.add_argument("--slurm-cluster-options")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")
    parser.add_argument("--skip-summary", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--nextflow-arg",
        action="append",
        default=[],
        help="Additional Nextflow argument; repeat for multiple arguments.",
    )
    return parser.parse_args()


def default_resource_paths() -> tuple[Path, Path, Path]:
    generated = PROJECT_DIR / "tests" / "generated_data"
    return (
        generated / "samplesheet.csv",
        generated / "databases" / "checkm2" / "uniref100.KO.1.dmnd",
        generated / "databases" / "gunc" / "gunc_db.dmnd",
    )


def build_command(
    *,
    nextflow: str,
    variant: str,
    environment: str,
    runtime: str | None,
    mode: str,
    input_path: Path,
    checkm2_db: Path,
    gunc_db: Path,
    output_dir: Path,
    work_dir: Path,
    run_name: str,
    resume: bool,
    slurm_account: str | None,
    slurm_queue: str | None,
    slurm_qos: str | None,
    slurm_cluster_options: str | None,
    extra_args: list[str],
) -> list[str]:
    assembler, binner = split_variant(variant)
    profiles = [environment]
    if runtime:
        profiles.append(runtime)
    if mode == "stub":
        profiles.append("stub")

    command = [
        nextflow,
        "run",
        str(PROJECT_DIR / "tests" / "workflows" / "benchmark.nf"),
        "-c",
        str(PROJECT_DIR / "tests" / "config" / "benchmark.config"),
        "-profile",
        ",".join(profiles),
        "--pipeline_root",
        str(PROJECT_DIR),
        "--benchmark_assembler",
        assembler,
        "--benchmark_binner",
        binner,
        "--input",
        str(input_path),
        "--checkm2_db",
        str(checkm2_db),
        "--gunc_db",
        str(gunc_db),
        "--outdir",
        str(output_dir),
        "-work-dir",
        str(work_dir),
        "-ansi-log",
        "false",
    ]
    if mode == "stub":
        command.append("-stub-run")
    if resume:
        command.extend(["-resume", run_name])
    else:
        command.extend(["-name", run_name])
    for parameter, value in (
        ("--slurm_account", slurm_account),
        ("--slurm_queue", slurm_queue),
        ("--slurm_qos", slurm_qos),
        ("--slurm_cluster_options", slurm_cluster_options),
    ):
        if value:
            command.extend([parameter, value])
    command.extend(extra_args)
    return command


def run_one(
    variant: str,
    command: list[str],
    output_dir: Path,
    run_name: str,
    dry_run: bool,
) -> tuple[str, int]:
    rendered = shlex.join(command)
    if dry_run:
        print(rendered)
        return variant, 0

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "launch_command.txt").write_text(rendered + "\n", encoding="utf-8")
    (output_dir / "nextflow_run_name.txt").write_text(run_name + "\n", encoding="utf-8")
    status_path = output_dir / "benchmark_status.json"
    status = {
        "variant": variant,
        "run_name": run_name,
        "state": "running",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "finished_at": None,
        "returncode": None,
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    log_path = output_dir / "nextflow.log"
    print(f"Starting {variant}")
    try:
        with log_path.open("w", encoding="utf-8", newline="\n") as log_handle:
            completed = subprocess.run(
                command,
                cwd=PROJECT_DIR,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=False,
                text=True,
            )
    except KeyboardInterrupt:
        status.update(
            state="interrupted",
            finished_at=datetime.now(timezone.utc).isoformat(),
            returncode=130,
        )
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise
    status.update(
        state="complete" if completed.returncode == 0 else "failed",
        finished_at=datetime.now(timezone.utc).isoformat(),
        returncode=completed.returncode,
    )
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    state = "completed" if completed.returncode == 0 else "failed"
    print(f"{variant}: {state} (log: {log_path})")
    return variant, completed.returncode


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")
    if args.mode == "run" and not args.runtime:
        raise SystemExit("A software runtime is required for --run")
    if args.stop_on_error and args.jobs > 1:
        raise SystemExit("--stop-on-error requires --jobs 1")
    if not args.dry_run and shutil.which(args.nextflow) is None:
        raise SystemExit(f"Nextflow executable was not found: {args.nextflow}")

    default_input, default_checkm2, default_gunc = default_resource_paths()
    input_path = (args.input or default_input).resolve()
    checkm2_db = (args.checkm2_db or default_checkm2).resolve()
    gunc_db = (args.gunc_db or default_gunc).resolve()
    if not args.dry_run:
        missing = [path for path in (input_path, checkm2_db, gunc_db) if not path.exists()]
        if missing:
            formatted = "\n".join(f"  - {path}" for path in missing)
            raise SystemExit(
                "Required benchmark inputs are missing:\n"
                f"{formatted}\n"
                "Generate local fixtures with tests/scripts/generate_synthetic_data.py "
                "or supply explicit production paths."
            )

    variants = args.variant or variant_names()
    results_root = args.results_root.resolve()
    work_root = args.work_root.resolve()
    jobs: list[tuple[str, list[str], Path, str]] = []
    run_token = datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
    for variant in variants:
        output_dir = results_root / variant
        run_name_path = output_dir / "nextflow_run_name.txt"
        resume_existing = args.resume and run_name_path.is_file()
        if resume_existing:
            run_name = run_name_path.read_text(encoding="utf-8").strip()
        else:
            run_name = f"benchmark_{variant}_{run_token}"
        command = build_command(
            nextflow=args.nextflow,
            variant=variant,
            environment=args.environment,
            runtime=args.runtime,
            mode=args.mode,
            input_path=input_path,
            checkm2_db=checkm2_db,
            gunc_db=gunc_db,
            output_dir=output_dir,
            work_dir=work_root / variant,
            run_name=run_name,
            resume=resume_existing,
            slurm_account=args.slurm_account,
            slurm_queue=args.slurm_queue,
            slurm_qos=args.slurm_qos,
            slurm_cluster_options=args.slurm_cluster_options,
            extra_args=args.nextflow_arg,
        )
        jobs.append((variant, command, output_dir, run_name))

    failures: list[str] = []
    if args.jobs == 1:
        for variant, command, output_dir, run_name in jobs:
            _, returncode = run_one(variant, command, output_dir, run_name, args.dry_run)
            if returncode:
                failures.append(variant)
                if args.stop_on_error:
                    break
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(run_one, variant, command, output_dir, run_name, args.dry_run): variant
                for variant, command, output_dir, run_name in jobs
            }
            for future in concurrent.futures.as_completed(futures):
                variant, returncode = future.result()
                if returncode:
                    failures.append(variant)

    if not args.skip_summary and not args.dry_run:
        summary_command = [
            sys.executable,
            str(PROJECT_DIR / "tests" / "scripts" / "summarize_variants.py"),
            "--results-root",
            str(results_root),
        ]
        summary = subprocess.run(summary_command, cwd=PROJECT_DIR, check=False)
        if summary.returncode:
            failures.append("variant_summary")

    if failures:
        print("Failed: " + ", ".join(sorted(failures)), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
