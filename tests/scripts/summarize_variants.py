#!/usr/bin/env python3
"""Extract native benchmark metrics and produce a transparent ranking."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable


ASSEMBLERS = ("megahit", "spades", "both")
BINNERS = ("comebin", "metabat2", "semibin2", "vamb", "all")
HQ_COMPLETENESS = 90.0
HQ_CONTAMINATION = 5.0

COMPARISON_COLUMNS = [
    "variant",
    "assembler",
    "binner",
    "status",
    "megahit_contigs",
    "megahit_total_length",
    "megahit_largest_contig",
    "megahit_n50",
    "megahit_l50",
    "spades_contigs",
    "spades_total_length",
    "spades_largest_contig",
    "spades_n50",
    "spades_l50",
    "assembly_contigs",
    "total_assembly_length",
    "largest_contig",
    "n50",
    "l50",
    "raw_bin_count",
    "comebin_bin_count",
    "metabat2_bin_count",
    "semibin2_bin_count",
    "vamb_bin_count",
    "dastool_refined_bin_count",
    "raw_mag_count",
    "total_mag_count",
    "high_quality_mag_count",
    "mean_completeness",
    "median_completeness",
    "minimum_completeness",
    "mean_contamination",
    "median_contamination",
    "maximum_contamination",
    "gunc_pass_count",
    "gunc_fail_count",
    "gunc_mean_clade_separation_score",
    "derep_input_mag_count",
    "derep_99_representative_count",
    "species_95_group_count",
    "species_95_representative_count",
    "final_mag_count",
    "final_hq_mag_count",
    "distinct_gtdb_taxa",
    "distinct_gtdb_species",
    "wall_clock_seconds",
    "total_task_seconds",
    "cpus_requested_sum",
    "peak_rss_bytes",
    "failed_tasks",
    "retried_tasks",
    "multiqc_report",
]

RANKING_DESCRIPTION = [
    "higher high-quality MAG count (>90% completeness and <5% contamination)",
    "higher median completeness",
    "lower median contamination",
    "fewer GUNC failures",
    "higher final non-redundant MAG count",
    "lower wall-clock runtime as a tie-breaker",
]


def variant_names() -> list[str]:
    return [f"{assembler}_{binner}" for assembler in ASSEMBLERS for binner in BINNERS]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize the 15 benchmark variants.")
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def as_number(value: str) -> float | int | None:
    cleaned = value.strip().replace(",", "")
    if not cleaned or cleaned.lower() in {"na", "n/a", "nan", "none", "-"}:
        return None
    try:
        number = float(cleaned)
    except ValueError:
        return None
    if math.isfinite(number) and number.is_integer():
        return int(number)
    return number if math.isfinite(number) else None


def read_delimited(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            return [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error):
        return []


def first_column(row: dict[str, str], *names: str) -> str | None:
    lookup = {normalized(key): value for key, value in row.items() if key is not None}
    for name in names:
        value = lookup.get(normalized(name))
        if value not in (None, ""):
            return value
    return None


def parse_metaquast(path: Path) -> dict[str, int | float | None]:
    metrics: dict[str, str] = {}
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) >= 2:
                    metrics[row[0].strip()] = row[1].strip()
    except (OSError, UnicodeDecodeError, csv.Error):
        return {}

    def metric(*names: str) -> int | float | None:
        for name in names:
            if name in metrics:
                return as_number(metrics[name])
        return None

    return {
        "contigs": metric("# contigs", "# contigs (>= 0 bp)"),
        "total_length": metric("Total length", "Total length (>= 0 bp)"),
        "largest_contig": metric("Largest contig"),
        "n50": metric("N50"),
        "l50": metric("L50"),
    }


def scan_files(root: Path) -> list[Path]:
    files: list[Path] = []
    pending = [root]
    seen_directories: set[tuple[int, int]] = set()
    while pending:
        directory = pending.pop()
        try:
            directory_stat = directory.stat()
            identity = (directory_stat.st_dev, directory_stat.st_ino)
            if identity in seen_directories:
                continue
            seen_directories.add(identity)
            entries = list(os.scandir(directory))
        except OSError:
            continue
        for entry in entries:
            path = Path(entry.path)
            try:
                if entry.is_dir(follow_symlinks=True):
                    pending.append(path)
                elif entry.is_file(follow_symlinks=True):
                    files.append(path)
            except OSError:
                continue
    return files


def find_metaquast_reports(files: Iterable[Path]) -> dict[str, Path]:
    file_list = list(files)
    candidates = sorted(
        path for path in file_list if path.name.endswith(".metaquast.report.tsv")
    )
    if not candidates:
        candidates = sorted(path for path in file_list if path.name == "report.tsv")
    reports: dict[str, Path] = {}
    for path in candidates:
        lower = str(path).lower()
        for assembler in ("megahit", "spades"):
            if assembler in lower and assembler not in reports:
                reports[assembler] = path
    return reports


def fasta_paths_in_named_directory(files: Iterable[Path], marker: str) -> set[Path]:
    paths: set[Path] = set()
    for path in files:
        if (
            path.suffix.lower() in {".fa", ".fna", ".fasta"}
            and marker.lower() in path.parent.name.lower()
        ):
            paths.add(path)
    return paths


def choose_stage_file(paths: Iterable[Path], stages: list[str]) -> Path | None:
    candidates = list(paths)
    for stage in stages:
        normalized_stage = normalized(stage)
        matches = [path for path in candidates if normalized_stage in normalized(str(path))]
        if matches:
            return sorted(matches)[0]
    return sorted(candidates)[0] if candidates else None


def parse_checkm2(path: Path | None) -> list[dict[str, float | str]]:
    if path is None:
        return []
    records: list[dict[str, float | str]] = []
    for row in read_delimited(path):
        name = first_column(row, "Name", "genome", "mag_id")
        completeness = as_number(first_column(row, "Completeness") or "")
        contamination = as_number(first_column(row, "Contamination") or "")
        if name is None or completeness is None or contamination is None:
            continue
        records.append(
            {
                "name": name,
                "completeness": float(completeness),
                "contamination": float(contamination),
            }
        )
    return records


def quality_statistics(records: list[dict[str, float | str]]) -> dict[str, int | float | None]:
    if not records:
        return {
            "count": None,
            "hq_count": None,
            "mean_completeness": None,
            "median_completeness": None,
            "minimum_completeness": None,
            "mean_contamination": None,
            "median_contamination": None,
            "maximum_contamination": None,
        }
    completeness = [float(row["completeness"]) for row in records]
    contamination = [float(row["contamination"]) for row in records]
    return {
        "count": len(records),
        "hq_count": sum(
            complete > HQ_COMPLETENESS and contaminate < HQ_CONTAMINATION
            for complete, contaminate in zip(completeness, contamination)
        ),
        "mean_completeness": statistics.fmean(completeness),
        "median_completeness": statistics.median(completeness),
        "minimum_completeness": min(completeness),
        "mean_contamination": statistics.fmean(contamination),
        "median_contamination": statistics.median(contamination),
        "maximum_contamination": max(contamination),
    }


def parse_gunc(path: Path | None) -> dict[str, int | float | None]:
    if path is None:
        return {"pass": None, "fail": None, "mean_css": None}
    rows = read_delimited(path)
    if not rows:
        return {"pass": None, "fail": None, "mean_css": None}
    statuses: list[bool] = []
    scores: list[float] = []
    for row in rows:
        status = first_column(row, "pass.GUNC", "pass_gunc", "pass")
        if status is not None:
            lowered = status.strip().lower()
            if lowered in {"true", "pass", "passed", "yes", "1"}:
                statuses.append(True)
            elif lowered in {"false", "fail", "failed", "no", "0"}:
                statuses.append(False)
        score = as_number(first_column(row, "clade_separation_score", "CSS") or "")
        if score is not None:
            scores.append(float(score))
    return {
        "pass": sum(statuses) if statuses else None,
        "fail": sum(not status for status in statuses) if statuses else None,
        "mean_css": statistics.fmean(scores) if scores else None,
    }


def parse_drep_clusters(path: Path | None) -> dict[str, int | None]:
    if path is None:
        return {"input": None, "groups": None}
    rows = read_delimited(path, delimiter=",")
    genomes = {
        first_column(row, "genome")
        for row in rows
        if first_column(row, "genome") is not None
    }
    groups = {
        first_column(row, "secondary_cluster")
        for row in rows
        if first_column(row, "secondary_cluster") is not None
    }
    return {
        "input": len(genomes) if genomes else None,
        "groups": len(groups) if groups else None,
    }


def parse_gtdb(files: Iterable[Path]) -> tuple[int | None, int | None]:
    classifications: set[str] = set()
    species: set[str] = set()
    for path in files:
        if not path.name.endswith("summary.tsv"):
            continue
        if "gtdb" not in str(path).lower():
            continue
        for row in read_delimited(path):
            classification = first_column(row, "classification")
            if not classification:
                continue
            classifications.add(classification)
            for taxon in classification.split(";"):
                if taxon.startswith("s__") and taxon != "s__":
                    species.add(taxon)
    return (
        len(classifications) if classifications else None,
        len(species) if species else None,
    )


DURATION_FACTORS = {
    "ms": 0.001,
    "s": 1.0,
    "m": 60.0,
    "h": 3600.0,
    "d": 86400.0,
}
MEMORY_FACTORS = {
    "b": 1,
    "kb": 1000,
    "mb": 1000**2,
    "gb": 1000**3,
    "tb": 1000**4,
    "kib": 1024,
    "mib": 1024**2,
    "gib": 1024**3,
    "tib": 1024**4,
}


def parse_scaled(value: str, factors: dict[str, float | int]) -> float | None:
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)\s*([A-Za-z]+)\s*", value or "")
    if not match:
        return None
    unit = match.group(2).lower()
    if unit not in factors:
        return None
    return float(match.group(1)) * float(factors[unit])


def parse_duration(value: str) -> float | None:
    cleaned = (value or "").strip()
    if not cleaned:
        return None
    matches = re.findall(r"([0-9]+(?:\.[0-9]+)?)\s*(ms|[smhd])", cleaned.lower())
    if not matches:
        return None
    return sum(float(number) * DURATION_FACTORS[unit] for number, unit in matches)


def parse_timestamp(value: str) -> datetime | None:
    cleaned = (value or "").strip()
    if not cleaned:
        return None
    try:
        return datetime.fromisoformat(cleaned.replace("Z", "+00:00"))
    except ValueError:
        return None


def parse_trace(path: Path | None) -> dict[str, int | float | None]:
    empty = {
        "wall_clock": None,
        "task_seconds": None,
        "cpus": None,
        "peak_rss": None,
        "failed": None,
        "retried": None,
    }
    if path is None:
        return empty
    rows = read_delimited(path)
    if not rows:
        return empty
    starts: list[datetime] = []
    completes: list[datetime] = []
    submits: list[datetime] = []
    derived_completes: list[datetime] = []
    runtimes: list[float] = []
    cpus: list[int] = []
    rss_values: list[float] = []
    failed = 0
    retried = 0
    status_seen = False
    attempt_seen = False
    for row in rows:
        submit = parse_timestamp(first_column(row, "submit") or "")
        start = parse_timestamp(first_column(row, "start") or "")
        complete = parse_timestamp(first_column(row, "complete") or "")
        duration = parse_duration(first_column(row, "duration") or "")
        if submit:
            submits.append(submit)
            if duration is not None:
                derived_completes.append(submit + timedelta(seconds=duration))
        if start:
            starts.append(start)
        if complete:
            completes.append(complete)
        runtime = parse_duration(first_column(row, "realtime", "duration") or "")
        if runtime is not None:
            runtimes.append(runtime)
        requested = as_number(first_column(row, "cpus") or "")
        if requested is not None:
            cpus.append(int(requested))
        rss = parse_scaled(first_column(row, "peak_rss", "rss") or "", MEMORY_FACTORS)
        if rss is not None:
            rss_values.append(rss)
        status = (first_column(row, "status") or "").upper()
        if status:
            status_seen = True
            if status not in {"COMPLETED", "CACHED"}:
                failed += 1
        attempt = as_number(first_column(row, "attempt") or "")
        if attempt is not None:
            attempt_seen = True
            if int(attempt) > 1:
                retried += int(attempt) - 1
    wall_clock = None
    if starts and completes:
        wall_clock = (max(completes) - min(starts)).total_seconds()
    elif submits and derived_completes:
        wall_clock = (max(derived_completes) - min(submits)).total_seconds()
    return {
        "wall_clock": wall_clock,
        "task_seconds": sum(runtimes) if runtimes else None,
        "cpus": sum(cpus) if cpus else None,
        "peak_rss": max(rss_values) if rss_values else None,
        "failed": failed if status_seen else None,
        "retried": retried if attempt_seen else None,
    }


def collect_variant_metrics(root: Path, variant: str) -> dict[str, Any]:
    assembler, binner = variant.split("_", maxsplit=1)
    metrics: dict[str, Any] = {column: None for column in COMPARISON_COLUMNS}
    metrics.update(
        {
            "variant": variant,
            "assembler": assembler,
            "binner": binner,
            "status": "missing" if not root.is_dir() else "incomplete",
        }
    )
    if not root.is_dir():
        return metrics

    files = scan_files(root)
    reports = find_metaquast_reports(files)
    assembly_values: dict[str, dict[str, int | float | None]] = {}
    for assembly_name, report in reports.items():
        parsed = parse_metaquast(report)
        assembly_values[assembly_name] = parsed
        for field in ("contigs", "total_length", "largest_contig", "n50", "l50"):
            metrics[f"{assembly_name}_{field}"] = parsed.get(field)

    selected_assemblies = [name for name in ("megahit", "spades") if name in assembly_values]
    if len(selected_assemblies) == 1:
        parsed = assembly_values[selected_assemblies[0]]
        metrics.update(
            {
                "assembly_contigs": parsed.get("contigs"),
                "total_assembly_length": parsed.get("total_length"),
                "largest_contig": parsed.get("largest_contig"),
                "n50": parsed.get("n50"),
                "l50": parsed.get("l50"),
            }
        )
    elif len(selected_assemblies) == 2:
        contigs = [assembly_values[name].get("contigs") for name in selected_assemblies]
        lengths = [assembly_values[name].get("total_length") for name in selected_assemblies]
        largest = [assembly_values[name].get("largest_contig") for name in selected_assemblies]
        if all(value is not None for value in contigs):
            metrics["assembly_contigs"] = sum(contigs)  # type: ignore[arg-type]
        if all(value is not None for value in lengths):
            metrics["total_assembly_length"] = sum(lengths)  # type: ignore[arg-type]
        if any(value is not None for value in largest):
            metrics["largest_contig"] = max(value for value in largest if value is not None)
        # N50 and L50 cannot be combined correctly without the original contig lengths.

    bin_markers = {
        "comebin": ".comebin.bins",
        "metabat2": ".metabat2.bins",
        "semibin2": ".semibin2.bins",
        "vamb": ".vamb.bins",
    }
    bin_counts: dict[str, int] = {}
    for name, marker in bin_markers.items():
        count = len(fasta_paths_in_named_directory(files, marker))
        bin_counts[name] = count
        metrics[f"{name}_bin_count"] = count if count else None
    dastool_count = len(fasta_paths_in_named_directory(files, "_dastool_bins"))
    metrics["dastool_refined_bin_count"] = dastool_count if dastool_count else None
    if binner == "all":
        measured = [bin_counts[name] for name in BINNERS if name != "all"]
        metrics["raw_bin_count"] = sum(measured) if any(measured) else None
    else:
        metrics["raw_bin_count"] = bin_counts[binner] or None

    quality_paths = sorted(
        path for path in files if path.name.endswith(".checkm2.quality_report.tsv")
    )
    if assembler == "both":
        final_quality_path = choose_stage_file(quality_paths, ["CHECKM2_FINAL", "CHECKM2_CLEAN"])
    else:
        final_quality_path = choose_stage_file(quality_paths, ["CHECKM2_CLEAN", "CHECKM2_FINAL"])
    raw_quality_path = choose_stage_file(quality_paths, ["CHECKM2_RAW"])
    final_quality = quality_statistics(parse_checkm2(final_quality_path))
    raw_quality = quality_statistics(parse_checkm2(raw_quality_path))
    metrics["raw_mag_count"] = raw_quality["count"]
    metrics["total_mag_count"] = final_quality["count"]
    metrics["high_quality_mag_count"] = final_quality["hq_count"]
    for field in (
        "mean_completeness",
        "median_completeness",
        "minimum_completeness",
        "mean_contamination",
        "median_contamination",
        "maximum_contamination",
    ):
        metrics[field] = final_quality[field]

    gunc_paths = sorted(path for path in files if path.name.endswith(".gunc.summary.tsv"))
    if assembler == "both":
        gunc_path = choose_stage_file(gunc_paths, ["GUNC_FINAL", "GUNC_CLEAN"])
    else:
        gunc_path = choose_stage_file(gunc_paths, ["GUNC_CLEAN", "GUNC_FINAL"])
    gunc = parse_gunc(gunc_path)
    metrics["gunc_pass_count"] = gunc["pass"]
    metrics["gunc_fail_count"] = gunc["fail"]
    metrics["gunc_mean_clade_separation_score"] = gunc["mean_css"]

    cluster_paths = sorted(path for path in files if path.name.endswith(".clusters.csv"))
    if assembler == "both":
        ani_path = choose_stage_file(cluster_paths, ["DREP_FINAL_99", "final_catalog_ani99"])
        species_path = choose_stage_file(cluster_paths, ["DREP_FINAL_SPECIES_95", "final_catalog_species95"])
    else:
        ani_path = choose_stage_file(cluster_paths, ["DREP_ANI_99", "_ani99"])
        species_path = choose_stage_file(cluster_paths, ["DREP_SPECIES_95", "_species95"])
    ani = parse_drep_clusters(ani_path)
    species = parse_drep_clusters(species_path)
    metrics["derep_input_mag_count"] = ani["input"]
    metrics["derep_99_representative_count"] = ani["groups"]
    metrics["species_95_group_count"] = species["groups"]
    metrics["species_95_representative_count"] = species["groups"]
    metrics["final_mag_count"] = (
        final_quality["count"]
        if final_quality["count"] is not None
        else ani["groups"]
    )
    metrics["final_hq_mag_count"] = final_quality["hq_count"]

    distinct_taxa, distinct_species = parse_gtdb(files)
    metrics["distinct_gtdb_taxa"] = distinct_taxa
    metrics["distinct_gtdb_species"] = distinct_species

    trace_candidates = sorted(path for path in files if path.name == "execution_trace.txt")
    trace = parse_trace(trace_candidates[0] if trace_candidates else None)
    metrics["wall_clock_seconds"] = trace["wall_clock"]
    metrics["total_task_seconds"] = trace["task_seconds"]
    metrics["cpus_requested_sum"] = trace["cpus"]
    metrics["peak_rss_bytes"] = trace["peak_rss"]
    metrics["failed_tasks"] = trace["failed"]
    metrics["retried_tasks"] = trace["retried"]

    multiqc = sorted(
        path
        for path in files
        if path.name == "multiqc_report.html" or path.name.endswith(".multiqc.html")
    )
    metrics["multiqc_report"] = str(multiqc[0].relative_to(root)) if multiqc else None
    has_native_metrics = bool(reports or quality_paths or cluster_paths)
    run_state = None
    status_files = [path for path in files if path.name == "benchmark_status.json"]
    if status_files:
        try:
            run_state = json.loads(status_files[0].read_text(encoding="utf-8")).get("state")
        except (OSError, json.JSONDecodeError):
            run_state = None
    if run_state in {"failed", "interrupted", "running"}:
        metrics["status"] = run_state
    elif trace["failed"]:
        metrics["status"] = "failed"
    elif has_native_metrics:
        metrics["status"] = "complete"
    return metrics


def descending_key(value: Any) -> tuple[bool, float]:
    return (value is None, -float(value) if value is not None else 0.0)


def ascending_key(value: Any) -> tuple[bool, float]:
    return (value is None, float(value) if value is not None else 0.0)


def ranking_key(metrics: dict[str, Any]) -> tuple[Any, ...]:
    return (
        descending_key(metrics["high_quality_mag_count"]),
        descending_key(metrics["median_completeness"]),
        ascending_key(metrics["median_contamination"]),
        ascending_key(metrics["gunc_fail_count"]),
        descending_key(metrics["final_mag_count"]),
        ascending_key(metrics["wall_clock_seconds"]),
        metrics["variant"],
    )


def format_value(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(format_value(row.get(column)) for column in columns)


def write_markdown(path: Path, ranked: list[dict[str, Any]]) -> None:
    lines = [
        "# Benchmark variant summary",
        "",
        "The ranking is lexicographic and does not use an opaque weighted score:",
        "",
    ]
    lines.extend(f"{index}. {description}." for index, description in enumerate(RANKING_DESCRIPTION, start=1))
    lines.extend(
        [
            "",
            "Missing measurements rank after measured values at each comparison step. "
            "For `both`, contig counts and assembly lengths are summed and the largest contig is the maximum; "
            "N50 and L50 remain missing because they cannot be combined correctly from summary statistics.",
            "",
            "| Rank | Variant | Status | HQ MAGs | Median completeness | Median contamination | GUNC failures | Final MAGs | Wall clock (s) |",
            "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in ranked:
        lines.append(
            "| {rank} | {variant} | {status} | {hq} | {completeness} | {contamination} | {gunc} | {mags} | {runtime} |".format(
                rank=row["rank"],
                variant=row["variant"],
                status=row["status"],
                hq=format_value(row["high_quality_mag_count"]),
                completeness=format_value(row["median_completeness"]),
                contamination=format_value(row["median_contamination"]),
                gunc=format_value(row["gunc_fail_count"]),
                mags=format_value(row["final_mag_count"]),
                runtime=format_value(row["wall_clock_seconds"]),
            )
        )
    measured = [row for row in ranked if row["high_quality_mag_count"] is not None]
    lines.extend(["", "## Interpretation", ""])
    if measured:
        top = measured[0]
        lines.append(
            f"`{top['variant']}` ranks first under the documented measured criteria. "
            "This ranking describes the recovered outputs and does not establish broader biological superiority."
        )
    else:
        lines.append("No variant has sufficient measured MAG-quality data for a biological comparison.")
    lines.extend(["", "`NA` denotes a metric that was not present in the available native outputs.", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    results_root = args.results_root.resolve()
    output_dir = (args.output_dir or results_root / "comparison").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    comparison = [
        collect_variant_metrics(results_root / variant, variant)
        for variant in variant_names()
    ]
    ordered = sorted(comparison, key=ranking_key)
    ranked = [dict(row, rank=index) for index, row in enumerate(ordered, start=1)]

    write_tsv(output_dir / "variant_comparison.tsv", comparison, COMPARISON_COLUMNS)
    write_tsv(output_dir / "variant_ranking.tsv", ranked, ["rank", *COMPARISON_COLUMNS])
    write_markdown(output_dir / "variant_summary.md", ranked)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "quality_criteria": {
            "completeness_strictly_greater_than": HQ_COMPLETENESS,
            "contamination_strictly_less_than": HQ_CONTAMINATION,
        },
        "ranking_strategy": RANKING_DESCRIPTION,
        "variants": ranked,
    }
    (output_dir / "variant_comparison.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote benchmark comparison to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
