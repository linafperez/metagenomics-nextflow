#!/usr/bin/env python3
"""Unit tests for benchmark orchestration and metric extraction."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_variants import build_command, variant_names  # noqa: E402
from summarize_variants import collect_variant_metrics, ranking_key  # noqa: E402


class BenchmarkMatrixTests(unittest.TestCase):
    def test_matrix_contains_exactly_fifteen_unique_variants(self) -> None:
        expected = {
            f"{assembler}_{binner}"
            for assembler in ("megahit", "spades", "both")
            for binner in ("comebin", "metabat2", "semibin2", "vamb", "all")
        }
        self.assertEqual(set(variant_names()), expected)
        self.assertEqual(len(variant_names()), 15)

    def test_stub_command_uses_test_entrypoint_and_variant_profiles(self) -> None:
        command = build_command(
            nextflow="nextflow",
            variant="both_all",
            environment="local",
            runtime="apptainer",
            mode="stub",
            input_path=Path("input.csv"),
            checkm2_db=Path("checkm2.dmnd"),
            gunc_db=Path("gunc.dmnd"),
            output_dir=Path("results"),
            work_dir=Path("work"),
            run_name="benchmark_both_all_saved",
            resume=True,
            slurm_account=None,
            slurm_queue=None,
            slurm_qos=None,
            slurm_cluster_options=None,
            extra_args=[],
        )
        self.assertIn("tests/workflows/benchmark.nf", "/".join(command).replace("\\", "/"))
        self.assertIn("local,apptainer,stub", command)
        self.assertIn("-stub-run", command)
        self.assertIn("-resume", command)
        self.assertNotIn("-name", command)
        self.assertEqual(command[command.index("-resume") + 1], "benchmark_both_all_saved")
        self.assertEqual(command[command.index("--benchmark_assembler") + 1], "both")
        self.assertEqual(command[command.index("--benchmark_binner") + 1], "all")


class BenchmarkMetricTests(unittest.TestCase):
    def write_fixture(self, variant_root: Path) -> None:
        metaquast = (
            variant_root
            / "native"
            / "BENCHMARK"
            / "MEGAHIT_ASSEMBLY"
            / "METAQUAST"
            / "megahit_coassembly.metaquast.report.tsv"
        )
        metaquast.parent.mkdir(parents=True)
        metaquast.write_text(
            "Assembly\tmegahit.fa\n"
            "# contigs\t4\n"
            "Total length\t10000\n"
            "Largest contig\t4000\n"
            "N50\t3000\n"
            "L50\t2\n",
            encoding="utf-8",
        )

        bin_dir = variant_root / "native" / "COMEBIN" / "megahit.comebin.bins"
        bin_dir.mkdir(parents=True)
        for index in (1, 2):
            (bin_dir / f"bin_{index}.fa").write_text(">contig\nACGT\n", encoding="utf-8")

        raw_quality = variant_root / "native" / "CHECKM2_RAW" / "raw.checkm2.quality_report.tsv"
        raw_quality.parent.mkdir(parents=True)
        raw_quality.write_text(
            "Name\tCompleteness\tContamination\n"
            "raw_1\t95\t2\n"
            "raw_2\t88\t3\n"
            "raw_3\t91\t6\n",
            encoding="utf-8",
        )
        clean_quality = variant_root / "native" / "CHECKM2_CLEAN" / "clean.checkm2.quality_report.tsv"
        clean_quality.parent.mkdir(parents=True)
        clean_quality.write_text(
            "Name\tCompleteness\tContamination\n"
            "mag_1\t95\t2\n"
            "mag_2\t88\t3\n",
            encoding="utf-8",
        )

        gunc = variant_root / "native" / "GUNC_CLEAN" / "clean.gunc.summary.tsv"
        gunc.parent.mkdir(parents=True)
        gunc.write_text(
            "genome\tclade_separation_score\tpass.GUNC\n"
            "mag_1\t0.05\tTrue\n"
            "mag_2\t0.15\tFalse\n",
            encoding="utf-8",
        )

        ani = variant_root / "native" / "DREP_ANI_99" / "catalog_ani99.clusters.csv"
        ani.parent.mkdir(parents=True)
        ani.write_text(
            "genome,primary_cluster,secondary_cluster\n"
            "a.fa,primary_1,secondary_1\n"
            "b.fa,primary_1,secondary_1\n"
            "c.fa,primary_1,secondary_2\n",
            encoding="utf-8",
        )
        species = variant_root / "native" / "DREP_SPECIES_95" / "catalog_species95.clusters.csv"
        species.parent.mkdir(parents=True)
        species.write_text(
            "genome,primary_cluster,secondary_cluster\n"
            "a.fa,primary_1,species_1\n"
            "c.fa,primary_1,species_1\n",
            encoding="utf-8",
        )

        pipeline_info = variant_root / "pipeline_info"
        pipeline_info.mkdir(parents=True)
        (pipeline_info / "execution_trace.txt").write_text(
            "task_id\tstatus\tstart\tcomplete\trealtime\tcpus\tpeak_rss\tattempt\n"
            "1\tCOMPLETED\t2026-08-22T10:00:00+00:00\t2026-08-22T10:00:02+00:00\t2s\t1\t10 MB\t1\n"
            "2\tCOMPLETED\t2026-08-22T10:00:02+00:00\t2026-08-22T10:00:05+00:00\t3s\t2\t12 MB\t2\n",
            encoding="utf-8",
        )
        (pipeline_info / "megahit_comebin.multiqc.html").write_text(
            "<!doctype html><title>MultiQC</title>\n",
            encoding="utf-8",
        )

    def test_native_metric_parsing_and_missing_value_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            variant_root = Path(temporary) / "megahit_comebin"
            self.write_fixture(variant_root)
            metrics = collect_variant_metrics(variant_root, "megahit_comebin")

            self.assertEqual(metrics["status"], "complete")
            self.assertEqual(metrics["assembly_contigs"], 4)
            self.assertEqual(metrics["comebin_bin_count"], 2)
            self.assertEqual(metrics["raw_bin_count"], 2)
            self.assertEqual(metrics["raw_mag_count"], 3)
            self.assertEqual(metrics["high_quality_mag_count"], 1)
            self.assertEqual(metrics["median_completeness"], 91.5)
            self.assertEqual(metrics["median_contamination"], 2.5)
            self.assertEqual(metrics["gunc_pass_count"], 1)
            self.assertEqual(metrics["gunc_fail_count"], 1)
            self.assertAlmostEqual(metrics["gunc_mean_clade_separation_score"], 0.1)
            self.assertEqual(metrics["derep_input_mag_count"], 3)
            self.assertEqual(metrics["derep_99_representative_count"], 2)
            self.assertEqual(metrics["species_95_group_count"], 1)
            self.assertEqual(metrics["wall_clock_seconds"], 5)
            self.assertEqual(metrics["total_task_seconds"], 5)
            self.assertEqual(metrics["cpus_requested_sum"], 3)
            self.assertEqual(metrics["peak_rss_bytes"], 12_000_000)
            self.assertEqual(metrics["retried_tasks"], 1)
            self.assertIsNone(metrics["distinct_gtdb_species"])
            self.assertEqual(
                metrics["multiqc_report"],
                "pipeline_info/megahit_comebin.multiqc.html",
            )

    def test_ranking_prefers_biological_quality_before_runtime(self) -> None:
        stronger = {
            "variant": "megahit_comebin",
            "high_quality_mag_count": 2,
            "median_completeness": 95,
            "median_contamination": 2,
            "gunc_fail_count": 0,
            "final_mag_count": 2,
            "wall_clock_seconds": 100,
        }
        faster = {
            "variant": "spades_vamb",
            "high_quality_mag_count": 1,
            "median_completeness": 99,
            "median_contamination": 1,
            "gunc_fail_count": 0,
            "final_mag_count": 3,
            "wall_clock_seconds": 1,
        }
        self.assertEqual(sorted([faster, stronger], key=ranking_key)[0], stronger)

    def test_runner_state_preserves_incomplete_execution_status(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            variant_root = Path(temporary) / "megahit_comebin"
            self.write_fixture(variant_root)
            (variant_root / "benchmark_status.json").write_text(
                json.dumps({"state": "running"}) + "\n",
                encoding="utf-8",
            )
            metrics = collect_variant_metrics(variant_root, "megahit_comebin")
            self.assertEqual(metrics["status"], "running")

    def test_summary_command_writes_all_structured_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            results_root = Path(temporary) / "results"
            self.write_fixture(results_root / "megahit_comebin")
            output_dir = Path(temporary) / "comparison"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_DIR / "summarize_variants.py"),
                    "--results-root",
                    str(results_root),
                    "--output-dir",
                    str(output_dir),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            for filename in (
                "variant_comparison.tsv",
                "variant_ranking.tsv",
                "variant_summary.md",
                "variant_comparison.json",
            ):
                self.assertTrue((output_dir / filename).is_file())
            payload = json.loads((output_dir / "variant_comparison.json").read_text(encoding="utf-8"))
            self.assertEqual(len(payload["variants"]), 15)


if __name__ == "__main__":
    unittest.main()
