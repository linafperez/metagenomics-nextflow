from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "bin" / "summarize_resources.py"
FIXTURES = Path(__file__).resolve().parent / "fixtures" / "resource_accounting"
SPEC = importlib.util.spec_from_file_location("summarize_resources", SCRIPT)
assert SPEC and SPEC.loader
summarize_resources = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = summarize_resources
SPEC.loader.exec_module(summarize_resources)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class UnitParserTests(unittest.TestCase):
    def test_duration_memory_percent_and_timestamp_normalization(self) -> None:
        self.assertEqual(summarize_resources.parse_duration("1500"), 1.5)
        self.assertEqual(
            summarize_resources.parse_duration("1h 2m 3.5s", raw_numeric=False),
            3723.5,
        )
        self.assertEqual(summarize_resources.parse_duration("01:02:03"), 3723)
        self.assertEqual(summarize_resources.parse_bytes("1.5 GB"), 1610612736)
        self.assertEqual(summarize_resources.parse_bytes("4096"), 4096)
        self.assertEqual(summarize_resources.parse_percent("150%"), 150)
        self.assertEqual(
            summarize_resources.timestamp_text("1704067200000"),
            "2024-01-01T00:00:00.000Z",
        )

    def test_local_manifest_reports_zero_sra_runs_without_false_warning(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            manifest = Path(temporary) / "samplesheet.csv"
            manifest.write_text(
                "sample,fastq_1,fastq_2\nS1,R1.fastq.gz,R2.fastq.gz\n",
                encoding="utf-8",
            )
            context = summarize_resources.BuildContext()
            counts = summarize_resources.load_manifest_counts(
                manifest, context, "local"
            )
            self.assertEqual(counts, {"biological_samples": 1, "sra_runs": 0})
            self.assertFalse(any("run-accession" in item for item in context.limitations))

    def test_gpu_metrics_are_bound_to_session_and_attempt(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            metrics_dir = Path(temporary)
            header = (
                "timestamp\tprocess\tsample_id\tsession_id\tattempt\tgpu_index\t"
                "gpu_uuid\tgpu_name\tutilization_gpu_percent\tmemory_used_mib\t"
                "memory_total_mib\n"
            )
            (metrics_dir / "first.gpu_metrics.tsv").write_text(
                header
                + "2024-01-01T00:00:00Z\tVAMB\tS1\tsession-a\t1\t0\tGPU-a\tA100\t10\t100\t1000\n",
                encoding="utf-8",
            )
            (metrics_dir / "retry.gpu_metrics.tsv").write_text(
                header
                + "2024-01-01T01:00:00Z\tVAMB\tS1\tsession-b\t2\t0\tGPU-b\tA100\t80\t800\t1000\n",
                encoding="utf-8",
            )
            context = summarize_resources.BuildContext()
            metrics = summarize_resources.load_gpu_metrics(metrics_dir, context)

            def task(session: str, attempt: int) -> object:
                return summarize_resources.TaskRecord(
                    {
                        "session_id": session,
                        "module": "VAMB",
                        "tag": "S1",
                        "attempt": attempt,
                        "gpu_models": "",
                        "gpu_utilization_mean_percent": None,
                        "gpu_utilization_max_percent": None,
                        "peak_gpu_memory_bytes": None,
                        "gpu_metric_samples": 0,
                    },
                    {},
                )

            first = task("session-a", 1)
            retry = task("session-b", 2)
            unrelated = task("session-c", 1)
            summarize_resources.attach_gpu_metrics(
                [first, retry, unrelated], metrics
            )
            self.assertEqual(first.values["gpu_utilization_mean_percent"], 10)
            self.assertEqual(retry.values["gpu_utilization_mean_percent"], 80)
            self.assertEqual(unrelated.values["gpu_metric_samples"], 0)


class ResourceSummaryTests(unittest.TestCase):
    def run_summary(self, output: Path) -> int:
        return summarize_resources.main(
            [
                "--registry",
                str(FIXTURES / "trace_registry.tsv"),
                "--output-dir",
                str(output),
                "--task-peaks",
                str(FIXTURES / "task_workdir_peaks.tsv"),
                "--task-work-timeseries",
                str(FIXTURES / "task_workdir_timeseries.tsv"),
                "--storage-timeseries",
                str(FIXTURES / "storage_usage_timeseries.tsv"),
                "--sra-manifest",
                str(FIXTURES / "sra_manifest.tsv"),
            ]
        )

    def test_outputs_formulas_cached_semantics_scopes_and_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "pipeline_info"
            self.assertEqual(self.run_summary(output), 0)
            required = [
                output / "execution_trace.tsv",
                output / "resources" / "resource_usage_by_task.tsv",
                output / "resources" / "resource_usage_by_process.tsv",
                output / "resources" / "resource_usage_by_subworkflow.tsv",
                output / "resources" / "resource_usage_summary.tsv",
                output / "resources" / "resource_usage_summary.json",
                output / "resources" / "resource_usage_summary.html",
            ]
            for path in required:
                self.assertTrue(path.exists(), path)

            tasks = read_tsv(required[1])
            self.assertEqual(len(tasks), 4)
            cached = next(row for row in tasks if row["status"] == "CACHED")
            self.assertEqual(cached["realtime_seconds"], "3600")
            self.assertEqual(cached["allocated_cpu_hours"], "0")
            self.assertEqual(cached["observed_cpu_hours"], "0")
            self.assertEqual(cached["accounted"], "false")

            fastp = next(
                row
                for row in tasks
                if row["status"] == "COMPLETED" and row["module"] == "FASTP"
            )
            self.assertEqual(fastp["allocated_cpu_hours"], "2")
            self.assertEqual(fastp["observed_cpu_hours"], "1.5")
            self.assertEqual(fastp["peak_rss_bytes"], "1073741824")
            gtdb = next(row for row in tasks if row["module"] == "GTDBTK_CLASSIFY")
            self.assertEqual(gtdb["realtime_seconds"], "1800")
            self.assertEqual(gtdb["requested_memory_bytes"], "17179869184")
            self.assertEqual(gtdb["task_peak_work_bytes"], "20000")

            summary = json.loads(required[5].read_text(encoding="utf-8"))
            self.assertEqual(summary["counts"]["nextflow_invocations"], 2)
            self.assertEqual(summary["counts"]["trace_rows"], 4)
            self.assertEqual(summary["counts"]["executed_task_runs"], 3)
            self.assertEqual(summary["counts"]["cached_task_rows_excluded"], 1)
            self.assertEqual(summary["counts"]["biological_samples"], 2)
            self.assertEqual(summary["counts"]["sra_runs"], 3)
            self.assertAlmostEqual(summary["cpu"]["allocated_cpu_hours"], 8.0)
            self.assertAlmostEqual(summary["cpu"]["observed_cpu_hours"], 4.5)
            self.assertEqual(summary["memory"]["max_peak_rss_bytes"], 3221225472)
            self.assertEqual(summary["memory"]["median_peak_rss_bytes"], 2147483648)
            self.assertEqual(summary["memory"]["mean_peak_rss_bytes"], 2147483648)
            self.assertEqual(summary["time"]["cumulative_task_realtime_seconds"], 7200)
            self.assertEqual(summary["time"]["workflow_wall_time_seconds"], 14400)
            self.assertEqual(summary["storage"]["peak_total_measured_bytes"], 1000)
            self.assertEqual(summary["storage"]["peak_total_dynamic_bytes"], 900)
            self.assertEqual(summary["storage"]["peaks"]["work_bytes"], 800)
            self.assertEqual(summary["storage"]["max_task_peak_work_bytes"], 30000)
            self.assertEqual(
                summary["largest_consumers"][
                    "subworkflow_by_sampled_concurrent_work"
                ],
                "mag_construction",
            )
            self.assertEqual(summary["gpu"]["tasks_requesting_accelerator"], 0)
            self.assertEqual(summary["gpu"]["accelerator_hours"], 0)

            subworkflows = {
                row["subworkflow"] for row in read_tsv(required[3])
            }
            self.assertEqual(
                subworkflows,
                {
                    "quality_control_and_filtering",
                    "mag_construction",
                    "taxonomic_classification_and_phylogenomics",
                },
            )
            mag_scope = next(
                row
                for row in read_tsv(required[3])
                if row["subworkflow"] == "mag_construction"
            )
            self.assertEqual(mag_scope["sampled_peak_concurrent_work_bytes"], "30000")
            merged_header = required[0].read_text(encoding="utf-8").splitlines()[0]
            self.assertIn("raw_%cpu", merged_header.split("\t"))
            self.assertTrue(
                any("never summed" in item for item in summary["limitations"])
            )

    def test_generation_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "pipeline_info"
            self.assertEqual(self.run_summary(output), 0)
            first = {
                path.relative_to(output): path.read_bytes()
                for path in output.rglob("*")
                if path.is_file()
            }
            self.assertEqual(self.run_summary(output), 0)
            second = {
                path.relative_to(output): path.read_bytes()
                for path in output.rglob("*")
                if path.is_file()
            }
            self.assertEqual(first, second)

    def test_missing_optional_metrics_are_reported_not_fabricated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace = root / "minimal.tsv"
            trace.write_text(
                "task_id\tprocess\tstatus\trealtime\tcpus\n"
                "1\tMETAGENOMICS:CHECK_SAMPLESHEET\tCOMPLETED\t1000\t1\n",
                encoding="utf-8",
            )
            registry = root / "registry.tsv"
            registry.write_text(
                "invocation_id\ttrace_path\ttrace_raw\n"
                "only\tminimal.tsv\ttrue\n",
                encoding="utf-8",
            )
            output = root / "out"
            self.assertEqual(
                summarize_resources.main(
                    ["--registry", str(registry), "--output-dir", str(output)]
                ),
                0,
            )
            summary = json.loads(
                (output / "resources" / "resource_usage_summary.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertIsNone(summary["memory"]["max_peak_rss_bytes"])
            joined = " ".join(summary["limitations"])
            self.assertIn("RSS statistics exclude", joined)
            self.assertIn("input manifest was not provided", joined.lower())
            task = read_tsv(output / "resources" / "resource_usage_by_task.tsv")[0]
            self.assertEqual(task["subworkflow"], "input_validation")


if __name__ == "__main__":
    unittest.main()
