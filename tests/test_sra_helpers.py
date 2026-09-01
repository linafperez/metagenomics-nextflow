from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
RESOLVER = REPOSITORY / "bin" / "resolve_sra_project.py"
ACQUIRE = REPOSITORY / "bin" / "acquire_sra_sample.py"
CHECKPOINTS = REPOSITORY / "bin" / "manage_sra_checkpoints.py"
FIXTURES = Path(__file__).with_name("fixtures")


def run_command(arguments: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *arguments],
        text=True,
        capture_output=True,
        check=False,
        **kwargs,
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ResolverTests(unittest.TestCase):
    def test_offline_resolution_is_deterministic_and_preserves_identity(self) -> None:
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            fixture = FIXTURES / "runinfo_valid.csv"
            for output in (first, second):
                result = run_command(
                    [
                        str(RESOLVER),
                        "PRJNA123456",
                        "--runinfo-file",
                        str(fixture),
                        "--output-dir",
                        output,
                    ]
                )
                self.assertEqual(result.returncode, 0, result.stderr)

            report_names = (
                "sra_project_runinfo.csv",
                "sra_project_manifest.tsv",
                "sra_sample_manifest.tsv",
                "sra_project_exclusions.tsv",
                "sra_project_summary.json",
            )
            for name in report_names:
                self.assertEqual((Path(first) / name).read_bytes(), (Path(second) / name).read_bytes())
            self.assertEqual(
                (Path(first) / "sra_project_runinfo.csv").read_bytes(), fixture.read_bytes()
            )

            rows = read_tsv(Path(first) / "sra_project_manifest.tsv")
            self.assertEqual(
                [row["run_accession"] for row in rows],
                ["DRR000004", "ERR000003", "SRR000001", "SRR000002"],
            )
            by_run = {row["run_accession"]: row for row in rows}
            self.assertEqual(by_run["SRR000001"]["sample_id"], "SAMN000001")
            self.assertEqual(by_run["SRR000002"]["sample_id"], "SAMN000001")
            self.assertEqual(by_run["ERR000003"]["sample_id"], "ERX000003")
            self.assertEqual(by_run["ERR000003"]["identity_source"], "Experiment")
            self.assertIn("missing_biosample_used_experiment", by_run["ERR000003"]["metadata_warnings"])
            self.assertEqual(by_run["DRR000004"]["sample_id"], "DRR000004")
            self.assertEqual(by_run["DRR000004"]["identity_source"], "Run")
            self.assertIn("missing_experiment_used_run", by_run["DRR000004"]["metadata_warnings"])
            samples = read_tsv(Path(first) / "sra_sample_manifest.tsv")
            self.assertEqual(len(samples), 3)
            biosample = next(row for row in samples if row["sample_id"] == "SAMN000001")
            self.assertEqual(biosample["run_count"], "2")
            self.assertEqual(biosample["run_accessions"], "SRR000001;SRR000002")
            summary = json.loads((Path(first) / "sra_project_summary.json").read_text())
            self.assertTrue(summary["valid"])
            self.assertEqual(summary["eligible_run_count"], 4)
            self.assertEqual(summary["excluded_run_count"], 0)

    def test_invalid_project_writes_reports_without_contacting_ncbi(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = run_command(
                [
                    str(RESOLVER),
                    "SRP123",
                    "--output-dir",
                    temporary,
                    "--write-invalid-and-succeed",
                ]
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((Path(temporary) / "sra_project_summary.json").read_text())
            self.assertFalse(summary["valid"])
            self.assertIn("invalid BioProject accession", summary["validation_errors"][0])

    def test_incompatible_and_duplicate_runs_are_reported_and_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = run_command(
                [
                    str(RESOLVER),
                    "--project",
                    "PRJNA123456",
                    "--runinfo-file",
                    str(FIXTURES / "runinfo_invalid.csv"),
                    "--output-dir",
                    temporary,
                ]
            )
            self.assertEqual(result.returncode, 2)
            exclusions = read_tsv(Path(temporary) / "sra_project_exclusions.tsv")
            self.assertEqual(len(exclusions), 8)
            reasons = ";".join(row["exclusion_reason"] for row in exclusions)
            for reason in (
                "library_layout_not_paired",
                "library_strategy_not_wgs",
                "library_source_not_metagenomic",
                "unsupported_platform",
                "controlled_or_restricted_access",
                "spots_with_mates_does_not_equal_spots",
                "duplicate_run_accession",
            ):
                self.assertIn(reason, reasons)
            summary = json.loads((Path(temporary) / "sra_project_summary.json").read_text())
            self.assertFalse(summary["valid"])
            self.assertEqual(summary["eligible_run_count"], 1)
            self.assertEqual(summary["excluded_run_count"], 8)

    def test_write_invalid_and_succeed_and_validate_existing_modes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            invalid = run_command(
                [
                    str(RESOLVER),
                    "PRJNA123456",
                    "--runinfo-file",
                    str(FIXTURES / "runinfo_invalid.csv"),
                    "--output-dir",
                    str(output),
                    "--write-invalid-and-succeed",
                ]
            )
            self.assertEqual(invalid.returncode, 0, invalid.stderr)
            validation = run_command([str(RESOLVER), "--validate-existing", str(output)])
            self.assertEqual(validation.returncode, 2)
            self.assertIn("marked invalid", validation.stderr)

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            resolved = run_command(
                [
                    str(RESOLVER),
                    "PRJNA123456",
                    "--runinfo-file",
                    str(FIXTURES / "runinfo_valid.csv"),
                    "--output-dir",
                    str(output),
                ]
            )
            self.assertEqual(resolved.returncode, 0, resolved.stderr)
            validation = run_command([str(RESOLVER), "--validate-existing", str(output)])
            self.assertEqual(validation.returncode, 0, validation.stderr)
            manifest = output / "sra_project_manifest.tsv"
            manifest.write_text(manifest.read_text() + "\n", encoding="utf-8")
            tampered = run_command([str(RESOLVER), "--validate-existing", str(output)])
            self.assertEqual(tampered.returncode, 2)
            self.assertIn("SHA-256 mismatch", tampered.stderr)


class AcquisitionTests(unittest.TestCase):
    def resolve_manifest(self, root: Path) -> Path:
        report_dir = root / "reports"
        result = run_command(
            [
                str(RESOLVER),
                "PRJNA123456",
                "--runinfo-file",
                str(FIXTURES / "runinfo_valid.csv"),
                "--output-dir",
                str(report_dir),
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return report_dir / "sra_project_manifest.tsv"

    def make_tool(self, directory: Path, name: str, body: str) -> Path:
        path = directory / name
        path.write_text("#!/usr/bin/env python3\n" + textwrap.dedent(body), encoding="utf-8")
        path.chmod(0o755)
        return path

    def make_fake_tools(self, root: Path) -> dict[str, Path]:
        tools = root / "fake-tools"
        tools.mkdir()
        prefetch = self.make_tool(
            tools,
            "prefetch",
            r'''
            import os
            import pathlib
            import sys
            args = sys.argv[1:]
            output = pathlib.Path(args[args.index("--output-directory") + 1])
            run = args[-1]
            target = output / run
            target.mkdir(parents=True)
            (target / f"{run}.sra").write_bytes(b"synthetic-sra")
            with open(os.environ["FAKE_TOOL_LOG"], "a") as log:
                log.write(f"prefetch {run}\n")
            ''',
        )
        validate = self.make_tool(
            tools,
            "vdb-validate",
            r'''
            import os
            import sys
            with open(os.environ["FAKE_TOOL_LOG"], "a") as log:
                log.write(f"validate {sys.argv[-1]}\n")
            ''',
        )
        fasterq = self.make_tool(
            tools,
            "fasterq-dump",
            r'''
            import os
            import pathlib
            import sys
            args = sys.argv[1:]
            output = pathlib.Path(args[args.index("--outdir") + 1])
            run = pathlib.Path(args[-1]).name
            output.mkdir(parents=True, exist_ok=True)
            for mate in (1, 2):
                (output / f"{run}_{mate}.fastq").write_text(
                    f"@{run}/{mate}\nACGT\n+\nIIII\n", encoding="ascii"
                )
            if os.environ.get("FAKE_SINGLETON") == "1":
                (output / f"{run}.fastq").write_text(
                    f"@{run}\nACGT\n+\nIIII\n", encoding="ascii"
                )
            with open(os.environ["FAKE_TOOL_LOG"], "a") as log:
                log.write(f"fasterq {run}\n")
            ''',
        )
        pigz = self.make_tool(
            tools,
            "pigz",
            r'''
            import gzip
            import pathlib
            import shutil
            import sys
            args = sys.argv[1:]
            decompress = "-d" in args or "-dc" in args
            to_stdout = "-c" in args or "-dc" in args
            positional = []
            skip = False
            for index, value in enumerate(args):
                if skip:
                    skip = False
                    continue
                if value in ("-p", "--processes"):
                    skip = True
                elif not value.startswith("-"):
                    positional.append(value)
            if decompress:
                for name in positional:
                    with gzip.open(name, "rb") as source:
                        shutil.copyfileobj(source, sys.stdout.buffer)
            elif to_stdout:
                with gzip.GzipFile(fileobj=sys.stdout.buffer, mode="wb", mtime=0) as target:
                    shutil.copyfileobj(sys.stdin.buffer, target)
            else:
                source_path = pathlib.Path(positional[-1])
                with source_path.open("rb") as source, gzip.GzipFile(
                    filename=str(source_path) + ".gz", mode="wb", mtime=0
                ) as target:
                    shutil.copyfileobj(source, target)
                source_path.unlink()
            ''',
        )
        return {
            "prefetch": prefetch,
            "validate": validate,
            "fasterq": fasterq,
            "pigz": pigz,
        }

    def acquisition_arguments(
        self, root: Path, manifest: Path, tools: dict[str, Path]
    ) -> list[str]:
        return [
            str(ACQUIRE),
            "--manifest",
            str(manifest),
            "--sample-id",
            "SAMN000001",
            "--output-dir",
            str(root / "output"),
            "--scratch-dir",
            str(root / "scratch"),
            "--prefetch-dir",
            str(root / "prefetch"),
            "--temp-dir",
            str(root / "temp"),
            "--threads",
            "2",
            "--prefetch-executable",
            tools["prefetch"].as_posix(),
            "--vdb-validate-executable",
            tools["validate"].as_posix(),
            "--fasterq-dump-executable",
            tools["fasterq"].as_posix(),
            "--pigz-executable",
            tools["pigz"].as_posix(),
        ]

    def test_dry_run_prints_sequential_plan_without_calling_tools(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.resolve_manifest(root)
            tools = self.make_fake_tools(root)
            result = run_command([*self.acquisition_arguments(root, manifest, tools), "--dry-run"])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertLess(result.stderr.index("SRR000001"), result.stderr.index("SRR000002"))
            self.assertIn("|", result.stderr)
            self.assertIn("SAMN000001_R1.fastq.gz", result.stderr)
            self.assertFalse((root / "scratch").exists())
            self.assertFalse((root / "output").exists())

    def test_injected_tools_merge_runs_and_cleanup_all_task_scratch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.resolve_manifest(root)
            tools = self.make_fake_tools(root)
            log = root / "tools.log"
            environment = os.environ.copy()
            environment["FAKE_TOOL_LOG"] = str(log)
            result = run_command(
                self.acquisition_arguments(root, manifest, tools), env=environment
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            output_1 = root / "output" / "SAMN000001_R1.fastq.gz"
            output_2 = root / "output" / "SAMN000001_R2.fastq.gz"
            with gzip.open(output_1, "rt", encoding="ascii") as handle:
                mate_1 = handle.read()
            with gzip.open(output_2, "rt", encoding="ascii") as handle:
                mate_2 = handle.read()
            self.assertEqual(mate_1.count("@SRR"), 2)
            self.assertEqual(mate_2.count("@SRR"), 2)
            self.assertLess(mate_1.index("SRR000001"), mate_1.index("SRR000002"))
            self.assertLess(mate_2.index("SRR000001"), mate_2.index("SRR000002"))
            commands = log.read_text().splitlines()
            self.assertEqual(
                [line.split()[0] for line in commands],
                ["prefetch", "validate", "fasterq", "prefetch", "validate", "fasterq"],
            )
            self.assertFalse((root / "scratch" / "SAMN000001.acquisition").exists())
            self.assertFalse((root / "prefetch" / "SAMN000001.prefetch").exists())
            self.assertFalse((root / "temp" / "SAMN000001.fasterq-temp").exists())
            self.assertEqual(list(root.rglob("*.fastq")), [])

    def test_singleton_output_is_rejected_and_uncompressed_files_are_removed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.resolve_manifest(root)
            tools = self.make_fake_tools(root)
            log = root / "tools.log"
            environment = os.environ.copy()
            environment["FAKE_TOOL_LOG"] = str(log)
            environment["FAKE_SINGLETON"] = "1"
            result = run_command(
                self.acquisition_arguments(root, manifest, tools), env=environment
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("unexpected/singleton", result.stderr)
            self.assertFalse((root / "output" / "SAMN000001_R1.fastq.gz").exists())
            self.assertFalse((root / "output" / "SAMN000001_R2.fastq.gz").exists())
            self.assertEqual(list(root.rglob("*.fastq")), [])


class CheckpointTests(unittest.TestCase):
    def make_manifest(self, root: Path) -> Path:
        resolved = root / "resolved"
        result = run_command(
            [
                str(RESOLVER),
                "PRJNA123456",
                "--runinfo-file",
                str(FIXTURES / "runinfo_valid.csv"),
                "--output-dir",
                str(resolved),
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        source = resolved / "sra_project_manifest.tsv"
        rows = [row for row in read_tsv(source) if row["sample_id"] == "SAMN000001"]
        manifest = root / "single_sample_manifest.tsv"
        with source.open("r", encoding="utf-8", newline="") as handle:
            fields = next(csv.reader(handle, delimiter="\t"))
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return manifest

    @staticmethod
    def make_pair(root: Path) -> tuple[Path, Path]:
        paths = (root / "nonhost_R1.fastq.gz", root / "nonhost_R2.fastq.gz")
        for mate, path in enumerate(paths, start=1):
            with gzip.open(path, "wt", encoding="ascii") as handle:
                handle.write(f"@read-a/{mate}\nACGT\n+\nIIII\n@read-b/{mate}\nTGCA\n+\nIIII\n")
        return paths

    @staticmethod
    def make_scientific_results(root: Path) -> tuple[Path, dict[str, Path]]:
        results = root / "results"
        artifacts = {
            "qc": results
            / "01_quality_control_and_filtering"
            / "fastqc"
            / "sample.fastqc.html",
            "final_mag": results
            / "02_mag_construction"
            / "final_catalog"
            / "final_catalog"
            / "MAG_1.fa",
            "provenance": results
            / "02_mag_construction"
            / "final_catalog"
            / "final_catalog.provenance.tsv",
            "quality": results
            / "02_mag_construction"
            / "final_catalog"
            / "final_catalog.quality.tsv",
            "checkm2": results
            / "02_mag_construction"
            / "final_catalog"
            / "evaluation"
            / "checkm2"
            / "final.checkm2.quality_report.tsv",
            "gunc": results
            / "02_mag_construction"
            / "final_catalog"
            / "evaluation"
            / "gunc"
            / "final.gunc.summary.tsv",
            "species": results
            / "02_mag_construction"
            / "final_catalog"
            / "species_95"
            / "final_species.representatives"
            / "MAG_1.fa",
            "gtdb": results
            / "03_taxonomic_classification_and_phylogenomics"
            / "gtdbtk"
            / "final.gtdbtk.bac120.summary.tsv",
            "tree": results
            / "03_taxonomic_classification_and_phylogenomics"
            / "phylophlan"
            / "final.phylophlan.tree.nwk",
            "annotations": results
            / "04_gene_prediction_and_functional_annotation"
            / "integrated"
            / "MAG_1"
            / "MAG_1.functional_annotations.tsv",
            "abundance": results
            / "05_mag_abundance_estimation"
            / "final_catalog.mag_abundance.long.tsv",
            "multiqc": results
            / "06_global_processing_evaluation"
            / "global_processing_evaluation.multiqc.html",
            "versions": results / "pipeline_info" / "software_versions.tsv",
        }
        for label, path in artifacts.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"synthetic {label}\n", encoding="utf-8")
        artifacts["abundance"].write_text(
            "sample\tmag_id\trelative_abundance_percent\tmean_coverage\t"
            "covered_fraction\tgenome_length\n"
            "SAMN000001\tMAG_1\t25.5\t3.25\t0.75\t2048\n",
            encoding="utf-8",
        )
        return results, artifacts

    def test_persist_reconcile_detect_tamper_and_cleanup_after_success(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.make_manifest(root)
            read_1, read_2 = self.make_pair(root)
            report = root / "sample.fastp.json"
            report.write_text("{}\n", encoding="utf-8")
            checkpoint = root / "checkpoint"
            record = root / "checkpoint_record.json"
            persisted = run_command(
                [
                    str(CHECKPOINTS), "persist",
                    "--run-manifest", str(manifest),
                    "--sample-id", "SAMN000001",
                    "--read-1", str(read_1),
                    "--read-2", str(read_2),
                    "--checkpoint-dir", str(checkpoint),
                    "--report", str(report),
                    "--output-record", str(record),
                ]
            )
            self.assertEqual(persisted.returncode, 0, persisted.stderr)
            owner = json.loads(
                (checkpoint / "sra_checkpoint_owner.json").read_text(encoding="utf-8")
            )
            self.assertEqual(owner["project_accession"], "PRJNA123456")
            self.assertEqual(
                owner["run_manifest_sha256"],
                hashlib.sha256(manifest.read_bytes()).hexdigest(),
            )
            checkpoint_manifest = root / "checkpoints.tsv"
            pending = root / "pending.tsv"
            status = root / "status.json"
            reconciled = run_command(
                [
                    str(CHECKPOINTS), "reconcile",
                    "--run-manifest", str(manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--output-manifest", str(checkpoint_manifest),
                    "--pending-output", str(pending),
                    "--status-output", str(status),
                    "--require-complete",
                ]
            )
            self.assertEqual(reconciled.returncode, 0, reconciled.stderr)
            self.assertEqual(json.loads(status.read_text())["complete_samples"], 1)
            checkpoint_rows = read_tsv(checkpoint_manifest)
            self.assertEqual(checkpoint_rows[0]["paired_fastq_records"], "2")
            validated_sample = run_command(
                [
                    str(CHECKPOINTS), "validate-sample",
                    "--run-manifest", str(manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--sample-id", "SAMN000001",
                ]
            )
            self.assertEqual(validated_sample.returncode, 0, validated_sample.stderr)
            self.assertEqual(json.loads(validated_sample.stdout)["sample_id"], "SAMN000001")
            persisted_read_1 = Path(checkpoint_rows[0]["read_1"])
            original = persisted_read_1.read_bytes()
            persisted_read_1.write_bytes(original + b"tamper")
            refused_sample = run_command(
                [
                    str(CHECKPOINTS), "validate-sample",
                    "--run-manifest", str(manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--sample-id", "SAMN000001",
                ]
            )
            self.assertEqual(refused_sample.returncode, 2)
            tampered = run_command(
                [
                    str(CHECKPOINTS), "reconcile",
                    "--run-manifest", str(manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--output-manifest", str(root / "tampered.tsv"),
                    "--pending-output", str(root / "tampered_pending.tsv"),
                    "--require-complete",
                ]
            )
            self.assertEqual(tampered.returncode, 2)
            persisted_read_1.write_bytes(original)

            persisted_report = checkpoint / "reports" / "SAMN000001" / report.name
            original_report = persisted_report.read_bytes()
            persisted_report.write_bytes(original_report + b"tamper")
            report_tampered = run_command(
                [
                    str(CHECKPOINTS), "reconcile",
                    "--run-manifest", str(manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--output-manifest", str(root / "report_tampered.tsv"),
                    "--pending-output", str(root / "report_tampered_pending.tsv"),
                    "--require-complete",
                ]
            )
            self.assertEqual(report_tampered.returncode, 2)
            report_pending_rows = read_tsv(root / "report_tampered_pending.tsv")
            self.assertEqual(len(report_pending_rows), 1)
            self.assertIn("report", report_pending_rows[0]["reason"].lower())
            persisted_report.write_bytes(original_report)

            success = root / "success.json"
            results, artifacts = self.make_scientific_results(root)
            multiqc = artifacts["multiqc"]
            versions = artifacts["versions"]
            abundance = artifacts["abundance"]

            def description(path: Path) -> dict[str, object]:
                return {
                    "path": str(path.resolve()),
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }

            success_payload = {
                "schema_version": 1,
                "status": "complete",
                "project_accession": "PRJNA123456",
                "checkpoint_manifest": description(checkpoint_manifest),
                "outputs": {
                    "multiqc_report": description(multiqc),
                    "software_versions": description(versions),
                    "mag_abundance": description(abundance),
                },
            }
            success.write_text(
                json.dumps(success_payload),
                encoding="utf-8",
            )

            unsealed_cleanup = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(unsealed_cleanup.returncode, 2)
            self.assertIn("no sealed scientific output inventory", unsealed_cleanup.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())

            valid_abundance = abundance.read_text(encoding="utf-8")
            abundance_header = valid_abundance.splitlines()[0] + "\n"
            invalid_abundance_cases = (
                ("header only", abundance_header, "contains no data rows"),
                (
                    "foreign sample",
                    abundance_header
                    + "SAMN999999\tMAG_1\t25.5\t3.25\t0.75\t2048\n",
                    "row outside the checkpoint-sample",
                ),
                (
                    "duplicate pair",
                    valid_abundance + valid_abundance.splitlines()[1] + "\n",
                    "duplicate sample/MAG pair",
                ),
                (
                    "non-finite metric",
                    valid_abundance.replace("\t3.25\t", "\tNaN\t"),
                    "non-finite mean_coverage",
                ),
            )
            for label, invalid_content, expected_error in invalid_abundance_cases:
                with self.subTest(invalid_abundance=label):
                    abundance.write_text(invalid_content, encoding="utf-8")
                    invalid_payload = json.loads(json.dumps(success_payload))
                    invalid_payload["outputs"]["mag_abundance"] = description(abundance)
                    success.write_text(json.dumps(invalid_payload), encoding="utf-8")
                    refused_semantics = run_command(
                        [
                            str(CHECKPOINTS), "seal-global",
                            "--success-marker", str(success),
                            "--checkpoint-manifest", str(checkpoint_manifest),
                            "--results-dir", str(results),
                            "--project-accession", "PRJNA123456",
                        ]
                    )
                    self.assertEqual(refused_semantics.returncode, 2)
                    self.assertIn(expected_error, refused_semantics.stderr)
                    self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
                    self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())

            abundance.write_text(valid_abundance, encoding="utf-8")
            success_payload["outputs"]["mag_abundance"] = description(abundance)

            missing_sample_manifest = root / "missing_sample_abundance.checkpoints.tsv"
            second_checkpoint_row = dict(checkpoint_rows[0])
            second_checkpoint_row["sample_order"] = "2"
            second_checkpoint_row["sample_id"] = "SAMN000002"
            with missing_sample_manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(checkpoint_rows[0]),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows([checkpoint_rows[0], second_checkpoint_row])
            missing_sample_payload = json.loads(json.dumps(success_payload))
            missing_sample_payload["checkpoint_manifest"] = description(
                missing_sample_manifest
            )
            success.write_text(json.dumps(missing_sample_payload), encoding="utf-8")
            refused_missing_sample = run_command(
                [
                    str(CHECKPOINTS), "seal-global",
                    "--success-marker", str(success),
                    "--checkpoint-manifest", str(missing_sample_manifest),
                    "--results-dir", str(results),
                    "--project-accession", "PRJNA123456",
                ]
            )
            self.assertEqual(refused_missing_sample.returncode, 2)
            self.assertIn("exact checkpoint-sample x final-MAG matrix", refused_missing_sample.stderr)

            extra_mag = artifacts["final_mag"].with_name("MAG_2.fa")
            extra_annotation = artifacts["annotations"].parent.parent / "MAG_2" / "MAG_2.functional_annotations.tsv"
            extra_mag.write_text(">MAG_2\nACGT\n", encoding="utf-8")
            extra_annotation.parent.mkdir(parents=True, exist_ok=True)
            extra_annotation.write_text("synthetic annotations\n", encoding="utf-8")
            success.write_text(json.dumps(success_payload), encoding="utf-8")
            refused_missing_mag = run_command(
                [
                    str(CHECKPOINTS), "seal-global",
                    "--success-marker", str(success),
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--results-dir", str(results),
                    "--project-accession", "PRJNA123456",
                ]
            )
            self.assertEqual(refused_missing_mag.returncode, 2)
            self.assertIn("exact checkpoint-sample x final-MAG matrix", refused_missing_mag.stderr)
            extra_mag.unlink()
            extra_annotation.unlink()
            extra_annotation.parent.rmdir()
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())
            success.write_text(json.dumps(success_payload), encoding="utf-8")

            sealed = run_command(
                [
                    str(CHECKPOINTS), "seal-global",
                    "--success-marker", str(success),
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--results-dir", str(results),
                    "--project-accession", "PRJNA123456",
                ]
            )
            self.assertEqual(sealed.returncode, 0, sealed.stderr)
            success_payload = json.loads(success.read_text(encoding="utf-8"))
            scientific_outputs = success_payload["scientific_outputs"]
            self.assertEqual(scientific_outputs["file_count"], len(artifacts))
            self.assertEqual(
                [item["relative_path"] for item in scientific_outputs["files"]],
                sorted(item["relative_path"] for item in scientific_outputs["files"]),
            )
            self.assertEqual(
                scientific_outputs["total_bytes"],
                sum(item["bytes"] for item in scientific_outputs["files"]),
            )
            self.assertEqual(
                scientific_outputs["required_artifacts"]["final_mags"],
                ["02_mag_construction/final_catalog/final_catalog/MAG_1.fa"],
            )
            marker_after_first_seal = success.read_bytes()
            sealed_again = run_command(
                [
                    str(CHECKPOINTS), "seal-global",
                    "--success-marker", str(success),
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--results-dir", str(results),
                    "--project-accession", "PRJNA123456",
                ]
            )
            self.assertEqual(sealed_again.returncode, 0, sealed_again.stderr)
            self.assertIn("already sealed and validated", sealed_again.stdout)
            self.assertEqual(success.read_bytes(), marker_after_first_seal)

            retained = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                    "--keep",
                ]
            )
            self.assertEqual(retained.returncode, 0, retained.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())

            missing_read = Path(checkpoint_rows[0]["read_2"])
            hidden_read = missing_read.with_suffix(".temporarily-hidden")
            missing_read.replace(hidden_read)
            refused_missing = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_missing.returncode, 2)
            self.assertIn("checkpoint read is missing", refused_missing.stderr)
            hidden_read.replace(missing_read)

            success_payload["checkpoint_manifest"]["sha256"] = "0" * 64
            success.write_text(json.dumps(success_payload), encoding="utf-8")
            refused_manifest = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_manifest.returncode, 2)
            self.assertIn("manifest SHA-256 differs", refused_manifest.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())

            success_payload["checkpoint_manifest"] = description(checkpoint_manifest)
            success.write_text(json.dumps(success_payload), encoding="utf-8")

            original_provenance = artifacts["provenance"].read_bytes()
            artifacts["provenance"].write_bytes(original_provenance + b"tamper")
            marker_before_refused_reseal = success.read_bytes()
            refused_reseal = run_command(
                [
                    str(CHECKPOINTS), "seal-global",
                    "--success-marker", str(success),
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--results-dir", str(results),
                    "--project-accession", "PRJNA123456",
                ]
            )
            self.assertEqual(refused_reseal.returncode, 2)
            self.assertIn("inventory differs", refused_reseal.stderr)
            self.assertEqual(success.read_bytes(), marker_before_refused_reseal)
            refused_scientific_tamper = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_scientific_tamper.returncode, 2)
            self.assertIn("inventory differs", refused_scientific_tamper.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())
            artifacts["provenance"].write_bytes(original_provenance)

            extra_artifact = (
                results / "01_quality_control_and_filtering" / "unexpected.tsv"
            )
            extra_artifact.write_text("unexpected\n", encoding="utf-8")
            refused_extra = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_extra.returncode, 2)
            self.assertIn("inventory differs", refused_extra.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())
            extra_artifact.unlink()

            hidden_qc = artifacts["qc"].with_suffix(".temporarily-hidden")
            artifacts["qc"].replace(hidden_qc)
            refused_missing_scientific = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_missing_scientific.returncode, 2)
            self.assertIn("inventory differs", refused_missing_scientific.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())
            hidden_qc.replace(artifacts["qc"])

            original_abundance = abundance.read_bytes()
            abundance.write_bytes(original_abundance + b"tamper")
            refused_output = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_output.returncode, 2)
            self.assertIn("durable global output size differs", refused_output.stderr)
            self.assertTrue(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertTrue(Path(checkpoint_rows[0]["read_2"]).exists())
            abundance.write_bytes(original_abundance)

            # Simulate an interruption in the only unavoidably ambiguous
            # window: the durable deletion plan was committed and the first
            # unlink succeeded, but removed_files was not updated yet.  A
            # later invocation must trust the already validated plan, finish
            # the remaining removals, and atomically mark the journal complete.
            planned_files = [
                {
                    "sample_id": checkpoint_rows[0]["sample_id"],
                    "path": checkpoint_rows[0][field],
                    "bytes": int(checkpoint_rows[0][f"{field}_bytes"]),
                    "sha256": checkpoint_rows[0][f"{field}_sha256"],
                }
                for field in ("read_1", "read_2")
            ]
            cleanup_record = checkpoint / "sra_checkpoint_cleanup.json"
            cleanup_record.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "project_accession": "PRJNA123456",
                        "status": "in_progress",
                        "started_at_utc": "2026-01-01T00:00:00Z",
                        "updated_at_utc": "2026-01-01T00:00:00Z",
                        "checkpoint_manifest": success_payload["checkpoint_manifest"],
                        "planned_files": planned_files,
                        "removed_files": [],
                        "removed_bytes": 0,
                        "reports_and_provenance_retained": True,
                    }
                ),
                encoding="utf-8",
            )
            Path(planned_files[0]["path"]).unlink()

            refused_keep_after_start = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                    "--keep",
                ]
            )
            self.assertEqual(refused_keep_after_start.returncode, 2)
            self.assertIn("already started", refused_keep_after_start.stderr)

            cleaned = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(cleaned.returncode, 0, cleaned.stderr)
            self.assertFalse(Path(checkpoint_rows[0]["read_1"]).exists())
            self.assertFalse(Path(checkpoint_rows[0]["read_2"]).exists())
            self.assertTrue((checkpoint / "records" / "SAMN000001.checkpoint.json").exists())
            self.assertTrue((checkpoint / "reports" / "SAMN000001" / report.name).exists())
            completed_cleanup = json.loads(cleanup_record.read_text(encoding="utf-8"))
            self.assertEqual(completed_cleanup["status"], "complete")
            self.assertEqual(completed_cleanup["planned_files"], planned_files)
            self.assertEqual(completed_cleanup["removed_files"], planned_files)
            self.assertEqual(
                completed_cleanup["removed_bytes"],
                sum(item["bytes"] for item in planned_files),
            )

            tampered_cleanup = dict(completed_cleanup)
            tampered_cleanup["planned_files"] = list(reversed(planned_files))
            cleanup_record.write_text(json.dumps(tampered_cleanup), encoding="utf-8")
            refused_tampered_journal = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(refused_tampered_journal.returncode, 2)
            self.assertIn("cleanup plan does not match", refused_tampered_journal.stderr)
            cleanup_record.write_text(json.dumps(completed_cleanup), encoding="utf-8")

            cleaned_again = run_command(
                [
                    str(CHECKPOINTS), "cleanup",
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--checkpoint-dir", str(checkpoint),
                    "--success-marker", str(success),
                ]
            )
            self.assertEqual(cleaned_again.returncode, 0, cleaned_again.stderr)
            self.assertIn("already completed", cleaned_again.stdout)

    def test_seal_global_rejects_symbolic_links(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            checkpoint_manifest = root / "checkpoints.tsv"
            checkpoint_manifest.write_text("synthetic checkpoint manifest\n", encoding="utf-8")
            results, artifacts = self.make_scientific_results(root)
            linked = (
                results / "01_quality_control_and_filtering" / "linked-artifact.tsv"
            )
            try:
                linked.symlink_to(artifacts["provenance"])
            except (NotImplementedError, OSError) as exc:
                self.skipTest(f"symbolic links are unavailable in this environment: {exc}")

            def description(path: Path) -> dict[str, object]:
                return {
                    "path": str(path.resolve()),
                    "bytes": path.stat().st_size,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }

            success = root / "success.json"
            success.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "status": "complete",
                        "project_accession": "PRJNA123456",
                        "checkpoint_manifest": description(checkpoint_manifest),
                        "outputs": {
                            "multiqc_report": description(artifacts["multiqc"]),
                            "software_versions": description(artifacts["versions"]),
                            "mag_abundance": description(artifacts["abundance"]),
                        },
                    }
                ),
                encoding="utf-8",
            )
            refused = run_command(
                [
                    str(CHECKPOINTS), "seal-global",
                    "--success-marker", str(success),
                    "--checkpoint-manifest", str(checkpoint_manifest),
                    "--results-dir", str(results),
                    "--project-accession", "PRJNA123456",
                ]
            )
            self.assertEqual(refused.returncode, 2)
            self.assertIn("refuses symbolic link", refused.stderr)
            self.assertNotIn(
                "scientific_outputs",
                json.loads(success.read_text(encoding="utf-8")),
            )

    def test_checkpoint_root_is_owned_by_one_exact_frozen_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.make_manifest(root)
            read_1, read_2 = self.make_pair(root)
            report = root / "sample.fastp.json"
            report.write_text('{"cohort":"original"}\n', encoding="utf-8")
            checkpoint = root / "checkpoint"

            base_arguments = [
                str(CHECKPOINTS),
                "persist",
                "--sample-id",
                "SAMN000001",
                "--read-1",
                str(read_1),
                "--read-2",
                str(read_2),
                "--checkpoint-dir",
                str(checkpoint),
                "--report",
                str(report),
            ]
            first = run_command(
                [*base_arguments, "--run-manifest", str(manifest)]
            )
            self.assertEqual(first.returncode, 0, first.stderr)
            protected_paths = [
                checkpoint / "reads" / "SAMN000001_host_removed_R1.fastq.gz",
                checkpoint / "reads" / "SAMN000001_host_removed_R2.fastq.gz",
                checkpoint / "reports" / "SAMN000001" / report.name,
                checkpoint / "records" / "SAMN000001.checkpoint.json",
            ]
            protected_content = {path: path.read_bytes() for path in protected_paths}

            for mate, path in enumerate((read_1, read_2), start=1):
                with gzip.open(path, "wt", encoding="ascii") as handle:
                    handle.write(f"@replacement/{mate}\nAAAA\n+\nIIII\n")
            report.write_text('{"cohort":"replacement"}\n', encoding="utf-8")
            replayed = run_command(
                [*base_arguments, "--run-manifest", str(manifest)]
            )
            self.assertEqual(replayed.returncode, 0, replayed.stderr)
            self.assertIn("immutable checkpoint already exists", replayed.stdout)
            for path, original in protected_content.items():
                self.assertEqual(path.read_bytes(), original)

            foreign_manifest = root / "foreign_project_manifest.tsv"
            foreign_manifest.write_text(
                manifest.read_text(encoding="utf-8").replace(
                    "PRJNA123456", "PRJNA999999"
                ),
                encoding="utf-8",
            )
            foreign = run_command(
                [*base_arguments, "--run-manifest", str(foreign_manifest)]
            )
            self.assertEqual(foreign.returncode, 2)
            self.assertIn("owned by BioProject", foreign.stderr)

            changed_manifest = root / "changed_cohort_manifest.tsv"
            changed_rows = read_tsv(manifest)
            changed_rows[0]["model"] = "DIFFERENT_MODEL"
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                fields = next(csv.reader(handle, delimiter="\t"))
            with changed_manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(changed_rows)
            changed = run_command(
                [*base_arguments, "--run-manifest", str(changed_manifest)]
            )
            self.assertEqual(changed.returncode, 2)
            self.assertIn("different frozen manifest", changed.stderr)
            for path, original in protected_content.items():
                self.assertEqual(path.read_bytes(), original)

            occupied = root / "occupied_checkpoint"
            occupied.mkdir()
            unrelated = occupied / "do-not-overwrite.txt"
            unrelated.write_text("unrelated\n", encoding="utf-8")
            unclaimed = run_command(
                [
                    *base_arguments,
                    "--run-manifest",
                    str(manifest),
                    "--checkpoint-dir",
                    str(occupied),
                ]
            )
            self.assertEqual(unclaimed.returncode, 2)
            self.assertIn("unclaimed checkpoint root is not empty", unclaimed.stderr)
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "unrelated\n")
            self.assertFalse((occupied / "sra_checkpoint_owner.json").exists())

if __name__ == "__main__":
    unittest.main()
