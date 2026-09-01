from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "bin" / "collect_slurm_accounting.py"
SPEC = importlib.util.spec_from_file_location("collect_slurm_accounting", SCRIPT)
assert SPEC and SPEC.loader
accounting = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = accounting
SPEC.loader.exec_module(accounting)


class SlurmAccountingTests(unittest.TestCase):
    def test_collect_uses_timeout_and_preserves_records(self) -> None:
        completed = subprocess.CompletedProcess(
            args=["sacct"],
            returncode=0,
            stdout="123|COMPLETED|60|00:01:00|1G|4|cpu=4,gres/gpu=1|cpu=4\n",
            stderr="",
        )
        with mock.patch.object(accounting.subprocess, "run", return_value=completed) as run:
            rows = accounting.collect("sacct", {"123": {"invocation-a"}}, 7.5)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["invocation_id"], "invocation-a")
        self.assertIn("gres/gpu=1", rows[0]["allocated_tres"])
        self.assertGreater(run.call_args.kwargs["timeout"], 0)
        self.assertLessEqual(run.call_args.kwargs["timeout"], 7.5)

    def test_timeout_budget_is_shared_across_sacct_batches(self) -> None:
        completed = subprocess.CompletedProcess(
            args=["sacct"],
            returncode=0,
            stdout="123|COMPLETED|60|00:01:00|1G|4|cpu=4|cpu=4\n",
            stderr="",
        )
        with mock.patch.object(
            accounting,
            "chunks",
            return_value=[["123"], ["124"]],
        ), mock.patch.object(
            accounting.time,
            "monotonic",
            side_effect=[100.0, 100.0, 103.0],
        ), mock.patch.object(
            accounting.subprocess,
            "run",
            return_value=completed,
        ) as run:
            accounting.collect("sacct", {"123": {"invocation-a"}}, 7.5)
        self.assertEqual(
            [call.kwargs["timeout"] for call in run.call_args_list],
            [7.5, 4.5],
        )

    def test_timeout_writes_explicit_unavailable_status(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            registry = root / "registry.tsv"
            registry.write_text("invocation_id\ttrace_path\n", encoding="utf-8")
            output = root / "slurm_accounting.tsv"
            with mock.patch.object(
                accounting,
                "registered_jobs",
                return_value={"123": {"invocation-a"}},
            ), mock.patch.object(
                accounting,
                "collect",
                side_effect=subprocess.TimeoutExpired("sacct", 0.1),
            ):
                result = accounting.main(
                    [
                        "--registry",
                        str(registry),
                        "--output",
                        str(output),
                        "--timeout-seconds",
                        "0.1",
                    ]
                )
            self.assertEqual(result, 0)
            with output.open(encoding="utf-8", newline="") as handle:
                self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), [])
            status = json.loads(
                output.with_suffix(".status.json").read_text(encoding="utf-8")
            )
            self.assertFalse(status["available"])
            self.assertEqual(status["state"], "error")
            self.assertIn("TimeoutExpired", status["error"])


if __name__ == "__main__":
    unittest.main()
