from __future__ import annotations

import csv
import importlib.util
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "bin" / "monitor_storage.py"
SPEC = importlib.util.spec_from_file_location("monitor_storage", SCRIPT)
assert SPEC and SPEC.loader
monitor_storage = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = monitor_storage
SPEC.loader.exec_module(monitor_storage)


class StorageMonitorTests(unittest.TestCase):
    def test_once_samples_categories_task_peaks_and_never_follows_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work = root / "work"
            task = work / "ab" / "0123456789abcdef"
            task.mkdir(parents=True)
            (task / ".command.sh").write_text("echo test\n", encoding="utf-8")
            (task / "output.bin").write_bytes(b"x" * 4096)
            results = root / "results"
            results.mkdir()
            (results / "report.txt").write_text("result\n", encoding="utf-8")
            database = root / "database"
            database.mkdir()
            (database / "index.bin").write_bytes(b"d" * 8192)
            checkpoint = root / "checkpoint"
            checkpoint.mkdir()
            (checkpoint / "state.json").write_text("{}\n", encoding="utf-8")

            external = root / "external"
            external.mkdir()
            (external / "large.bin").write_bytes(b"z" * (2 * 1024 * 1024))
            symlink_created = False
            try:
                (task / "external-link").symlink_to(external, target_is_directory=True)
                symlink_created = True
            except OSError:
                # Unprivileged Windows runners cannot create symlinks.  The
                # category and peak assertions still run; POSIX/Developer Mode
                # runners additionally exercise the no-follow contract below.
                pass

            output = root / "telemetry"
            result = monitor_storage.main(
                [
                    "--output-dir",
                    str(output),
                    "--work-dir",
                    str(work),
                    "--checkpoint-dir",
                    str(checkpoint),
                    "--results-dir",
                    str(results),
                    "--database-dir",
                    str(database),
                    "--once",
                ]
            )
            self.assertEqual(result, 0)

            with (output / "storage_usage_timeseries.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["measurement_complete"], "true")
            self.assertGreater(int(rows[0]["database_bytes"]), 0)
            self.assertEqual(
                int(rows[0]["total_measured_bytes"]),
                int(rows[0]["total_dynamic_bytes"]) + int(rows[0]["database_bytes"]),
            )

            with (output / "task_workdir_peaks.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                peak_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(peak_rows), 1)
            self.assertEqual(peak_rows[0]["workdir"], os.path.abspath(task))
            task_size = int(peak_rows[0]["peak_work_bytes"])
            external_size = monitor_storage.scan_tree(external).allocated_bytes
            self.assertIsNotNone(external_size)
            if symlink_created:
                self.assertLess(task_size, int(external_size))

    def test_database_is_measured_once_and_task_peaks_survive_samples(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = root / "database"
            database.mkdir()
            (database / "before.bin").write_bytes(b"a" * 4096)
            work = root / "work"
            task = work / "ff" / "task"
            task.mkdir(parents=True)
            (task / "one.bin").write_bytes(b"b" * 4096)
            output = root / "telemetry"
            instance = monitor_storage.StorageMonitor(
                paths={
                    "work": work,
                    "checkpoint": None,
                    "sra_cache": None,
                    "sra_scratch": None,
                    "results": None,
                    "database": database,
                },
                timeseries_path=output / "storage_usage_timeseries.tsv",
                task_peaks_path=output / "task_workdir_peaks.tsv",
            )
            instance.measure_database()
            instance.sample()
            (database / "after.bin").write_bytes(b"c" * (1024 * 1024))
            (task / "two.bin").write_bytes(b"d" * (64 * 1024))
            instance.sample()

            with instance.timeseries_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["database_bytes"], rows[1]["database_bytes"])
            with instance.task_peaks_path.open(encoding="utf-8", newline="") as handle:
                peaks = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(peaks[0]["samples"], "2")
            self.assertGreater(int(peaks[0]["peak_work_bytes"]), 64 * 1024)

    def test_request_ack_forces_peak_persistence_before_cleanup(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work = root / "work"
            task = work / "aa" / "task"
            task.mkdir(parents=True)
            output = root / "telemetry"
            request = root / "sample-request"
            stop = root / "stop"
            process = subprocess.Popen(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--output-dir",
                    str(output),
                    "--work-dir",
                    str(work),
                    "--sample-request-file",
                    str(request),
                    "--stop-file",
                    str(stop),
                    "--interval-seconds",
                    "60",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            timeseries = output / "storage_usage_timeseries.tsv"
            deadline = time.monotonic() + 5
            while not timeseries.exists() and time.monotonic() < deadline:
                time.sleep(0.02)
            self.assertTrue(timeseries.exists())

            payload = task / "transient.bin"
            payload.write_bytes(b"x" * (128 * 1024))
            request.write_text("force-before-cleanup\n", encoding="utf-8")
            deadline = time.monotonic() + 5
            while request.exists() and time.monotonic() < deadline:
                time.sleep(0.02)
            self.assertFalse(request.exists(), "monitor did not acknowledge forced sample")
            shutil.rmtree(task)
            if os.name == "nt":
                stop.touch()
            else:
                process.send_signal(signal.SIGTERM)
            _stdout, stderr = process.communicate(timeout=5)
            self.assertEqual(process.returncode, 0, stderr)

            with (output / "task_workdir_peaks.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                peaks = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(peaks), 1)
            self.assertGreaterEqual(int(peaks[0]["peak_work_bytes"]), 128 * 1024)

    def test_sigterm_where_supported_and_stop_file_exit_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work = root / "work"
            work.mkdir()
            output = root / "telemetry"
            first_stop = root / "first-stop"
            process = subprocess.Popen(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--output-dir",
                    str(output),
                    "--work-dir",
                    str(work),
                    "--stop-file",
                    str(first_stop),
                    "--interval-seconds",
                    "10",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            deadline = time.monotonic() + 5
            timeseries = output / "storage_usage_timeseries.tsv"
            while not timeseries.exists() and time.monotonic() < deadline:
                time.sleep(0.02)
            self.assertTrue(timeseries.exists())
            if os.name == "nt":
                first_stop.touch()
            else:
                process.send_signal(signal.SIGTERM)
            _stdout, stderr = process.communicate(timeout=5)
            self.assertEqual(process.returncode, 0, stderr)

            stop = root / "stop"
            stop.touch()
            stopped = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--output-dir",
                    str(root / "stopped-output"),
                    "--work-dir",
                    str(work),
                    "--stop-file",
                    str(stop),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=5,
            )
            self.assertEqual(stopped.returncode, 0, stopped.stderr)


if __name__ == "__main__":
    unittest.main()
