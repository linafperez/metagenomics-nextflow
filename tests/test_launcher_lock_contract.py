from __future__ import annotations

import re
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]


class LauncherLockContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.launcher = (REPOSITORY / "metagenomics_pipeline.sh").read_text(
            encoding="utf-8"
        )

    def function_body(self, name: str) -> str:
        match = re.search(
            rf"(?ms)^{re.escape(name)}\(\) \{{\n(.*?)^\}}$", self.launcher
        )
        self.assertIsNotNone(match, f"missing shell function {name}")
        return match.group(1)

    def test_locks_are_atomic_fail_closed_and_safely_released(self) -> None:
        acquire = self.function_body("acquire_run_lock")
        release = self.function_body("release_run_locks")

        self.assertIn('mkdir -- "${lock_dir}"', acquire)
        self.assertIn("run lock already exists", acquire)
        self.assertIn("This is fail-closed", acquire)
        self.assertLess(
            acquire.index("printf '%s\\n' \"${run_lock_token}\" > \"${claim_file}\""),
            acquire.index('run_lock_dirs+=("${lock_dir}")'),
        )
        self.assertIn('rmdir -- "${lock_dir}"', acquire)
        for field in (
            "lock_token",
            "scope",
            "pid",
            "hostname",
            "started_at_utc",
            "results_root",
            "checkpoint_root",
        ):
            self.assertIn(f"printf '{field}\\t", acquire)

        self.assertIn('for ((index = ${#run_lock_dirs[@]} - 1;', release)
        self.assertIn('rmdir -- "${lock_dir}"', release)
        self.assertNotIn("rm -rf", release)
        self.assertIn("unexpected content", release)
        self.assertIn("ownership marker is missing or changed", release)
        self.assertIn('[[ ! -f "${claim_file}" || -L "${claim_file}" ]]', release)
        self.assertIn('[[ "${claim_value}" != "${run_lock_token}" ]]', release)

    def test_results_then_checkpoint_lock_precede_execution(self) -> None:
        results_lock = self.launcher.index(
            'acquire_run_lock "${resource_root}/.metagenomics_run.lock" results'
        )
        checkpoint_lock = self.launcher.index(
            'acquire_run_lock "${sra_checkpoint_dir}.metagenomics_run.lock" checkpoint'
        )
        checkpoint_root_creation = self.launcher.index(
            'mkdir -p -- "${sra_checkpoint_dir}"'
        )
        sra_telemetry = self.launcher.index(
            "initialize_project_telemetry", checkpoint_lock
        )
        storage_monitor = self.launcher.index(
            'start_storage_monitor "${work_dir}" "${sra_checkpoint_dir}"'
        )
        scientific_execution = self.launcher.index(
            "run_nextflow discovery sra-discovery", checkpoint_lock
        )

        self.assertLess(results_lock, checkpoint_lock)
        self.assertLess(checkpoint_lock, sra_telemetry)
        self.assertLess(checkpoint_lock, checkpoint_root_creation)
        self.assertLess(checkpoint_lock, storage_monitor)
        self.assertLess(checkpoint_lock, scientific_execution)

    def test_dry_run_does_not_acquire_or_create_locks(self) -> None:
        results_lock = self.launcher.index(
            'acquire_run_lock "${resource_root}/.metagenomics_run.lock" results'
        )
        checkpoint_lock = self.launcher.index(
            'acquire_run_lock "${sra_checkpoint_dir}.metagenomics_run.lock" checkpoint'
        )

        results_guard = self.launcher.rfind(
            'if [[ "${dry_run}" == false ]]; then', 0, results_lock
        )
        results_guard_end = self.launcher.index("\nfi", results_lock)
        checkpoint_guard = self.launcher.rfind(
            'if [[ "${dry_run}" == false ]]; then', 0, checkpoint_lock
        )
        checkpoint_guard_end = self.launcher.index("\nfi", checkpoint_lock)

        self.assertLess(results_guard, results_lock, results_guard_end)
        self.assertLess(checkpoint_guard, checkpoint_lock, checkpoint_guard_end)
        dry_plan = self.launcher.index(
            'if [[ "${dry_run}" == true ]]; then', checkpoint_lock
        )
        self.assertLess(results_lock, dry_plan)
        self.assertLess(checkpoint_lock, dry_plan)

    def test_exit_finalization_releases_both_locks(self) -> None:
        finish = self.function_body("finish_telemetry")
        self.assertIn("release_run_locks", finish)
        final_release = finish.rindex("release_run_locks")
        self.assertLess(finish.index("summarize_resources"), final_release)
        self.assertLess(final_release, finish.index('exit "${requested_exit}"'))
        self.assertIn("trap 'release_run_locks' EXIT", self.launcher)
        self.assertIn("trap 'finish_telemetry $?' EXIT", self.launcher)


if __name__ == "__main__":
    unittest.main()
