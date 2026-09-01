from __future__ import annotations

import json
import re
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]


class StaticProductionContractTests(unittest.TestCase):
    def read(self, relative: str) -> str:
        return (REPOSITORY / relative).read_text(encoding="utf-8")

    def configured_process_names(self) -> set[str]:
        """Return module process names plus aliases visible to configuration selectors."""
        names: set[str] = set()
        for path in (REPOSITORY / "modules").rglob("main.nf"):
            names.update(
                re.findall(
                    r"(?m)^process\s+([A-Z][A-Z0-9_]*)\s*\{",
                    path.read_text(encoding="utf-8"),
                )
            )
        for root in (REPOSITORY / "workflows", REPOSITORY / "subworkflows"):
            for path in root.rglob("*.nf"):
                source = path.read_text(encoding="utf-8")
                for include in re.findall(r"include\s*\{([^}]+)\}", source):
                    for item in include.split(";"):
                        item = item.strip()
                        alias = re.search(r"\bas\s+([A-Z][A-Z0-9_]*)$", item)
                        names.add(alias.group(1) if alias else item)
        return names

    def test_schema_and_mutually_exclusive_dispatch(self) -> None:
        schema = json.loads(self.read("nextflow_schema.json"))
        self.assertEqual(len(schema["oneOf"]), 2)
        self.assertEqual(schema["oneOf"][0]["required"], ["input"])
        self.assertEqual(schema["oneOf"][0]["properties"]["sraProject"]["type"], "null")
        self.assertEqual(schema["oneOf"][1]["required"], ["sraProject"])
        self.assertEqual(schema["oneOf"][1]["properties"]["input"]["type"], "null")
        launcher = self.read("metagenomics_pipeline.sh")
        self.assertIn("--input and --sra-project are mutually exclusive", launcher)
        main = self.read("main.nf")
        for stage in ("local", "sra-discovery", "sra-checkpoints", "sra-preprocess", "sra-global"):
            self.assertIn(stage, main)

    def test_schema_covers_every_declared_parameter(self) -> None:
        schema = json.loads(self.read("nextflow_schema.json"))
        config = self.read("nextflow.config")
        params_match = re.search(r"(?ms)^params\s*\{(.*?)^\}", config)
        self.assertIsNotNone(params_match)
        configured = set(
            re.findall(
                r"(?m)^\s{4}([A-Za-z_][A-Za-z0-9_]*)\s*=",
                params_match.group(1),
            )
        )
        self.assertEqual(configured, set(schema["properties"]))
        self.assertEqual(schema["properties"]["gpuAccelerators"]["minimum"], 1)
        self.assertEqual(schema["properties"]["gpuAccelerators"]["maximum"], 1)

    def test_all_exact_process_selectors_resolve_to_a_process_or_alias(self) -> None:
        known_names = self.configured_process_names()
        for config_name in ("resources.config", "modules.config", "gpu.config"):
            source = self.read(f"conf/{config_name}")
            raw_selectors: list[str] = []
            for match in re.finditer(
                r"(?m)^\s*withName:\s*(?:'([^']+)'|([A-Z][A-Z0-9_]*))\s*\{",
                source,
            ):
                selector = match.group(1) or match.group(2)
                raw_selectors.append(selector)
                for process_name in selector.split("|"):
                    self.assertIn(
                        process_name,
                        known_names,
                        f"ghost selector {process_name!r} in conf/{config_name}",
                    )
            self.assertEqual(
                len(raw_selectors),
                len(set(raw_selectors)),
                f"duplicate exact withName block in conf/{config_name}",
            )

        resource_regexes = set(
            re.findall(r"withName:\s*/([^/]+)/\s*\{", self.read("conf/resources.config"))
        )
        self.assertEqual(
            resource_regexes,
            {
                ".*SPADES.*:VAMB",
                ".*SPADES.*:DREP_(ANI_99|SPECIES_95)",
            },
        )

    def test_resource_blocks_use_valid_units_and_preserve_metaspades_request(self) -> None:
        resources = self.read("conf/resources.config")
        blocks = re.findall(r"(?ms)^\s*withName:[^{]+\{\s*(.*?)^\s*\}", resources)
        self.assertGreater(len(blocks), 0)
        for block in blocks:
            self.assertRegex(block, r"(?m)^\s*cpus\s*=\s*[1-9][0-9]*\s*$")
            self.assertRegex(block, r"(?m)^\s*memory\s*=\s*[1-9][0-9]*(?:\.[0-9]+)?\.(?:KB|MB|GB|TB)\s*$")
            self.assertRegex(block, r"(?m)^\s*time\s*=\s*[1-9][0-9]*(?:\.[0-9]+)?\.(?:m|h|d)\s*$")
        spades = re.search(r"withName:\s*SPADES\s*\{([^}]+)\}", resources)
        self.assertIsNotNone(spades)
        self.assertRegex(spades.group(1), r"cpus\s*=\s*32")
        self.assertRegex(spades.group(1), r"memory\s*=\s*1800\.GB")
        self.assertRegex(spades.group(1), r"time\s*=\s*14\.d")

    def test_sra_process_selectors_are_not_ghosts(self) -> None:
        declarations = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (REPOSITORY / "modules").rglob("main.nf")
        )
        configs = self.read("conf/resources.config") + self.read("conf/modules.config")
        for process in (
            "RESOLVE_SRA_PROJECT",
            "VALIDATE_SRA_PROJECT",
            "SRA_ACQUIRE",
            "PERSIST_SRA_CHECKPOINT",
            "CHECK_SRA_CHECKPOINTS",
            "FINALIZE_SRA_GLOBAL_RUN",
        ):
            self.assertRegex(declarations, rf"(?m)^process\s+{process}\s*\{{")
            self.assertIn(process, configs)

    def test_compression_chain_and_checkpoint_contract(self) -> None:
        self.assertIn(".fastq.gz", self.read("modules/core/fastp/main.nf"))
        bowtie = self.read("modules/core/bowtie2/main.nf")
        self.assertIn("--un-conc-gz", bowtie)
        checkpoint = self.read("bin/manage_sra_checkpoints.py")
        for evidence in (
            "validate_pair",
            "run_manifest_sha256",
            "paired_fastq_records",
            "sha256_file",
            "completion_record_missing",
        ):
            self.assertIn(evidence, checkpoint)
        self.assertIn("--force", self.read("modules/local/sra_preprocessing/main.nf"))

    def test_both_input_modes_share_the_complete_global_scientific_workflow(self) -> None:
        local = self.read("workflows/metagenomics.nf")
        sra = self.read("workflows/sra_global.nf")
        for wrapper in (local, sra):
            self.assertIn("include { METAGENOMICS_GLOBAL }", wrapper)
            self.assertIn("METAGENOMICS_GLOBAL(", wrapper)

        global_workflow = self.read("subworkflows/local/metagenomics_global/main.nf")
        for stage in (
            "MAG_CONSTRUCTION",
            "TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS",
            "GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION",
            "MAG_ABUNDANCE_ESTIMATION",
            "GLOBAL_PROCESSING_EVALUATION",
        ):
            self.assertIn(f"include {{ {stage} }}", global_workflow)
            self.assertIn(f"{stage}(", global_workflow)
        self.assertIn(
            "MAG_ABUNDANCE_ESTIMATION(\n        MAG_CONSTRUCTION.out.final_mags,\n        ch_filtered_reads",
            global_workflow,
        )
        self.assertIn("ch_preprocessing_reports", global_workflow)
        self.assertIn("GLOBAL_PROCESSING_EVALUATION(ch_global_reports", global_workflow)

        mag_construction = self.read("subworkflows/local/mag_construction/main.nf")
        self.assertIn("MEGAHIT_BRANCH", mag_construction)
        self.assertIn("SPADES_BRANCH", mag_construction)
        self.assertIn("FINAL_MAG_CATALOG", mag_construction)

    def test_sequential_launcher_orders_complete_sample_lifecycle(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        sample_loop = launcher.index('for sample_id in "${pending_samples[@]}"')
        positions = [
            launcher.rfind("checkpoints_initial", 0, sample_loop),
            sample_loop,
            launcher.index('validate_sample_checkpoint "${sample_id}"', sample_loop),
            launcher.index('safe_remove_sample_work "${sample_work}"', sample_loop),
            launcher.index("checkpoints_final", sample_loop),
            launcher.index("sra-global", launcher.index("checkpoints_final", sample_loop)),
        ]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("validate-sample", launcher)
        self.assertIn("frozen_platforms", launcher)
        monitor_start = launcher.index(
            'start_storage_monitor "${work_dir}" "${sra_checkpoint_dir}"'
        )
        frozen_validation = launcher.index("validate_frozen_sra_state", monitor_start)
        global_fast_path = launcher.index(
            'if [[ -f "${state_dir}/sra_global_success.json"'
        )
        self.assertLess(frozen_validation, global_fast_path)

        checkpoint_module = self.read("modules/local/sra_preprocessing/main.nf")
        self.assertEqual(checkpoint_module.count("cache false"), 2)
        for process_name in ("PERSIST_SRA_CHECKPOINT", "CHECK_SRA_CHECKPOINTS"):
            process_start = checkpoint_module.index(f"process {process_name}")
            next_process = checkpoint_module.find("\nprocess ", process_start + 1)
            process_body = checkpoint_module[
                process_start : next_process if next_process >= 0 else None
            ]
            self.assertIn("cache false", process_body)

    def test_sra_automatically_enables_disk_efficient_profile(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        gate = launcher.index(
            'if [[ "${storage_constrained}" == true || -n "${sra_project}"'
        )
        profile = launcher.index("profiles+=(disk_efficient)", gate)
        end = launcher.index("fi", profile)
        self.assertLess(gate, profile)
        self.assertLess(profile, end)
        self.assertIn("SRA mode always does this", launcher)
        self.assertIn("executor.queueSize = 1", self.read("conf/disk_efficient.config"))
        self.assertIn(
            "local Conda/Apptainer/Singularity GPU mode requires CUDA_VISIBLE_DEVICES",
            launcher,
        )
        self.assertIn('--gpuTelemetryInterval "${gpu_telemetry_interval}"', launcher)
        self.assertIn("--gpu-telemetry-interval", self.read("README.md"))

    def test_launcher_requires_supported_host_python(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        self.assertIn('MINIMUM_PYTHON_VERSION="3.10.0"', launcher)
        self.assertIn("check_python", launcher)
        self.assertIn(
            'version_at_least "${detected_version}" "${MINIMUM_PYTHON_VERSION}"',
            launcher,
        )
        self.assertIn("Python 3.10 or newer", self.read("README.md"))

    def test_launcher_locks_results_and_shared_checkpoints_fail_closed(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        self.assertIn("acquire_run_lock()", launcher)
        self.assertIn("release_run_locks()", launcher)
        self.assertIn(
            "another launcher may be active. This is fail-closed",
            launcher,
        )
        self.assertIn("inspect owner.tsv", launcher)
        self.assertNotIn("stale lock", launcher.lower())
        self.assertIn(
            'acquire_run_lock "${resource_root}/.metagenomics_run.lock" results',
            launcher,
        )
        self.assertIn(
            'acquire_run_lock "${sra_checkpoint_dir}.metagenomics_run.lock" checkpoint',
            launcher,
        )
        self.assertIn("rmdir -- \"${lock_dir}\"", launcher)
        self.assertNotIn('rm -rf -- "${lock_dir}"', launcher)

        production_gate = launcher.index('if [[ "${dry_run}" == false ]]')
        results_lock = launcher.index(
            'acquire_run_lock "${resource_root}/.metagenomics_run.lock" results',
            production_gate,
        )
        telemetry_initialization = launcher.index(
            "initialize_project_telemetry", results_lock
        )
        self.assertLess(results_lock, telemetry_initialization)

        checkpoint_validation = launcher.index(
            'require_interpolation_safe_path "--outdir" "${outdir}"'
        )
        checkpoint_lock = launcher.index(
            'acquire_run_lock "${sra_checkpoint_dir}.metagenomics_run.lock" checkpoint',
            checkpoint_validation,
        )
        checkpoint_creation = launcher.index(
            'mkdir -p -- "${sra_checkpoint_dir}"', checkpoint_lock
        )
        storage_monitor = launcher.index(
            'start_storage_monitor "${work_dir}" "${sra_checkpoint_dir}"',
            checkpoint_lock,
        )
        self.assertLess(checkpoint_validation, checkpoint_lock)
        self.assertLess(checkpoint_lock, checkpoint_creation)
        self.assertLess(checkpoint_lock, storage_monitor)

    def test_cleanup_is_bound_to_hashed_global_outputs_and_forced_storage_samples(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        for cleanup_function in (
            "safe_remove_sample_work()",
            "safe_remove_global_work()",
            "cleanup_sra_checkpoints()",
        ):
            start = launcher.index(cleanup_function)
            body_end = launcher.index("\n}", start)
            self.assertIn("force_storage_sample", launcher[start:body_end])
        self.assertIn(
            'cleanup --checkpoint-manifest "${checkpoint_manifest}"',
            launcher,
        )
        self.assertIn('--success-marker "${success_marker}"', launcher)
        self.assertIn('"status": "in_progress"', self.read("bin/manage_sra_checkpoints.py"))
        self.assertIn("planned_files", self.read("bin/manage_sra_checkpoints.py"))

        finalizer = self.read("modules/local/sra_preprocessing/main.nf")
        for artifact in (
            "checkpoint_manifest",
            "multiqc_report",
            "software_versions",
            "mag_abundance",
        ):
            self.assertIn(f"'{artifact}': describe(", finalizer)
        self.assertIn("'sha256': digest.hexdigest()", finalizer)

        checkpoint_manager = self.read("bin/manage_sra_checkpoints.py")
        for refusal in (
            "checkpoint manifest path differs",
            "checkpoint manifest size differs",
            "checkpoint manifest SHA-256 differs",
            "durable global output is missing",
            "durable global output size differs",
            "durable global output SHA-256 differs",
            "cleanup refused size-mismatched checkpoint",
            "cleanup refused hash-mismatched checkpoint",
        ):
            self.assertIn(refusal, checkpoint_manager)
        cleanup_start = checkpoint_manager.index("def cleanup(")
        self.assertLess(
            checkpoint_manager.index("if args.keep:", cleanup_start),
            checkpoint_manager.index("path.unlink()", cleanup_start),
        )

    def test_cleanup_requires_an_exact_complete_scientific_output_seal(self) -> None:
        checkpoint_manager = self.read("bin/manage_sra_checkpoints.py")
        self.assertIn("SCIENTIFIC_RESULT_ROOTS", checkpoint_manager)
        for root in (
            "01_quality_control_and_filtering",
            "02_mag_construction",
            "03_taxonomic_classification_and_phylogenomics",
            "04_gene_prediction_and_functional_annotation",
            "05_mag_abundance_estimation",
            "06_global_processing_evaluation",
        ):
            self.assertIn(f'"{root}"', checkpoint_manager)

        required_start = checkpoint_manager.index("REQUIRED_SCIENTIFIC_ARTIFACTS")
        required_end = checkpoint_manager.index("BASELINE_OUTPUT_PATHS", required_start)
        required_contract = checkpoint_manager[required_start:required_end]
        for label in (
            "final_mags",
            "final_catalog_provenance",
            "final_catalog_quality",
            "final_checkm2",
            "final_gunc",
            "species_representatives",
            "gtdbtk_summaries",
            "phylogenomic_tree",
            "functional_annotations",
            "mag_abundance",
            "multiqc_report",
            "software_versions",
        ):
            self.assertIn(f'"{label}"', required_contract)

        self.assertIn('commands.add_parser(\n        "seal-global"', checkpoint_manager)
        for field in (
            "schema_version",
            "results_root",
            "file_count",
            "total_bytes",
            "files",
            "required_artifacts",
            "sealed_at_utc",
            "relative_path",
            "bytes",
            "sha256",
        ):
            self.assertIn(f'"{field}"', checkpoint_manager)
        self.assertIn("functional annotation MAG identifiers do not exactly match", checkpoint_manager)
        self.assertIn("scientific result inventory refuses symbolic link", checkpoint_manager)
        self.assertIn("missing, tampered, or extra artifact", checkpoint_manager)

        cleanup_start = checkpoint_manager.index("def cleanup(")
        seal_validation = checkpoint_manager.index(
            "validate_scientific_inventory(success)", cleanup_start
        )
        first_unlink = checkpoint_manager.index("path.unlink()", cleanup_start)
        self.assertLess(seal_validation, first_unlink)

        launcher = self.read("metagenomics_pipeline.sh")
        seal_helper_start = launcher.index("seal_sra_global_outputs()")
        seal_helper_end = launcher.index("\n}", seal_helper_start)
        seal_helper = launcher[seal_helper_start:seal_helper_end]
        for token in (
            "seal-global",
            "--success-marker",
            "--checkpoint-manifest",
            "--results-dir",
            "--project-accession",
        ):
            self.assertIn(token, seal_helper)

        fast_path = launcher.index(
            'if [[ -f "${state_dir}/sra_global_success.json" ]]'
        )
        fast_seal = launcher.index("seal_sra_global_outputs", fast_path)
        fast_cleanup = launcher.index("cleanup_sra_checkpoints", fast_seal)
        fast_work_cleanup = launcher.index("safe_remove_global_work", fast_cleanup)
        self.assertLess(fast_seal, fast_cleanup)
        self.assertLess(fast_cleanup, fast_work_cleanup)

        global_run = launcher.rindex("run_nextflow global sra-global")
        marker_gate = launcher.index(
            '[[ -s "${success_marker}" ]]', global_run
        )
        normal_seal = launcher.index("seal_sra_global_outputs", marker_gate)
        normal_cleanup = launcher.index("cleanup_sra_checkpoints", normal_seal)
        normal_work_cleanup = launcher.index("safe_remove_global_work", normal_cleanup)
        self.assertLess(global_run, marker_gate)
        self.assertLess(marker_gate, normal_seal)
        self.assertLess(normal_seal, normal_cleanup)
        self.assertLess(normal_cleanup, normal_work_cleanup)

    def test_gpu_requests_are_limited_to_verified_candidates(self) -> None:
        gpu_config = self.read("conf/gpu.config")
        selector = re.search(r"withName:\s*'([^']+)'\s*\{\s*accelerator", gpu_config)
        self.assertIsNotNone(selector)
        self.assertEqual(set(selector.group(1).split("|")), {"COMEBIN", "SEMIBIN2", "VAMB"})
        self.assertNotIn("MEGAHIT|", selector.group(1))
        self.assertIn("params.enableGpu = true", gpu_config)
        self.assertIn("includeConfig 'conf/gpu.config'", self.read("nextflow.config"))
        self.assertIn("--engine ${compute_engine}", self.read("modules/core/semibin2/main.nf"))
        self.assertIn("--cuda", self.read("modules/core/vamb/main.nf"))
        self.assertIn("CUDA_VISIBLE_DEVICES=''", self.read("modules/core/comebin/main.nf"))

    def test_gpu_metrics_are_attempt_and_session_scoped_at_source_and_join(self) -> None:
        for module in ("comebin", "semibin2", "vamb"):
            source = self.read(f"modules/core/{module}/main.nf")
            self.assertIn("__SESSION_ID__.__ATTEMPT__", source)
            self.assertIn("session_id\\tattempt\\tgpu_index", source)
            self.assertIn("workflow.sessionId.toString()", source)
            self.assertIn("task.attempt.toString()", source)
            self.assertIn("torch.cuda.device_count() != 1", source)
            self.assertIn("path('*.gpu_metrics.tsv'), optional: true", source)

        module_config = self.read("conf/modules.config")
        self.assertEqual(module_config.count("pattern: '*.gpu_metrics.tsv'"), 3)

        summarizer = self.read("bin/summarize_resources.py")
        self.assertIn("(session_id, process, sample_id, attempt)", summarizer)
        self.assertIn("metrics.get((session_id, process, tag, int(attempt)))", summarizer)
        trace = self.read("conf/base.config")
        self.assertIn("accelerator,accelerator_type", trace)
        self.assertIn("attempt", trace)

    def test_staged_resume_uses_the_last_uuid_for_the_same_invocation_key(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        self.assertIn("latest_resume_session()", launcher)
        self.assertIn('latest_resume_session "${invocation_id}"', launcher)
        self.assertIn('command+=(-resume "${resume_session}")', launcher)
        self.assertIn('record_resume_session "${invocation_id}"', launcher)
        self.assertIn("resume_key\\tsession_id\\trecorded_at_utc", launcher)
        self.assertNotIn("command+=(-resume)\n", launcher)

    def test_sra_publication_and_internal_launcher_controls_are_cleanup_safe(self) -> None:
        launcher = self.read("metagenomics_pipeline.sh")
        self.assertIn("--publish_dir_mode copy", launcher)
        self.assertIn("SRA mode requires --publish_dir_mode copy", launcher)
        self.assertIn("--env HOME=${sra_container_home}", launcher)
        self.assertIn('"${sra_container_home}" "${preprocessing_work_root}"', launcher)
        self.assertIn(
            "runOptions = '-u $(id -u):$(id -g)'",
            self.read("conf/docker.config"),
        )
        for reserved in ("-work-dir", "-profile", "--executionStage", "--telemetryDir"):
            self.assertIn(reserved, launcher)
        self.assertIn("the '--' passthrough delimiter is not supported", launcher)

        checkpoint_manager = self.read("bin/manage_sra_checkpoints.py")
        self.assertIn("durable regular copy, not a symlink", checkpoint_manager)
        self.assertIn("checkpoint report SHA-256 differs", checkpoint_manager)
        self.assertIn('"sha256": sha256_file(path)', checkpoint_manager)

    def test_heavy_native_temporaries_are_removed_only_after_outputs(self) -> None:
        megahit = self.read("modules/core/megahit/main.nf")
        self.assertLess(
            megahit.index('test -s "${prefix}.contigs.fa"'),
            megahit.index('rm -rf -- "${prefix}.megahit"'),
        )
        spades = self.read("modules/core/spades/main.nf")
        self.assertLess(
            spades.index('test -s "${prefix}.contigs.fa"'),
            spades.index('rm -rf -- "${prefix}.spades"'),
        )
        for module, suffix in (
            ("comebin", "comebin"),
            ("semibin2", "semibin2"),
            ("vamb", "vamb"),
        ):
            source = self.read(f"modules/core/{module}/main.nf")
            self.assertIn("keep_native_outputs", source)
            self.assertIn(f"rm -rf -- '${{prefix}}.{suffix}'", source)

        for module, suffix, required_output in (
            ("metaquast", "metaquast", 'test -s "${prefix}.metaquast.report.tsv"'),
            ("checkm2", "checkm2", 'test -s "${prefix}.checkm2.quality_report.tsv"'),
            ("gunc", "gunc", 'test -s "${prefix}.gunc.summary.tsv"'),
            ("drep", "drep", 'test -s "${prefix}.clusters.csv"'),
        ):
            source = self.read(f"modules/core/{module}/main.nf")
            self.assertIn("optional: true", source)
            self.assertLess(
                source.index(required_output),
                source.index(f"rm -rf -- '${{prefix}}.{suffix}'"),
            )

        gunc = self.read("modules/core/gunc/main.nf")
        self.assertLess(gunc.index('mkdir -p input_bins "${prefix}.gunc"'), gunc.index("gunc run"))
        self.assertIn('--temp_dir "${prefix}.gunc.tmp"', gunc)

        gtdbtk = self.read("modules/core/gtdbtk/main.nf")
        self.assertIn('--scratch_dir "\\$PWD/gtdbtk_pplacer_scratch"', gtdbtk)
        self.assertIn('--tmpdir "\\$PWD/gtdbtk_tmp"', gtdbtk)
        self.assertLess(
            gtdbtk.index('test -s "${prefix}.gtdbtk.ar53.summary.tsv"'),
            gtdbtk.index("rm -rf -- genomes gtdbtk_pplacer_scratch gtdbtk_tmp"),
        )


if __name__ == "__main__":
    unittest.main()
