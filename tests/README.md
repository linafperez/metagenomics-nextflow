# Testing and benchmark

Local tests use deterministic synthetic reads and a synthetic host reference
created under ignored `tests/generated_data/` and `tests/generated_reference/`
directories. They never download the study dataset or a production human
reference.

Generate the fixtures directly when needed:

```bash
python3 tests/scripts/generate_synthetic_data.py
```

The production `--stub` mode traverses the complete scientific graph with
structured module stubs. The scoped `--test-local` mode selects
`tests/workflows/synthetic_real.nf`, builds a genuine tiny Bowtie2 index from the
generated FASTA, and runs real FastQC, fastp, host removal, MEGAHIT, SPAdes, and
MetaQUAST processes. It stops before database-heavy and licensed stages.

```bash
./metagenomics_pipeline.sh --local --conda --stub
./metagenomics_pipeline.sh --local --conda --test-local
```

The launcher regenerates missing fixtures automatically. See
`tests/data/README.md` and `tests/reference/README.md` for their locations and
contracts.

## Assembler and binner benchmark

`tests/workflows/benchmark.nf` is isolated from the production entrypoint. It
evaluates MEGAHIT, SPAdes, or both with COMEBin, MetaBAT2, SemiBin2, Vamb, or
all four binners. Individual-binner variants pass that binner's bins directly
to MAG refinement. Only the `all` strategy runs DAS Tool.

Generate ignored local fixtures and validate the complete stub matrix:

```bash
python3 tests/scripts/generate_synthetic_data.py
python3 tests/scripts/run_variants.py --local --stub
```

An optional runtime can be included for backend validation:

```bash
python3 tests/scripts/run_variants.py --local --docker --stub
python3 tests/scripts/run_variants.py --local --conda --stub
python3 tests/scripts/run_variants.py --local --apptainer --stub
```

For a later real HPC run, provide filtered non-host paired reads and prepared
CheckM2 and GUNC databases:

```bash
python3 tests/scripts/run_variants.py \
    --hpc --apptainer --run \
    --input /shared/project/filtered_reads.csv \
    --checkm2-db /shared/databases/checkm2/uniref100.KO.1.dmnd \
    --gunc-db /shared/databases/gunc/gunc_db.dmnd \
    --slurm-account project_account --slurm-queue compute \
    --jobs 3 --resume
```

Each variant is written below `tests/results/<assembler>_<binner>/`, with an
independent work directory below `tests/work/variants/`. The runner collects
Nextflow trace, report, timeline, and DAG files. It then invokes
`tests/scripts/summarize_variants.py`, which creates comparison TSV, ranking
TSV, Markdown, and JSON outputs under `tests/results/comparison/`.

Ranking is lexicographic: high-quality MAG count, median completeness, median
contamination, GUNC failures, final non-redundant MAG count, then wall-clock
runtime. Missing measurements are retained as `NA`/JSON `null` and sort after
measured values. No weighted score is used.
