# Metagenomics Nextflow

This project is a modular Nextflow workflow for paired-end shotgun metagenomics.
It is being developed incrementally for local WSL2/Docker execution and future
Linux/SLURM execution.

## Development status

Phase 1 implements only `QUALITY_CONTROL_AND_FILTERING`:

```text
paired-end raw FASTQ
  -> FastQC (raw)
  -> fastp
  -> FastQC (clean)
  -> Bowtie2 host removal
  -> paired non-host FASTQ
```

FastQC is implemented once and reused through the `FASTQC_RAW` and
`FASTQC_CLEAN` aliases. The remaining scientific subworkflows are architecture
skeletons and are not executed.

## Scientific hierarchy

```text
METAGENOMICS
|-- QUALITY_CONTROL_AND_FILTERING
|-- MAG_CONSTRUCTION
|   |-- ASSEMBLY
|   |-- BINNING
|   `-- MAG_REFINEMENT
|-- TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS
|-- GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION
|-- MAG_ABUNDANCE_ESTIMATION
`-- GLOBAL_PROCESSING_EVALUATION
```

## Input

The pipeline accepts a CSV samplesheet:

```csv
sample,fastq_1,fastq_2
sample_A,/data/sample_A_R1.fastq.gz,/data/sample_A_R2.fastq.gz
```

Phase 1 requires paired-end reads, one row per unique sample, and an already
built Bowtie2 index. Relative FASTQ paths are resolved from the samplesheet
directory. Supported extensions are `.fastq`, `.fastq.gz`, `.fq`, and `.fq.gz`.

## Local execution

```bash
nextflow run . \
    -profile local \
    --input /path/to/samplesheet.csv \
    --host_bowtie2_index /path/to/index/GRCh38_p14 \
    --outdir results
```

The `local` profile uses Docker. On the first real run, Nextflow will pull the
configured FastQC, fastp, and Bowtie2 images if they are not already available.
No container image was pulled while creating Phase 1.

## Test profile

The `test` profile points to user-supplied files under `tests/`:

```bash
nextflow run . -profile test
```

It will fail until the paired FASTQ files and the six Bowtie2 index files
described under `tests/` are provided.

## SLURM profile

`-profile slurm` currently selects the SLURM executor only. It is intentionally
preliminary: account, queue, QoS, time limits, container runtime, and external
resource paths must be configured for the target HPC before use.

## Data and database policy

The pipeline does not download FASTQ files, human references, Bowtie2 indexes,
or biological databases. All required resources must be provided by the user.
