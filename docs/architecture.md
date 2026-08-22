# Pipeline architecture

## 1. Entry point: `main.nf`

`main.nf` is intentionally small. It imports and invokes the main scientific
workflow, `METAGENOMICS`.

## 2. Main workflow: `workflows/metagenomics.nf`

The main workflow validates public parameters, validates the samplesheet,
converts its rows to Nextflow tuples, stages the Bowtie2 index, and invokes the
currently implemented scientific subworkflow.

The current tuple shape is:

```text
[meta, reads]
```

where `meta` is `[id: sample_id, single_end: false]` and `reads` is the ordered
pair `[fastq_1, fastq_2]`.

## 3. Scientific subworkflows

The top-level scientific hierarchy is:

```text
METAGENOMICS
|-- QUALITY_CONTROL_AND_FILTERING
|-- MAG_CONSTRUCTION
|-- TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS
|-- GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION
|-- MAG_ABUNDANCE_ESTIMATION
`-- GLOBAL_PROCESSING_EVALUATION
```

Only `QUALITY_CONTROL_AND_FILTERING` is executable in Phase 1. The other files
are intentionally empty workflow skeletons with no fake outputs.

## 4. Nested `MAG_CONSTRUCTION`

`MAG_CONSTRUCTION` preserves its scientific hierarchy by importing three nested
subworkflows:

```text
MAG_CONSTRUCTION
|-- ASSEMBLY
|-- BINNING
`-- MAG_REFINEMENT
```

These workflows are not invoked in Phase 1.

## 5. Reusable modules

`modules/core/` contains reusable bioinformatics wrappers written for this
project. A tool is implemented once and reused through channels, aliases,
metadata, and configuration.

`modules/local/` contains pipeline-specific helper processes. In Phase 1 it
contains only the samplesheet validation process, backed by
`bin/check_samplesheet.py`.

nf-core pipelines and components are architectural references only. No nf-core
module is downloaded, installed, or vendored in this repository.

## 6. FastQC reuse

`QUALITY_CONTROL_AND_FILTERING` imports the same process twice:

```text
FASTQC as FASTQC_RAW
FASTQC as FASTQC_CLEAN
```

The aliases create two invocation contexts without duplicating the process
implementation. Fully qualified `withName` selectors in `conf/modules.config`
give each context its prefix and publication path.

Keyed joins use the sample ID to preserve the requested order:

```text
raw FastQC completion -> fastp
clean FastQC completion -> Bowtie2
```

FastQC remains an assessment step; it does not create a pass/fail branch or
modify the reads.

## 7. Configuration layers

- `conf/base.config`: safe process behavior and execution reports.
- `conf/resources.config`: CPU and memory requests derived from legacy SLURM
  declarations.
- `conf/modules.config`: scientific arguments, contextual prefixes, and result
  publication.
- `conf/local.config`: local executor with Docker.
- `conf/test.config`: local test paths and limited concurrency.
- `conf/slurm.config`: preliminary SLURM executor selection only.

Scientific parameters use `--`, for example `--input`. Execution profiles use
one hyphen, for example `-profile local`.
