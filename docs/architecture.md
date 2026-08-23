# Pipeline architecture

This document describes the implemented Nextflow DSL2 organization, channel
contracts, dependency graph, and configuration boundaries. Scientific usage and
copy-pasteable execution commands are in the repository [README](../README.md).

## Design principles

- The production entrypoint always executes the complete scientific strategy.
- Scientific independence is preserved in the DAG; executor configuration
  controls physical concurrency.
- A bioinformatics tool has one maintained module implementation even when it
  is used in multiple contexts.
- Context is expressed through aliases, metadata, inputs, `task.ext`, and
  `withName` configuration rather than copied modules.
- Assembler and MAG provenance travel in metadata and explicit tables, not
  inferred result-directory names.
- Published outputs are for users. Downstream tasks consume upstream channels
  and never read files back from `results/`.
- Large databases, licensed files, generated test biology, caches, work files,
  and results remain outside version control.
- nf-core, CLL, and legacy Bash code are references; project modules are
  maintained under `modules/core/` and no nf-core module tree is vendored.

## Repository layout

```text
main.nf                         thin production entrypoint
workflows/metagenomics.nf       complete top-level workflow
subworkflows/local/             scientific and reusable orchestration
modules/core/                   one reusable wrapper per bioinformatics tool
modules/local/                  pipeline-specific transformation processes
bin/                            transparent Python and Bash helpers
conf/                           configuration layers and orthogonal profiles
containers/                     project container build definitions
assets/                         samplesheet, MultiQC, and PhyloPhlAn templates
tests/workflows/                isolated validation and benchmark entrypoints
tests/subworkflows/             benchmark-only orchestration
tests/scripts/                  fixture, matrix, metric, and unit-test tools
docs/                           architecture documentation
```

`main.nf` only imports and invokes `METAGENOMICS`. This keeps the production
entrypoint distinct from every test and benchmark entrypoint.

## Top-level dataflow

`workflows/metagenomics.nf` performs these orchestration duties:

1. validates required production input, database, and license parameters plus
   numeric ranges and ANI relationships;
2. stages and validates the samplesheet through `CHECK_SAMPLESHEET`;
3. creates paired-read tuples and a six-file Bowtie2 index channel;
4. invokes the six scientific stages;
5. mixes genuine report inputs for global MultiQC evaluation;
6. collects module version records into `software_versions.tsv`;
7. emits the principal final deliverables.

`nextflow_schema.json` is the machine-readable production parameter contract. It
defines types, enumerations, defaults, and numeric ranges. Runtime validation in
`METAGENOMICS` repeats the scientific range checks and enforces
`derep_ani > species_ani`, so direct Nextflow execution does not depend on an
editor or external schema client for critical validation.

The complete hierarchy is:

```text
METAGENOMICS
|-- QUALITY_CONTROL_AND_FILTERING
|-- MAG_CONSTRUCTION
|   |-- MEGAHIT_BRANCH
|   |   |-- MEGAHIT_ASSEMBLY
|   |   |-- MEGAHIT_BINNING
|   |   `-- MEGAHIT_MAG_REFINEMENT
|   |-- SPADES_BRANCH
|   |   |-- SPADES_ASSEMBLY
|   |   |-- SPADES_BINNING
|   |   `-- SPADES_MAG_REFINEMENT
|   `-- FINAL_MAG_CATALOG
|-- TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS
|-- GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION
|-- MAG_ABUNDANCE_ESTIMATION
`-- GLOBAL_PROCESSING_EVALUATION
```

The last four consumers start only after `FINAL_MAG_CATALOG` emits final MAGs,
but are otherwise independent. MultiQC additionally waits for the report
channels it collects.

## Data contracts and metadata

### Read-level contract

After samplesheet normalization, paired reads use:

```text
tuple(meta, reads)

meta  = [id: sample_id, single_end: false]
reads = [fastq_1, fastq_2]
```

Read-level metadata contains only sample identity and pairing state. Assembly,
binner, and MAG fields are added only when they become relevant.

The samplesheet helper validates the exact header, unique safe sample IDs,
paired files, supported FASTQ extensions, distinct mates, and file existence.
It resolves relative paths against the samplesheet directory and publishes the
normalized CSV under `pipeline_info/`.

### Assembly contract

Each assembler subworkflow collects all filtered sample tuples, sorts them by
sample ID for deterministic input order, and emits:

```text
tuple(meta, contigs)

meta = [
    id: '<assembler>_coassembly',
    assembler: 'megahit' | 'spades',
    branch: 'megahit' | 'spades',
    sample_ids: [...]
]
```

The assembler metadata persists through contig filtering, coverage, binning,
DAS Tool, and branch refinement.

### Bin and catalog contracts

Binner and evaluation modules conventionally emit:

```text
tuple(meta, paths)
```

where `paths` is one or more bin/MAG FASTA files, a structured table, a native
result directory, or a log. When records from independent processes must be
rejoined, the subworkflow maps the tuple to an explicit `meta.id` key, joins on
that key, and then removes the temporary key.

The final merge helper assigns deterministic collision-free identifiers:

```text
MAG_MEGAHIT_0001
MAG_SPADES_0001
```

It produces provenance and quality tables containing assembler, original MAG
identifier, original file, completeness, and contamination. The finalization
helper copies only selected dRep representatives and subsets both metadata
tables to exactly the FASTA files present. It fails if FASTA, provenance, and
quality counts differ.

### Version contract

Every practical module emits a three-value record:

```text
tuple(process_name, tool_name, version)
```

The top-level workflow mixes these channels. `COLLECT_VERSIONS` normalizes,
deduplicates, sorts, and writes `software_versions.tsv`.

## Stage implementation

### Quality control and filtering

`QUALITY_CONTROL_AND_FILTERING` imports the single FastQC module twice:

```nextflow
include { FASTQC as FASTQC_RAW }   from '../../../modules/core/fastqc/main'
include { FASTQC as FASTQC_CLEAN } from '../../../modules/core/fastqc/main'
```

Fully qualified `withName` selectors distinguish prefixes and publication
paths. Keyed joins deliberately make fastp wait for raw FastQC and Bowtie2 wait
for clean FastQC while retaining the original/cleaned read tuples. This models
the requested stage order without treating FastQC as a read-transforming or
pass/fail process.

The host index is a value channel containing the six `.bt2` or `.bt2l` files
resolved from the configured prefix. Bowtie2 receives the prefix basename
separately and emits paired non-host reads plus its stderr summary for MultiQC.

### Assembly branches

`MEGAHIT_BRANCH` and `SPADES_BRANCH` are separate invocations from the same
filtered-read source. There is no edge between them.

MEGAHIT and SPAdes have distinct tool modules because they are different tools.
Both assembly subworkflows reuse the single MetaQUAST implementation through
the aliases `METAQUAST_MEGAHIT` and `METAQUAST_SPADES`. Metadata controls the
assembler-specific publication directory.

The MEGAHIT process uses all forward reads as one comma-separated input and all
reverse reads as the paired counterpart, `meta-large`, a 1,000 bp minimum, and
90% of available process memory. SPAdes builds a deterministic multi-library
dataset and invokes `spades.py --meta` with process CPUs and memory.

### Shared binning core

Both assembler-specific binning wrappers call the same `BINNING` subworkflow:

1. `FILTER_CONTIGS` applies `params.min_contig_length` before reconstruction.
2. All filtered read pairs are collected in deterministic sample order.
3. `COVERM_CONTIG` maps reads once and emits cached per-sample BAMs, MetaBAT2
   depth, and the dense mean-depth matrix required by Vamb.
4. COMEBin and SemiBin2 consume contigs plus BAMs; MetaBAT2 consumes contigs
   plus MetaBAT depth; Vamb consumes contigs plus its abundance matrix.
5. The four binners are invoked without edges between one another.
6. DAS Tool waits for all four contig-to-bin mappings and the filtered contigs.

The binning core emits every individual bin set, mapping table, coverage
representation, DAS Tool bin set, native summary, log, and version channel.
Assembler wrappers retain obvious MEGAHIT/SPAdes scientific branches while
avoiding duplicated tool code.

### Shared MAG-refinement core

Both branch refinement wrappers call `MAG_REFINEMENT`:

```text
DAS Tool MAGs
  |-- CheckM2 raw
  `-- GUNC raw
        -> CheckM2-based strict HQ selection (>90 completeness, <5 contamination)
        -> dRep at 99% ANI
        -> dRep at 95% ANI for species groups
        |-- CheckM2 clean on 99% representatives
        `-- GUNC clean on 99% representatives
```

GUNC can run alongside raw CheckM2 because its result does not gate the
scientific HQ definition. It remains a preserved chimerism/taxonomic-consistency
assessment. The selection helper writes the selected FASTA set, an audit TSV,
and dRep-compatible `genomeInfo.csv`.

The reusable dRep process receives ANI, coverage, and stage values. Both 99%
and 95% operations use:

- primary ANI (`-pa`) 0.90;
- secondary ANI (`-sa`) 0.99 or 0.95;
- minimum comparable coverage (`-nc`) 0.30;
- larger-coverage comparison mode;
- fastANI secondary comparisons;
- multiround primary clustering.

`-comp 0 -con 100` prevents dRep from silently replacing the already applied
strict external quality filter; the validated `genomeInfo.csv` supplies the
quality values used for representative scoring.

### Final combined catalog

`FINAL_MAG_CATALOG` is owned by the higher-level MAG construction stage, not an
assembler branch. It joins each branch's 99% representative set with that
branch's clean CheckM2 table, creates unique IDs, runs cross-assembler dRep at
99%, and separately creates 95% species groups from those representatives.

`FINALIZE_MAG_CATALOG` defines the production catalog from the final 99%
representatives and enforces exact provenance/quality-table membership. Final
CheckM2 and GUNC runs assess this emitted catalog. The final 95% representatives
and both cluster tables are exposed separately for species-level interpretation.

### Taxonomy and phylogenomics

GTDB-Tk and PhyloPhlAn independently consume the final MAG tuple. GTDB-Tk uses
`classify_wf`, a staged batch file, GTDB release 226, and separate bacterial and
archaeal summary outputs. PhyloPhlAn stages collision-checked genome IDs,
conserved marker resources, and the versioned configuration in
`assets/phylophlan_iqtree.cfg`.

The PhyloPhlAn process checks that its backend is IQ-TREE 3.0.1 before running.
The Conda environment contains the exact PhyloPhlAn/IQ-TREE pair. Container
profiles use `params.phylophlan_container` to select the repository-built image
containing both versions. IQ-TREE is a backend inside this process, not an
independent scientific subworkflow.

### Gene prediction and functional annotation

The functional subworkflow expands the final catalog to one tuple per MAG and
derives collision-checked MAG identifiers. GeneMarkS-2 emits proteins, coding
nucleotides, and GFF features with MAG-aware protein identifiers. The process
stages the externally licensed installation and key in its isolated task and
does not publish or copy the key to results.

eggNOG-mapper and InterProScan consume the same predicted-protein channel and
can run concurrently. The integration join is keyed by `mag_id` and combines:

```text
proteins + GeneMark GFF + eggNOG annotations + InterProScan TSV
```

The helper produces one row per predicted protein and a JSON summary. It
retains identifiers, coordinates, sequence length, preferred name,
description, ortholog and COG fields, GO, EC, KEGG, CAZy, PFAM, InterPro
accessions/descriptions, member databases, pathways, and contributing sources.

The eggNOG database contract intentionally validates mapper database 5.0.2,
including its annotation, taxonomy, and DIAMOND files. The preparer uses the
current official download host because the URL embedded in mapper 2.1.13 is
obsolete. eggNOG-mapper 2.1.13 is not compatible with an eggNOG 6.0 mapper
database.

### MAG abundance

The abundance subworkflow combines final MAGs with all filtered paired reads.
`COVERM_GENOME` uses `minimap2-sr`, properly paired reads only, 95% minimum read
identity, and 75% minimum aligned-read percentage. CoverM 0.7.0 reports covered
bases and genome length rather than a native `covered_fraction` method, so the
wrapper derives covered fraction as covered bases divided by length.

The native wide table is retained. `NORMALIZE_ABUNDANCE` validates every sample
metric group and creates a long table with:

```text
sample, mag_id, relative_abundance_percent, mean_coverage,
covered_fraction, genome_length
```

### Global processing evaluation

The top-level workflow mixes report and log outputs only after their producers
finish. `GLOBAL_PROCESSING_EVALUATION` recursively extracts path values,
deduplicates them, collects them into one staged input, and invokes the single
MultiQC 1.35 process with `assets/multiqc_config.yml`.

The configuration deliberately enables supported parsers for FastQC (separate
raw and clean sections), fastp, Bowtie2, MEGAHIT, QUAST/MetaQUAST, CheckM2, and
GTDB-Tk. It does not pretend that CoverM, the four binners, DAS Tool, GUNC,
dRep, PhyloPhlAn, GeneMarkS-2, eggNOG-mapper, or InterProScan are parsed when no
reliable native parser is configured. Those native outputs remain published.

## Reusable modules and stubs

Every practical process has a `stub:` block. Stubs create minimal structured
FASTA, FASTQ, TSV, CSV, JSON, GFF, Newick, HTML, or directory outputs as required
by downstream consumers. They do more than create empty files when a later
process parses the content.

The `stub` profile sets conservative resources and replaces tool runtime
definitions with a Python base runtime. `-stub-run` selects each process stub;
`--stub_run true` also switches production resource selectors to small values.
This dual setting is intentional.

Stubs validate:

- channel shapes and keyed joins;
- branch separation and convergence;
- metadata propagation;
- filename collision prevention;
- strict quality-table selection;
- dRep representative/group contracts;
- final catalog metadata subsetting;
- downstream annotation and abundance joins;
- MultiQC input collection and version aggregation.

They do not validate bioinformatics executables, packaging, databases, a
GeneMark license, scientific performance, or SLURM scheduling.

## Configuration composition

`nextflow.config` includes three global layers:

| File | Responsibility |
| --- | --- |
| `conf/base.config` | shell safety, retry behavior, report, timeline, trace, and DAG |
| `conf/resources.config` | process CPU, memory, and time requests |
| `conf/modules.config` | process arguments, prefixes, publication paths, and optional intermediates |

Profiles add independent dimensions:

| Profile | Responsibility |
| --- | --- |
| `local` | local executor, resource ceilings, local queue size |
| `hpc` | SLURM executor, account/queue/QoS/site options, queue size |
| `docker` | Docker enabled; other runtimes disabled |
| `conda` | Conda environments with required Mamba solving; containers disabled |
| `apptainer` | Apptainer with automatic mounts; other runtimes disabled |
| `singularity` | Singularity compatibility with automatic mounts |
| `test` | ignored generated fixture paths and conservative local ceilings |
| `stub` | structured stub resources and generic stub runtime |

The expected production composition is `-profile <environment>,<runtime>`, for
example `local,docker` or `hpc,apptainer`. Test and stub profiles are additional
layers, not substitutes for the environment/runtime distinction.

The Conda profile requires both Conda and the `mamba` executable. Its cache uses
`params.conda_cache_dir` when set and otherwise falls back to the ignored
repository directory `.conda/`.

Docker, Apptainer, and Singularity production runs require a combined image with
PhyloPhlAn 3.1.1 and IQ-TREE 3.0.1 for the `PHYLOPHLAN` process. The stock
PhyloPhlAn image does not satisfy that complete runtime contract. The repository
build definition is `containers/phylophlan-iqtree/Dockerfile`; users build or
publish it and pass the OCI reference through `params.phylophlan_container`.
Apptainer and Singularity consume the same published image. This prerequisite
does not affect structured stubs or the scoped local real-tool workflow, which
do not run production phylogenomics.

Local real-tool validation uses the separate
`tests/workflows/synthetic_real.nf` entrypoint and
`tests/config/synthetic_real.config`. It builds a genuine tiny Bowtie2 index and
runs FastQC, fastp, host removal, both assemblers, and both MetaQUAST evaluations
with generated reads. It intentionally omits database-heavy and licensed stages
rather than substituting fake databases in a real-tool run.

Resource selectors preserve the legacy scheduling intent. Examples include 32
CPUs/500 GB for MEGAHIT, 32 CPUs/1,800 GB for SPAdes, 24 CPUs/250 GB for DAS
Tool, 24 CPUs/100 GB for dRep and GTDB-Tk, 32 CPUs/120 GB for PhyloPhlAn, and
16 CPUs/128 GB for InterProScan. Local/HPC `resourceLimits` are configurable
ceilings; they do not alter scientific dependencies.

## Publication boundary

`conf/modules.config` maps outputs into the numbered scientific hierarchy. The
default publication mode is `copy`. Large optional products are controlled by:

- `params.save_clean_reads`;
- `params.save_host_removed_reads`;
- `params.save_bam`;
- `params.save_intermediates`.

All processes still exchange required intermediates in Nextflow work
directories even when an optional user-facing publication is disabled.

## Database preparation boundary

`bin/prepare_databases.sh` is intentionally outside the scientific DAG. It can
inspect, selectively prepare, resume, checksum, or explicitly rebuild external
resources without coupling normal runs to network access. It validates that
the database root is not `/`, `/tmp`, this repository, or another Git worktree,
uses staging and a lock, and writes a manifest plus a Nextflow configuration.

GeneMarkS-2 is validation-only in this mechanism. Its licensed installer and
key are never downloaded or copied. Normal and HPC test runs only consume the
prepared configuration; they do not invoke database preparation.

## Benchmark isolation

`tests/workflows/benchmark.nf` is a separate entrypoint. Its selectors exist
only in test configuration:

```text
benchmark_assembler = megahit | spades | both
benchmark_binner    = comebin | metabat2 | semibin2 | vamb | all
```

`tests/subworkflows/benchmark_binning` reuses production contig filtering,
CoverM, and binner modules. Single-binner strategies emit that bin set directly
to the shared production refinement subworkflow. Only `all` joins all four
mappings into DAS Tool. `both` invokes the normal final combined-catalog
subworkflow after independent refinement branches.

`tests/scripts/run_variants.py` creates isolated output/work directories and
launches one or all 15 combinations, optionally concurrently. The benchmark
configuration writes an execution trace, report, timeline, and DAG per variant.

`tests/scripts/summarize_variants.py` reads native structured outputs and never
uses MultiQC as its sole metric source. It records missing values explicitly and
uses the documented lexicographic ranking:

1. strict HQ MAG count;
2. median completeness;
3. median contamination;
4. GUNC failure count;
5. final non-redundant MAG count;
6. wall-clock runtime tie-breaker.

No weighted or hidden composite score is calculated.

## External validation boundary

Static checks, helper tests, and structured stubs can run without production
databases. Real backend validation is conditional on the environment. Docker
requires a running daemon; the Conda profile requires Conda, Mamba, and
environment solving; Apptainer/Singularity require their executables and image
access; container production requires the combined PhyloPhlAn/IQ-TREE image;
GeneMarkS-2 requires a valid license; database-heavy stages require their full
versioned resources; and HPC validation requires SLURM and shared storage.

A configured profile must not be described as runtime-tested until its command
has completed in that environment. The repository keeps pending environment,
database, license, and HPC checks distinct from successful static or stub
validation.
