# Pipeline architecture

This document describes the implemented Nextflow DSL2 organization, channel
contracts, dependency graph, and configuration boundaries. Scientific usage and
copy-pasteable execution commands are in the repository [README](../README.md).

## Design principles

- The launcher accepts exactly one input mode: a local paired-read samplesheet
  or a public SRA BioProject. Both modes converge on the same complete global
  scientific strategy.
- Scientific independence is preserved in the DAG; executor configuration
  controls physical concurrency.
- SRA cohort discovery is frozen before acquisition. A biological sample is
  considered complete only when its durable checkpoint record validates
  against that exact frozen manifest.
- A bioinformatics tool has one maintained module implementation even when it
  is used in multiple contexts.
- Context is expressed through aliases, metadata, inputs, `task.ext`, and
  `withName` configuration rather than copied modules.
- Assembler and MAG provenance travel in metadata and explicit tables, not
  inferred result-directory names.
- Published outputs are for users. Downstream tasks consume upstream channels
  and never read files back from `results/`.
- Large databases, licensed files, caches, work files, and results remain
  outside version control.
- nf-core, CLL, and legacy Bash code are references; project modules are
  maintained under `modules/core/` and no nf-core module tree is vendored.

## Repository layout

```text
metagenomics_pipeline.sh         production launcher and staged-run controller
main.nf                          internal execution-stage dispatcher
workflows/metagenomics.nf        local-samplesheet entry workflow
workflows/sra_*.nf               SRA discovery, reconciliation, sample, global
subworkflows/local/metagenomics_global/
                                 shared post-filtering scientific workflow
subworkflows/local/             scientific and reusable orchestration
modules/core/                   one reusable wrapper per bioinformatics tool
modules/local/                  pipeline-specific transformation processes
bin/                            transparent Python and Bash helpers
conf/                           configuration layers and orthogonal profiles
containers/                     project container build definitions
assets/                         samplesheet, MultiQC, and PhyloPhlAn templates
docs/                           architecture documentation
```

`main.nf` dispatches on the internal `params.executionStage` value. `auto`
accepts only a local `--input`; BioProject execution is deliberately rejected
there because it requires lifecycle control across several Nextflow sessions.
The launcher alone selects the internal stages `local`, `sra-discovery`,
`sra-checkpoints`, `sra-preprocess`, and `sra-global`.
Its host-side resolver, checkpoint, telemetry, and accounting helpers require
Python 3.10 or newer; this prerequisite is checked before any production stage.

## Top-level dataflow

The two mutually exclusive production paths are:

```text
local samplesheet
  -> METAGENOMICS
     -> QUALITY_CONTROL_AND_FILTERING
     -> METAGENOMICS_GLOBAL

SRA BioProject
  -> SRA_PROJECT_DISCOVERY                         (once; freeze cohort)
  -> SRA_CHECKPOINT_RECONCILIATION                (find pending samples)
  -> SRA_SAMPLE_PREPROCESSING, one sample at a time
     -> SRA_ACQUIRE -> QUALITY_CONTROL_AND_FILTERING
     -> PERSIST_SRA_CHECKPOINT
  -> SRA_CHECKPOINT_RECONCILIATION --require-complete
  -> SRA_GLOBAL
     -> METAGENOMICS_GLOBAL
     -> FINALIZE_SRA_GLOBAL_RUN                    (baseline success marker)
  -> seal the complete published scientific result inventory
  -> validated checkpoint-read cleanup, unless retention was requested
```

`workflows/metagenomics.nf` performs these local-input duties:

1. validates required production input, database, and license parameters plus
   numeric ranges and ANI relationships;
2. stages and validates the samplesheet through `CHECK_SAMPLESHEET`;
3. creates paired-read tuples and a six-file Bowtie2 index channel;
4. invokes preprocessing and then the shared global scientific workflow;
5. mixes genuine report inputs for global MultiQC evaluation;
6. collects module version records into `software_versions.tsv`;
7. emits the principal final deliverables.

`SRA_GLOBAL` reconstructs the same filtered-read cohort and preprocessing
report/version inputs from validated durable checkpoints before invoking
`METAGENOMICS_GLOBAL`. It does not re-run acquisition or read preprocessing.
The shared global workflow therefore preserves the same MEGAHIT and metaSPAdes
coassemblies, four-binner/DAS Tool paths, MAG refinement and catalog, taxonomy,
functional annotation, final-catalog abundance, and MultiQC strategy in both
input modes.

`nextflow_schema.json` is the machine-readable production parameter contract. It
defines types, enumerations, defaults, and numeric ranges. Runtime validation in
`METAGENOMICS_GLOBAL` repeats the scientific range checks and enforces
`derep_ani > species_ani`, so direct Nextflow execution does not depend on an
editor or external schema client for critical validation.

The shared scientific hierarchy is:

```text
METAGENOMICS
|-- QUALITY_CONTROL_AND_FILTERING
`-- METAGENOMICS_GLOBAL
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

In SRA mode, `SRA_GLOBAL` contains the same `METAGENOMICS_GLOBAL` subtree and
adds `FINALIZE_SRA_GLOBAL_RUN` after MultiQC, the version table, and the final
long abundance table are all non-empty.

The last four consumers start only after `FINAL_MAG_CATALOG` emits final MAGs,
but are otherwise independent. MultiQC additionally waits for the report
channels it collects.

## Data contracts and metadata

### Read-level contract

After local samplesheet normalization and during per-sample SRA preprocessing,
paired reads use:

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

### Frozen SRA cohort contract

`SRA_PROJECT_DISCOVERY` emits and publishes a run manifest, biological-sample
manifest, exclusion audit, JSON summary, raw RunInfo CSV, and validation
sentinel. Eligible runs must belong to the requested BioProject, be public,
paired, `WGS`, `METAGENOMIC`, and on the configured short-read platform
allowlist. The run accession must be valid; spot, paired-spot, and base counts
must be valid and positive when supplied, and paired spots must equal total
spots when both are present. Ineligible rows remain visible in the exclusion
audit rather than silently entering the cohort.

The biological identity key is the valid BioSample accession when present,
then the Experiment accession, and finally the Run accession as an explicitly
warned fallback. All eligible runs with the same key form one biological
sample. The frozen run manifest is ordered by integer `sample_order`, then
integer `run_order`; this order and its SHA-256 bind every later checkpoint.
When frozen state already exists, the launcher validates and reuses it instead
of querying NCBI again, and rejects a requested project that differs from the
frozen project.

### Durable SRA checkpoint contract

For each pending biological sample, `SRA_SAMPLE_PREPROCESSING` acquires all of
its runs in frozen order, creates one merged paired read set, and applies the
normal raw FastQC -> fastp -> clean FastQC -> Bowtie2 host-removal chain. The
Bowtie2 non-host pair is the only large read data to cross the durable
boundary; its small reports and completion record accompany it:

```text
<checkpoint>/reads/<sample_id>_host_removed_R1.fastq.gz
<checkpoint>/reads/<sample_id>_host_removed_R2.fastq.gz
<checkpoint>/reports/<sample_id>/<native small report files>
<checkpoint>/records/<sample_id>.checkpoint.json
```

The JSON record is the commit marker, not directory or FASTQ existence. It is
atomically written only after both gzip FASTQs have been copied, fully parsed,
shown to have equal non-zero paired record counts and matching read names, and
hashed. It records the project, sample order and identity provenance,
BioSample, experiment and run accessions, frozen-manifest SHA-256, read paths,
sizes and SHA-256 values, paired record count, durable report paths,
completion time, schema version, and `status=complete`.

`SRA_CHECKPOINT_RECONCILIATION` revalidates every completion record, read pair,
hash, count, report path, and frozen-manifest association. It emits:

```text
sra_checkpoint_manifest.tsv  validated completed rows in frozen sample order
sra_pending_samples.tsv      sample_order, sample_id, reason
sra_checkpoint_status.json   expected/complete/pending counts and completion
```

The checkpoint manifest carries these exact columns:

```text
schema_version, project_accession, sample_order, sample_id, identity_source,
biosample_accession, experiment_accessions, run_count, run_accessions,
run_manifest_sha256, read_1, read_1_bytes, read_1_sha256, read_2,
read_2_bytes, read_2_sha256, paired_fastq_records, reports_json,
completed_at_utc, status
```

For the global handoff, `SRA_GLOBAL` maps each row to:

```text
tuple(meta, reads)

meta = [
    id: sample_id,
    single_end: false,
    biosample_accession: biosample_accession,
    identity_source: identity_source,
    sample_order: sample_order as Integer,
    run_accessions: run_accessions split on ';'
]
reads = [read_1, read_2]
```

Non-version paths from `reports_json` form the preprocessing report channel.
Persisted `*_versions.tsv` rows are converted back to the exact
`tuple(process_name, tool_name, version)` contract and mixed with resolver and
checkpoint-manager provenance before global version collection.

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

## SRA lifecycle, failure, and resume semantics

The launcher requires checkpoint and scratch roots for SRA mode and resolves
all storage paths before execution. The repository, results, checkpoint, and
scratch roots must be safe, non-broad, and non-overlapping. Checkpoints are
durable input to the later global run; acquisition cache, fasterq temporary
space, and per-sample Nextflow work are explicitly disposable.
Docker SRA processes run as the invoking host UID:GID and use a writable home
under acquisition scratch. Consequently the non-root SRA image and the Python
checkpoint containers create files that the host-side validator and cleanup
controller can read and remove; Apptainer, Singularity, and Conda already
preserve the invoking identity.
BioProject mode forces `publish_dir_mode=copy`; cleanup validation rejects a
checkpoint manifest or bound global result that is a symlink rather than a
durable regular copy.

Before results state or telemetry is initialized, the launcher acquires an
atomic results-state lock. SRA mode normalizes and safety-checks every storage
path, then acquires a second, sibling lock for the external checkpoint root
before that root is created or touched, the storage monitor starts, or a
scientific stage runs. The sibling location leaves a previously unused
checkpoint directory empty for its immutable ownership claim. Acquisition is
ordered results first and checkpoint second, so competing launchers cannot
deadlock. If either atomic directory creation finds an existing lock, the new
launcher fails closed and releases only locks that it acquired itself. Lock
metadata identifies the process, host, start time, results root, and checkpoint
root. Locks are released through the launcher exit handler, but are never
automatically broken as "stale": a PID cannot establish liveness across HPC
hosts. An operator must inspect the metadata and prove that the owner is gone
before manually removing that exact lock directory. Dry runs acquire no lock
and mutate no lifecycle state.

Staged results and SRA paths reject quotes, backticks, dollar signs,
backslashes, and line breaks because they cross generated-script boundaries.
Container-backed SRA storage paths additionally reject whitespace, comma, and
colon because they are embedded in bind specifications. Separately measured
acquisition, cache, temporary, and work roots may not overlap, preventing
double-counted storage categories.

The launcher performs an initial reconciliation and iterates only the pending
sample IDs, in the order written by the frozen manifest. It waits for the
entire sample lifecycle and a durable completion record before starting the
next sample. Once that record exists, the host launcher independently reruns
the single-sample validator over the frozen-manifest binding, read hashes and
sizes, paired FASTQ structure, and retained reports. Only after that succeeds
does it delete that sample's exact validated work root. A failure stops the
loop: completed checkpoints remain,
the failed sample has no trusted commit marker, and a subsequent launch
validates the frozen state and checkpoints before resuming at the first pending
sample. Nextflow `-resume` is orthogonal: the append-only
`pipeline_info/resources/resume_sessions.tsv` associates each stable invocation
key with its most recent valid session UUID. The launcher supplies
`-resume <UUID>` only when rerunning that same key; it never lets a sample,
reconciliation, or global stage borrow another stage's cache identity.
Checkpoint reconciliation remains the cross-invocation source of truth.
`CHECK_SRA_CHECKPOINTS` and `PERSIST_SRA_CHECKPOINT` use `cache false` because
the external checkpoint root is mutable state not represented in a Nextflow
`path` input hash. This prevents `-resume` from replaying stale pending rows or
skipping a required repair while leaving upstream scientific caching intact.
Before reconciliation or persistence, an exclusive
`sra_checkpoint_owner.json` claims an otherwise empty checkpoint root for one
BioProject and one frozen-manifest SHA-256. Every later reconciliation,
single-sample validation, persistence, and cleanup validates that ownership;
cross-project or changed-cohort reuse fails before any managed file is copied.
The frozen state additionally binds the normalized platform allowlist; a
different `--sra-platforms` value requires a fresh results/state root.

The final reconciliation runs with `--require-complete`, so the global cohort
cannot start while any frozen sample is absent or invalid. After the shared
global analysis, `FINALIZE_SRA_GLOBAL_RUN` first requires non-empty MultiQC
HTML, `software_versions.tsv`, and final long MAG-abundance output. Its
atomically published baseline `sra_global_success.json` binds the BioProject
and checkpoint manifest plus those three outputs by absolute path, byte size,
and SHA-256. That baseline alone cannot authorize cleanup.

Once Nextflow has finished publishing, the host checkpoint controller runs the
`seal-global` gate. It rejects symbolic links and inventories every regular file
under these durable scientific roots:

```text
01_quality_control_and_filtering/
02_mag_construction/
03_taxonomic_classification_and_phylogenomics/
04_gene_prediction_and_functional_annotation/
05_mag_abundance_estimation/
06_global_processing_evaluation/
pipeline_info/software_versions.tsv
```

The deterministic `scientific_outputs` object added atomically to the marker
contains each relative path, byte count, and SHA-256, plus aggregate file/byte
counts and the paths satisfying every mandatory deliverable family:

| Required family | Minimum durable evidence |
| --- | --- |
| final MAG catalog | one or more `final_catalog/*.fa`, plus catalog provenance and quality TSVs |
| final quality screening | final CheckM2 quality report and GUNC summary |
| species catalog | one or more 95% species representative MAG FASTAs |
| taxonomy and phylogenomics | at least one GTDB-Tk bacterial/archaeal summary and the final Newick tree |
| functional annotation | integrated annotation table for every final MAG |
| final abundance | normalized long-form table with the exact completed-checkpoint-sample x final-MAG Cartesian matrix |
| project reporting | global MultiQC HTML and `software_versions.tsv` |

Cleanup independently rebuilds that inventory immediately before deletion and
requires an exact match. A missing or extra path, changed byte count or hash,
missing required family, symbolic link, incomplete marker, or path outside the
results root fails before the first unlink. Numbered scientific outputs and the
version table are sealed; mutable launcher state, invocation logs, resource
telemetry, and cleanup provenance are durable but deliberately excluded because
they continue changing during shutdown and cannot authorize checkpoint
deletion.

Default cleanup occurs only after the complete seal and every checkpoint
artifact have been revalidated. It removes only the two exact
`<sample_id>_host_removed_R[12].fastq.gz` paths listed by validated checkpoint
rows and writes `sra_checkpoint_cleanup.json`; checkpoint records, discovery
state, reports, and provenance remain. The cleanup record is first persisted as
an `in_progress` journal containing the complete manifest-bound deletion plan,
then atomically updated after each unlink and finalized as `complete`. A restart
can therefore finish the same plan even if interruption occurred between an
unlink and its journal update. `--keep-sra-checkpoints` validates the complete
seal and checkpoint pair but retains the reads. If a baseline success marker is
found on a later launch, no scientific stage is rerun, but the launcher still
calls `seal-global` before the guarded retention/cleanup path. This recovers an
interruption between marker publication and sealing: an unsealed baseline is
extended only after the complete current result tree passes every gate, while
an existing seal is revalidated and left byte-for-byte unchanged. Calling
`cleanup` directly with an unsealed marker remains a hard error. In either
retention mode, the exact SRA global work root is removed only after seal
validation; it remains intact whenever global analysis or sealing fails.

## Compression and disk boundary

SRA Toolkit's `fasterq-dump` cannot emit gzip. The acquisition environment and
runtime checks pin SRA Toolkit 3.4.1 and pigz 2.8. Accordingly, uncompressed run
FASTQs exist only inside the configured sample-owned acquisition/fasterq
scratch. Each run is prefetched, validated with `vdb-validate`, converted, and
compressed sequentially; that run's prefetched SRA directory and fasterq
temporary directory are removed before acquisition advances. When a
BioSample has multiple runs, its compressed mates are concatenated in frozen
run order by a streaming decompress/recompress operation into one gzip pair;
the final pair is fully checked for gzip integrity, FASTQ structure, mate-name
agreement, and record count.

FastQC reads gzip inputs, fastp writes gzip, and Bowtie2 writes the durable
host-removed pair through `--un-conc-gz`. The downstream MEGAHIT, metaSPAdes,
and CoverM wrappers consume those `.fastq.gz` paths directly. Thus no durable
raw, trimmed, or host-removed uncompressed FASTQ collection is created. In SRA
mode the launcher also disables the optional publication of cleaned and
host-removed read duplicates; the sole durable large read boundary is the
external checkpoint pair. Tool-native temporary files and assembly/binning
intermediates remain task-local. In the default `save_intermediates=false`
mode, large native assembly, MetaQUAST, neural-binner, dRep, CheckM2, and GUNC
trees are removed only after their normalized outputs validate; annotation
staging/temp trees and DAS Tool auxiliaries follow the same success-only rule.
Failed tasks retain their work directory for diagnosis.

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
sample\tmag_id\trelative_abundance_percent\tmean_coverage\tcovered_fraction\tgenome_length
```

CoverM may emit a pseudo-genome row named `unmapped` for reads that did not map
to the supplied catalog. The normalizer excludes that pseudo-row because it is
not a MAG; it does not renormalize or otherwise alter CoverM's values for real
MAGs.

Before the scientific output inventory can be sealed, the checkpoint controller
revalidates this table semantically. Its header must exactly match the six
columns above in that order. Every metric must be finite and numeric;
`relative_abundance_percent` must be in `[0, 100]`, `mean_coverage` must be
non-negative, `covered_fraction` must be in `[0, 1]`, and `genome_length` must
be positive. Sample IDs come from the complete checkpoint manifest and MAG IDs
come from the final-catalog FASTA filenames. Pairs must be unique, and the table
must contain exactly one row for every checkpoint sample x final MAG
combination; no missing, duplicate, or out-of-cohort row can authorize
checkpoint cleanup.

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
| `disk_efficient` | queue size one; serial physical submission without adding scientific DAG edges |
| `gpu` | GPU requests and verified GPU environments for COMEBin, SemiBin2, and Vamb only |

The base production composition is `-profile <environment>,<runtime>`, for
example `local,docker` or `hpc,apptainer`. The launcher always appends
`disk_efficient` for SRA mode, also appends it for local-FASTQ mode when
`--storage-constrained` is set, and appends `gpu` for `--enable-gpu`. Local GPU
mode also selects `disk_efficient` so GPU-capable tasks do not contend for the
same device. These profiles affect scheduling and execution environments, not
channel dependencies.

The Conda profile requires both Conda and the `mamba` executable. Its cache uses
`params.conda_cache_dir` when set and otherwise falls back to the ignored
repository directory `.conda/`.

Docker, Apptainer, and Singularity production runs require a combined image with
PhyloPhlAn 3.1.1 and IQ-TREE 3.0.1 for the `PHYLOPHLAN` process. The stock
PhyloPhlAn image does not satisfy that complete runtime contract. The repository
build definition is `containers/phylophlan-iqtree/Dockerfile`; users build or
publish it and pass the OCI reference through `params.phylophlan_container`.
Apptainer and Singularity consume the same published image.

Resource selectors preserve the legacy scheduling intent. Examples include 32
CPUs/500 GB for MEGAHIT, 32 CPUs/1,800 GB for SPAdes, 24 CPUs/250 GB for DAS
Tool, 24 CPUs/100 GB for dRep and GTDB-Tk, 32 CPUs/120 GB for PhyloPhlAn, and
16 CPUs/128 GB for InterProScan. Local/HPC `resourceLimits` are configurable
ceilings; they do not alter scientific dependencies.

### Optional GPU execution

CPU execution is the default and remains the reproducibility baseline. The GPU
profile changes only tools with an explicit supported path in the pinned
version: COMEBin 1.0.4 retains CUDA visibility, SemiBin 1.5.0 switches from
`--engine cpu` to `--engine gpu`, and Vamb 5.0.4 adds `--cuda`. COMEBin and
SemiBin validate PyTorch 2.1.2 with CUDA 11.8; Vamb validates PyTorch 2.6.0
with CUDA 12.4. All three require CUDA to be usable and exactly one device to
be visible before scientific execution. MEGAHIT's historical ignored GPU flag
is not enabled, and no GPU capability is claimed for metaSPAdes or the other
production tools.

Each enabled process requests exactly one Nextflow accelerator. Docker receives
`--gpus 1`; Apptainer and Singularity receive `--nv`. On SLURM, the launcher
also requires an explicit site-specific `--slurm-gpu-gres` and combines that
GRES request with the configured account, QoS, and cluster options. The profile
selects the three pinned GPU container/Conda definitions without changing the
remaining scientific arguments or thresholds. CPU and GPU numerical results
need not be bitwise identical, so the execution mode, exact environment, GPU
model, and metrics are retained as provenance.

Local Conda, Apptainer, and Singularity GPU launches additionally require a
single, non-comma-separated `CUDA_VISIBLE_DEVICES` selector. This makes the
wrappers' one-visible-device contract explicit on multi-GPU hosts; Docker uses
its `--gpus 1` runtime selection, while SLURM supplies visibility through the
scheduler allocation.

Each of the three wrappers writes a task-specific `*.gpu_metrics.tsv`, sampled
at `params.gpuTelemetryInterval` (10 seconds by default). Its
filename and rows include the Nextflow session ID and task attempt. When
`nvidia-smi` is available, it samples only the visible device where that
identity is exposed, recording index, UUID, model, utilization, and used/total
memory at the configured interval; otherwise the header-only file preserves
the fact that observed device telemetry was unavailable. These files are
published under `pipeline_info/resources/gpu_tasks/` and joined on session,
process, sample/tag, and attempt so retries and separate invocations cannot be
silently conflated.

## Multi-invocation telemetry and resource accounting

The launcher starts one storage monitor for the complete logical project, even
when SRA mode uses many Nextflow sessions. Every invocation gets a unique
`pipeline_info/invocations/<run>_<stage>/` directory containing its Nextflow
log, raw trace, timeline, execution report, and DAG. The append-only
`pipeline_info/resources/trace_registry.tsv` records invocation ID, trace path,
Nextflow session UUID, internal stage, start/finish times, launch directory,
status, and exit code. This makes failed and empty invocations explicit and
prevents later stages from overwriting earlier telemetry.
`pipeline_info/resources/resume_sessions.tsv` separately records the last
session UUID by stable invocation key and is consulted only when the user
requests `--resume`.

`monitor_storage.py` samples allocated filesystem bytes without following
symbolic links and writes:

- `storage_usage_timeseries.tsv` for work, checkpoint, SRA cache, SRA
  acquisition scratch, fasterq temporary, results, optional database,
  `total_dynamic_bytes` (all changing categories, excluding databases), and
  `total_measured_bytes` (dynamic plus database when available);
- `task_workdir_peaks.tsv` for the largest sampled allocated size of each
  Nextflow hash work directory, including hash directories nested below the
  SRA invocation work roots;
- `task_workdir_timeseries.tsv` for every observed task directory at each
  timestamp. The global series also records the active launcher invocation and
  stage, allowing SRA preprocessing storage to be reported across its separate
  sample invocations.

The first dynamic sample is persisted before the immutable database tree is
scanned once in a background thread. The launcher uses a bounded
request/acknowledgement sample immediately before deleting completed sample
work, validated checkpoint reads, or successful global work, and takes a final
sample when stopping the monitor. Measurements are best effort: sampling can
miss a short-lived peak, a database scan may still be incomplete at shutdown,
and overlapping or unreadable roots are marked incomplete rather than
presenting the summed total as exact.

`summarize_resources.py` merges all registered raw traces, deduplicates trace
rows by session/task identity where possible, and retains cached rows for
provenance while excluding them from consumption totals. It reports four
levels:

| Level | Output and aggregation |
| --- | --- |
| task | `resource_usage_by_task.tsv`; raw allocation/utilization/I/O plus sampled task disk and GPU observations |
| process | `resource_usage_by_process.tsv`; executed/cached counts, time, CPU, memory, I/O, largest individual task-work peak, sampled concurrent process work, and GPU summaries by fully qualified process |
| outer subworkflow | `resource_usage_by_subworkflow.tsv`; the same aggregate metrics by stable scientific scope, plus sampled concurrent work and the SRA-preprocessing whole-stage dynamic peak where attributable |
| project | `resource_usage_summary.{tsv,json,html}`; outcome, counts, wall/cumulative time, CPU, memory, I/O, storage, GPU, SLURM availability, largest consumers, warnings, and limitations |

The stable outer-scope mapping is:

| Trace identity | Reporting scope |
| --- | --- |
| `QUALITY_CONTROL_AND_FILTERING` | `quality_control_and_filtering` |
| `MAG_CONSTRUCTION` | `mag_construction` |
| `TAXONOMIC_CLASSIFICATION_AND_PHYLOGENOMICS` | `taxonomic_classification_and_phylogenomics` |
| `GENE_PREDICTION_AND_FUNCTIONAL_ANNOTATION` | `gene_prediction_and_functional_annotation` |
| `MAG_ABUNDANCE_ESTIMATION` | `mag_abundance_estimation` |
| `GLOBAL_PROCESSING_EVALUATION` | `global_processing_evaluation` |
| `SRA_ACQUISITION_AND_PREPROCESSING`, `SRA_ACQUISITION`, `SRA_SAMPLE_PREPROCESSING`, or `SRA_PREPROCESSING` | `sra_acquisition_and_preprocessing` |
| direct `CHECK_SAMPLESHEET`, `RESOLVE_SRA_PROJECT`, `RESOLVE_SRA_INPUT`, `VALIDATE_SRA_MANIFEST`, `VALIDATE_SRA_PROJECT`, or `CHECK_SRA_CHECKPOINTS` task | `input_validation` |
| direct `COLLECT_VERSIONS`, `SUMMARIZE_RESOURCES`, `MONITOR_STORAGE`, or `FINALIZE_SRA_GLOBAL_RUN` task | `pipeline_reporting` |
| otherwise, an invocation stage containing `SRA` | `sra_acquisition_and_preprocessing` |
| no recognized workflow, task, or stage | `unmapped` and an explicit limitation |

Allocated CPU hours use `requested_cpus * realtime_hours`; observed CPU hours
use `(%cpu / 100) * realtime_hours` only where `%cpu` exists. Requested GPU
hours use `requested_accelerators * realtime_hours`. Peak RSS is summarized by
maximum/median/mean per task and is never summed as a pipeline-RAM estimate.
`max_task_peak_work_bytes` remains the largest observed individual task, while
`sampled_peak_concurrent_work_bytes` sums all trace-attributed task directories
in the same process or subworkflow snapshot. Those concurrent values exclude
external SRA cache/scratch/temp roots. For
`sra_acquisition_and_preprocessing`, the separate
`sampled_stage_peak_dynamic_storage_bytes` covers all configured dynamic roots
sampled while `sra-preprocess` was active. Logical project wall time spans the
earliest to latest registered invocation timestamps, while cumulative task
runtime is reported separately.

On HPC, the launcher optionally resolves trace `native_id` values with `sacct`
and preserves state, elapsed time, TotalCPU, MaxRSS, allocated CPUs, AllocTRES,
and ReqTRES in `pipeline_info/resources/slurm_accounting.tsv`. One global
60-second budget is shared by all accounting query batches, and its
availability/state/error are written
to `pipeline_info/resources/slurm_accounting.status.json`, including timeout or
an empty job set. This is best-effort scheduler evidence and is summarized
independently; it never replaces portable Nextflow observations.

## Publication boundary

`conf/modules.config` maps outputs into the numbered scientific hierarchy. The
default publication mode is `copy`. Large optional products are controlled by:

- `params.save_clean_reads`;
- `params.save_host_removed_reads`;
- `params.save_bam`;
- `params.save_intermediates`.

All processes still exchange required intermediates in Nextflow work
directories even when an optional user-facing publication is disabled.

SRA discovery/reconciliation state and the global-success marker are published
under `pipeline_info/sra/`. Per-sample small checkpoint process outputs are
also copied to the invocation telemetry directory, while the durable large
host-removed reads and their completion records are written directly beneath
the external checkpoint root. This is an intentional exception to the usual
Nextflow-only channel boundary: checkpoints are a validated handoff between
separate Nextflow invocations, not an attempt to rediscover arbitrary files
from the numbered results tree.

## Database preparation boundary

`bin/prepare_databases.sh` is intentionally outside the scientific DAG. It can
inspect, selectively prepare, resume, checksum, or explicitly rebuild external
resources without coupling normal runs to network access. It validates that
the database root is not `/`, `/tmp`, this repository, or another Git worktree,
uses staging and a lock, and writes a manifest plus a Nextflow configuration.

GeneMarkS-2 is validation-only in this mechanism. Its licensed installer and
key are never downloaded or copied. Production runs only consume the prepared
configuration; they do not invoke database preparation.
