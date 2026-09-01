# Shotgun metagenomics pipeline

This repository implements a production-oriented shotgun metagenomics pipeline
in Nextflow DSL2. It accepts either an existing paired-end FASTQ samplesheet or
an NCBI SRA BioProject accession, removes human reads, reconstructs
metagenome-assembled genomes (MAGs) independently with MEGAHIT and metaSPAdes,
refines and combines both catalogs, and produces taxonomy, phylogenomics,
functional annotation, abundance, and global processing reports.

The normal production workflow always runs both assemblers, all four binners,
DAS Tool, both assembler-specific refinement paths, and final cross-assembler
catalog selection.

## Workflow

```text
local mode: paired-end FASTQ samplesheet -> raw pair ------------------|
                                                                      |
BioProject mode: frozen manifest -> one biological sample             |
  -> acquire/validate each run -> immediate gzip -> merged raw pair --|
                                                                      v
             FastQC (raw) -> fastp -> FastQC (clean) -> Bowtie2 host removal
                                                                      |
             BioProject mode: atomic durable checkpoint, then next sample
                                                                      v
             paired non-host reads (SRA: only after full reconciliation)
       |-- MEGAHIT coassembly -> MetaQUAST
       |     -> CoverM coverage
       |     -> COMEBin ---------|
       |     -> MetaBAT2 --------|-> DAS Tool
       |     -> SemiBin2 --------|
       |     -> Vamb ------------|
       |     -> CheckM2 + GUNC -> strict HQ selection
       |     -> dRep 99% ANI -> dRep 95% species groups
       |     -> CheckM2 + GUNC
       |
       `-- metaSPAdes coassembly -> MetaQUAST
             -> CoverM coverage
             -> COMEBin ---------|
             -> MetaBAT2 --------|-> DAS Tool
             -> SemiBin2 --------|
             -> Vamb ------------|
             -> CheckM2 + GUNC -> strict HQ selection
             -> dRep 99% ANI -> dRep 95% species groups
             -> CheckM2 + GUNC
                         |
                         `-> provenance-safe catalog merge
                             -> final dRep 99% ANI representatives
                             -> final dRep 95% species groups
                             -> CheckM2 + GUNC
                             -> FINAL MAG CATALOG
                                  |-- GTDB-Tk taxonomy
                                  |-- PhyloPhlAn + IQ-TREE phylogenomics
                                  |-- GeneMarkS-2
                                  |     |-- eggNOG-mapper
                                  |     `-- InterProScan
                                  |          -> integrated annotations
                                  `-- CoverM per-sample MAG abundance

all supported reports -> MultiQC
all module records     -> software_versions.tsv
```

Independent branches remain independent in the Nextflow graph. In particular,
MEGAHIT and SPAdes can run concurrently; COMEBin, MetaBAT2, SemiBin2, and Vamb
can run concurrently within each assembly branch; and final-catalog taxonomy,
phylogenomics, functional annotation, and abundance can run concurrently.
Executor limits, not artificial scientific dependencies, control actual
concurrency.

### Scientific behavior

1. **Quality control and filtering.** The same FastQC module evaluates raw and
   cleaned reads. fastp uses Q30 filtering, a 30 nt minimum length, paired-end
   adapter detection, poly-G and poly-X trimming, and right-side sliding-window
   cutting with a four-base window and Q30 mean. Bowtie2 uses
   `--very-sensitive-local --phred33` against an externally supplied
   GRCh38.p14 index and retains paired non-host reads. FastQC is diagnostic and
   does not create an automatic pass/fail branch.
2. **Assembly.** Reads from every sample are coassembled separately with
   MEGAHIT (`meta-large`, 1,000 bp assembly minimum) and metaSPAdes (`--meta`).
   The same MetaQUAST module evaluates both assemblies. Before binning, a common
   1,500 bp downstream contig threshold is applied by default.
3. **Binning.** CoverM creates reusable BAM, MetaBAT2 depth, and Vamb abundance
   representations. COMEBin, MetaBAT2, SemiBin2, and Vamb run independently;
   DAS Tool integrates their four bin sets.
4. **Assembler-specific refinement.** CheckM2 and GUNC evaluate raw DAS Tool
   bins. High-quality MAG selection is strict: completeness must be greater
   than 90% and contamination must be less than 5%. GUNC results are retained
   for inspection but do not impose an undocumented exclusion rule. dRep first
   selects strain-level representatives at 99% secondary ANI and then records
   species-level groups/representatives at 95% secondary ANI. Both operations
   use 90% primary clustering and at least 30% comparable alignment coverage.
   The 99%-representative set is evaluated again with CheckM2 and GUNC.
5. **Final catalog.** The branch catalogs are renamed with unique,
   assembler-aware MAG identifiers and joined to matching quality and provenance
   rows. Final 99% representative selection and explicit 95% species grouping
   are repeated across the combined input. The production final catalog is the
   final 99%-representative set; the 95% species representatives and cluster
   table are additional outputs. Final MAGs are evaluated again with CheckM2
   and GUNC.
6. **Taxonomy and phylogenomics.** GTDB-Tk classifies final MAGs against GTDB
   release 226. PhyloPhlAn uses conserved markers and IQ-TREE 3.0.1 as its tree
   backend, producing the alignment and Newick tree. iTOL upload is intentionally
   outside the pipeline.
7. **Functional annotation.** GeneMarkS-2 predicts coding sequences and
   proteins for each final MAG. eggNOG-mapper and InterProScan independently
   annotate the predicted proteins. A local integration step consolidates
   descriptions, orthologs, COG, GO, KEGG, EC, CAZy, PFAM, and InterPro fields
   when present.
8. **MAG abundance.** CoverM maps every filtered sample to the final catalog
   using properly paired reads, at least 95% read identity, and at least 75%
   aligned-read fraction. It emits wide and analysis-ready long tables with
   relative abundance, mean coverage, covered fraction, and genome length. The
   long-table normalizer omits CoverM's `unmapped` pseudo-genome row without
   renormalizing or otherwise changing the abundance values calculated for real
   MAGs.
9. **Global evaluation.** MultiQC collects genuine native parser inputs for
   FastQC, fastp, Bowtie2, MEGAHIT, QUAST/MetaQUAST, CheckM2, and GTDB-Tk.
   Outputs from tools without a dependable parser remain available in their
   native result directories and are not claimed as automatically parsed.

## Pinned software

| Component | Version | Role or compatibility note |
| --- | ---: | --- |
| Nextflow | 26.04.6 or newer | DSL2 execution engine |
| NCBI SRA Toolkit | 3.4.1 | BioProject run acquisition and validation |
| pigz | 2.8 | Immediate parallel compression and deterministic run merging |
| FastQC | 0.12.1 | Raw and cleaned read QC |
| fastp | 1.0.1 | Paired-read trimming and filtering |
| Bowtie2 | 2.5.4 | Human-read removal |
| MEGAHIT | 1.2.9 | Metagenomic coassembly |
| SPAdes | 4.2.0 | metaSPAdes coassembly |
| MetaQUAST / QUAST | 5.3.0 | Assembly assessment |
| CoverM | 0.7.0 | Contig coverage and final MAG abundance |
| COMEBin | 1.0.4 | The requested upstream release “1.04” is packaged as 1.0.4 |
| MetaBAT2 | 2.18 | Coverage-guided binning |
| SemiBin2 | 1.5.0 | Semi-supervised binning |
| Vamb | 5.0.4 | Variational autoencoder binning |
| DAS Tool | 1.1.7 | Multi-binner refinement |
| CheckM2 | 1.1.0 | Completeness and contamination |
| GUNC | 1.0.6 | Taxonomic consistency and chimerism assessment |
| dRep | 3.6.2 | 99% dereplication and 95% species grouping |
| GTDB-Tk | 2.6.1 | GTDB release 226 classification |
| PhyloPhlAn | 3.1.1 | Phylogenomic reconstruction |
| IQ-TREE | 3.0.1 | PhyloPhlAn tree backend |
| GeneMarkS-2 | 1.15 | Licensed prokaryotic gene prediction |
| eggNOG-mapper | 2.1.13 | Protein functional annotation |
| eggNOG mapper database | 5.0.2 | Supported database for mapper 2.1.13; see note below |
| InterProScan | 5.59-91.0 | Protein signatures, GO, pathway, and InterPro annotation |
| MultiQC | 1.35 | Global supported-report aggregation |

eggNOG-mapper 2.1.13 does not support an eggNOG 6.0 mapper database. The
pipeline therefore pins the compatible eggNOG-mapper database 5.0.2 rather
than silently presenting an unvalidated eggNOG 6 configuration. This is an
explicit compatibility substitution; the mapper version remains 2.1.13.

Each module records its declared pinned version in
`results/pipeline_info/software_versions.tsv`. Bioinformatics tools are pinned
in both module Conda environments and container definitions wherever packaging
permits; site-provided runtime overrides remain the operator's responsibility.

## Requirements

Common requirements:

- Linux or WSL2 with Bash, Python 3.10 or newer, Java 17 or newer, and
  Nextflow 26.04.6 or newer.
- Read access to input FASTQ files and production databases.
- Network access to NCBI E-utilities and SRA endpoints when BioProject mode is
  used; local samplesheet mode does not require SRA access.
- Separate external durable-checkpoint and disposable-scratch roots for
  BioProject mode.
- Sufficient storage for Nextflow work files, container/Conda caches, databases,
  coassemblies, alignments, and results.
- One supported software runtime: Docker, Conda with Mamba installed, Apptainer,
  or compatible Singularity. The `conda` profile requires the `mamba`
  executable; a Conda installation without Mamba is insufficient.

Local execution uses the Nextflow local executor. Default production resource
requests preserve the historical workflow, including 500 GB for MEGAHIT and
1,800 GB for SPAdes. A workstation that cannot satisfy those requests should
use an appropriately sized HPC environment; do not lower production resources
without assessing the dataset.

HPC execution additionally requires a SLURM cluster, a shared filesystem, a
site-approved Conda or container runtime, and configured account/queue/QoS
values where the site requires them. Obsolete hostnames are not embedded in the
pipeline.

## Input modes

Every production run must select exactly one input mode. `--input` and
`--sra-project` are mutually exclusive; the launcher rejects a command that
provides both or neither.

### Existing paired FASTQ files

Use `--input` with a CSV samplesheet:

```csv
sample,fastq_1,fastq_2
sample_A,/data/reads/sample_A_R1.fastq.gz,/data/reads/sample_A_R2.fastq.gz
sample_B,reads/sample_B_R1.fq.gz,reads/sample_B_R2.fq.gz
```

The header must be exactly `sample,fastq_1,fastq_2`. Sample identifiers must be
unique and may contain letters, digits, `.`, `_`, and `-`, but must begin with a
letter or digit. Both read files are required. Supported extensions are
`.fastq`, `.fastq.gz`, `.fq`, and `.fq.gz`. Relative paths are resolved from
the samplesheet directory. See [assets/samplesheet.csv](assets/samplesheet.csv)
for a template.

`nextflow_schema.json` defines production parameter types, allowed values, and
numeric ranges. The production workflow also validates required parameters,
numeric ranges, and the relationship `derep_ani > species_ani` before launching
scientific stages. `CHECK_SAMPLESHEET` then validates and normalizes the input
table and verifies every paired FASTQ path.

### NCBI SRA BioProject

Use `--sra-project` with a `PRJNA`, `PRJEB`, or `PRJDB` accession. Discovery
queries NCBI RunInfo metadata and freezes the returned cohort before any read
acquisition. A valid cohort must contain only public runs with a download path,
`PAIRED` layout, `WGS` strategy, `METAGENOMIC` source, and a platform in the
configured allowlist (`ILLUMINA,BGISEQ` by default). Missing, duplicate,
restricted, inconsistent, or incompatible records cause a fail-closed
validation error; runs are not silently discarded to manufacture a compatible
cohort.

BioSample accession is the primary biological-sample identity. All eligible
runs with the same BioSample are grouped and merged in deterministic accession
order. If RunInfo has no valid BioSample, the resolver uses the experiment
accession and, only if that is also unavailable, the run accession; these
fallbacks are explicit in `identity_source` and `metadata_warnings` rather than
being presented as BioSample identity.

Discovery writes the raw metadata, deterministic run and sample manifests,
exclusion table, and a summary with file sizes and SHA-256 hashes under
`<outdir>/pipeline_info/sra/`:

```text
sra_project_runinfo.csv
sra_project_manifest.tsv
sra_sample_manifest.tsv
sra_project_exclusions.tsv
sra_project_summary.json
```

Every later stage validates those frozen files and their hashes. A restart with
the same results directory reuses the cohort without re-querying NCBI. Use a
different results directory when a deliberate fresh metadata resolution is
required.

### Sequential SRA lifecycle and recovery

BioProject execution is staged by `metagenomics_pipeline.sh`:

```text
sra-discovery -> sra-checkpoints (initial) -> sra-preprocess (one per pending sample)
              -> sra-checkpoints (complete gate) -> sra-global
              -> seal every scientific output -> verified cleanup
```

These `executionStage` values are internal orchestration details; the supported
production command is the single launcher invocation shown below.

1. Resolve and freeze metadata, then reconcile existing durable checkpoints.
2. Select the next pending biological sample in deterministic manifest order.
3. For each of its runs, `prefetch` and `vdb-validate` the SRA object,
   `fasterq-dump --split-files` into the configured temporary root, immediately
   compress both mates with pigz, and remove that run's SRA and uncompressed
   temporary files before acquiring the next run.
4. Merge multi-run mates in deterministic order into one gzip pair, run raw
   FastQC, fastp, clean FastQC, and Bowtie2 host removal, and retain gzip
   throughout the persistent read path.
5. Atomically copy the non-host pair and small reports to the external
   checkpoint root. A completion record is committed only after a full gzip and
   paired-FASTQ scan, mate-name/count checks, sizes, SHA-256 hashes, and binding
   to the frozen-manifest hash all succeed. After the sample invocation returns,
   the launcher independently revalidates that record, both mates, every retained
   report, and the frozen-manifest binding. Only then is that sample's disposable
   Nextflow work directory eligible for removal. The launcher requests a storage
   sample and waits for a bounded acknowledgement before removing that exact
   directory and starting the next sample; failure of the optional monitor is
   warned but does not convert a scientifically complete sample into a failure.
6. Require a complete checkpoint reconciliation before one unchanged global
   MEGAHIT + metaSPAdes/MAG/abundance/MultiQC invocation consumes all pairs.
   Its work root is retained on failure and is eligible for safe deletion only
   after every durable scientific result has been inventoried and validated.

On restart, records and read pairs are fully revalidated and only missing or
invalid samples return to the pending loop. A partial pair without its atomic
completion record is never accepted. With `--resume`, the launcher records the
last valid Nextflow session UUID separately for each stable invocation key
(`discovery`, each sample ID, reconciliation, and `global`) in
`pipeline_info/resources/resume_sessions.tsv`; it passes `-resume <UUID>` only
back to that same key and starts a first-time key without `-resume`. If a
checkpoint disappears or changes, reconciliation and the external checkpoint
commit deliberately bypass Nextflow caching: upstream scientific tasks may
still resume, but mutable durable state is always observed and repersisted.
The frozen cohort also binds the normalized platform allowlist; changing
`--sra-platforms` requires a fresh results/state root rather than silently
reinterpreting an existing manifest. After a successful global invocation, the
launcher seals the published result set before deleting any checkpoint read.
The seal records the relative path, byte count, and SHA-256 of every regular
file under the six numbered scientific result roots, plus
`pipeline_info/software_versions.tsv`. It also requires the final MAG catalog,
catalog provenance and quality, final CheckM2 and GUNC summaries, 95% species
representatives, GTDB-Tk summary, phylogenomic tree, integrated functional
annotations, final long-form abundance, global MultiQC report, and software
versions. The abundance gate requires the exact ordered header, finite numeric
metrics in their valid ranges, unique sample/MAG pairs, and exactly one row for
the Cartesian product of completed checkpoint samples and final-catalog MAGs.
A missing, additional, changed, empty required, invalid, or symbolic-link
artifact makes validation fail closed.

If a global-success marker is already present, a rerun skips all scientific
stages but still executes the explicit `seal-global` gate before cleanup. An
already sealed marker is revalidated without being rewritten; an unsealed
baseline left by interruption is completed only if the entire current result
tree passes the same required-artifact and hash scan. Cleanup itself always
rejects an unsealed marker. If the global workflow or sealing fails, checkpoint
reads and global work are retained. After revalidation, cleanup deletes only
the exact checkpoint FASTQ files; reports, records, manifests, and cleanup
provenance remain. Before the first unlink it atomically writes the complete
deletion plan to an `in_progress` cleanup journal, updates that journal after
every file, and marks it `complete` last; an interrupted cleanup therefore
resumes the same validated plan instead of treating an already removed mate as
corruption. Pass `--keep-sra-checkpoints` to retain the
FASTQ pairs as well. Whether checkpoint reads are deleted or explicitly
retained, their validation completes first; the launcher then forces one final
storage sample and safely removes the exact SRA global work root. Local FASTQ
mode retains its normal Nextflow work directory for resume.

The seal intentionally excludes `pipeline_info/resources/`, invocation logs,
and SRA lifecycle state: telemetry and cleanup provenance continue changing
during shutdown and are not scientific inputs. They remain durable outputs but
cannot authorize deletion of checkpoint reads.

Production launches use an atomic, fail-closed ownership lock. Every mode locks
the selected results state; BioProject mode additionally locks the external
checkpoint path with a sibling lock directory, so two different result roots
cannot mutate the same checkpoint cohort. The results lock precedes results
state and telemetry initialization; after SRA paths are normalized and safety
checked, the checkpoint lock precedes checkpoint creation, storage monitoring,
and scientific stages. Locks are acquired in deterministic order and released
by the exit handler. An existing lock is never assumed stale or removed
automatically, including across HPC hosts. Inspect its owner metadata, prove
that the recorded run is no longer active, and only then remove that exact lock
directory manually.

The lock directories are
`<outdir>/pipeline_info/.metagenomics_run.lock/` and, in BioProject mode,
`<sra-checkpoint-dir>.metagenomics_run.lock/`; each publishes `owner.tsv` for
operator inspection.

SRA checkpoint, scratch, and results roots must be distinct. Each checkpoint
root must be a dedicated empty directory on first use; it is atomically sealed
to one BioProject and one exact frozen-manifest SHA-256 before any sample copy.
Reusing it for another project or cohort is rejected before existing data can
be overwritten. Checkpoint and
scratch roots are required to be outside the Git repository; scratch is
disposable, while the checkpoint root is the recovery boundary. BioProject
mode forces `publish_dir_mode=copy` and rejects symlink publication; checkpoint
manifests and global outputs bound to cleanup must be durable regular copies.
The detailed
compressed-I/O contract and upstream support evidence are recorded in
[docs/compression_audit.tsv](docs/compression_audit.tsv).

### Compression-first and storage boundaries

The persistent read path remains gzip-compressed. The only unavoidable FASTQ
decompression is the per-run output of `fasterq-dump`; it is created under the
configured temporary root, compressed immediately, and removed before the next
run is acquired. The implemented read-I/O contract is:

| Tool | Pinned version | Compressed input? | Compressed output? | Streaming? | Temporary decompression? | Implementation strategy |
| --- | ---: | --- | --- | --- | --- | --- |
| SRA Toolkit | 3.4.1 | SRA object, not FASTQ gzip | No; `fasterq-dump` emits plain FASTQ | Split paired FASTQ is not produced through stdout | Yes, one run pair only | Convert in explicit scratch, immediately pigz both mates, validate, then remove the plain FASTQ and that run's SRA object |
| pigz | 2.8 | Plain FASTQ or gzip stream | Yes | stdin/stdout supported | No additional plain files | Compress each run immediately; stream-decompress/recompress multi-run mates in frozen order |
| FastQC | 0.12.1 | Yes | Reports only | Not used | No | Read raw and trimmed `.fastq.gz` directly |
| fastp | 1.0.1 | Yes | Yes | Supported, but paired stdout is interleaved | No | Read gzip and write two named gzip mates directly |
| Bowtie2 | 2.5.4 | Yes, through the wrapper | Yes, with `--un-conc-gz` | Supported but not used | No | Read trimmed gzip and write the non-host checkpoint pair as gzip |
| MEGAHIT | 1.2.9 | Yes | Contigs/reports, not reads | Not needed | No pipeline-created FASTQ | Consume all ordered non-host gzip pairs directly |
| metaSPAdes | 4.2.0 | Yes | Contigs/reports, not reads | Not needed | No pipeline-created FASTQ | Reference ordered gzip mates in one dataset YAML |
| CoverM | 0.7.0 | Yes | BAM/coverage/abundance, not reads | Not needed | No pipeline-created FASTQ | Consume paired non-host gzip directly for contig coverage and final MAG abundance |

Storage roles are deliberately separate:

- the SRA cache holds one prefetched run at a time and is disposable;
- acquisition scratch and the `fasterq-dump` temporary root hold run-local
  conversion data and are disposable;
- the external checkpoint root is durable across invocations and contains the
  validated non-host pairs, small reports, and completion records;
- Nextflow `work/` contains task intermediates; completed SRA sample work is
  removed after checkpoint commit, failed global work is retained, and
  successful SRA global work is removed only after durable-output validation;
- with `save_intermediates=false`, large native assembly, MetaQUAST, binner,
  dRep, CheckM2, and GUNC trees plus annotation/DAS Tool scratch are pruned
  inside a successful task only after required normalized outputs validate;
- databases are immutable external inputs, never mixed with cache or work;
- results contain durable scientific outputs, manifests, provenance, and
  telemetry. Default post-success cleanup removes checkpoint read pairs only
  after every final consumer and complete scientific-output seal validation;
  retention is
  controlled by `--keep-sra-checkpoints`.

## External resources

All production resources live outside Git. Normal pipeline execution never
downloads them.

| Parameter | Required resource | Expected value |
| --- | --- | --- |
| `--host_bowtie2_index` | GRCh38.p14, GENCODE release 44 | Bowtie2 index prefix without `.1.bt2` and related suffixes |
| `--checkm2_db` | CheckM2 1.1.0 database | Database directory or its DIAMOND file |
| `--gunc_db` | GUNC 1.0.6 ProGenomes 2.1 database | GUNC `.dmnd` database file |
| `--gtdbtk_db` | GTDB-Tk release 226 | Unpacked GTDB data directory |
| `--phylophlan_db` | PhyloPhlAn 3.1.1 markers | PhyloPhlAn database directory |
| `--phylophlan_config` | PhyloPhlAn backend configuration | Defaults to `assets/phylophlan_iqtree.cfg` |
| `--genemark_home` | GeneMarkS-2 1.15 | Licensed installation containing `gms2.pl` |
| `--genemark_key` | GeneMark license key | Readable key file; never commit it |
| `--eggnog_db` | eggNOG-mapper database 5.0.2 | Directory containing `eggnog.db`, `eggnog.taxa.db`, and `eggnog_proteins.dmnd` |
| `--interproscan_data` | InterProScan 5.59-91.0 data | Data directory containing the bundled analysis databases |

### Prepare databases separately

The dedicated preparer checks existing resources, resumes HTTP downloads,
validates expected files and published checksums where available, and avoids
redownloading valid databases. It stores recoverable backups when an explicit
redownload replaces an existing resource. The database root must be outside
this repository and outside any Git worktree.

```bash
./metagenomics_pipeline.sh --hpc --prepare-databases \
    --db-root /shared/databases/metagenomics \
    --genemark-home /opt/genemarks2-1.15 \
    --genemark-key /secure/licenses/gm_key
```

Use a local GRCh38.p14 FASTA instead of downloading it:

```bash
./metagenomics_pipeline.sh --hpc --prepare-databases \
    --db-root /shared/databases/metagenomics \
    --human-reference /shared/references/GRCh38.p14.genome.fa.gz
```

Inspect without downloading, or prepare selected resources only:

```bash
./metagenomics_pipeline.sh --hpc --prepare-databases \
    --db-root /shared/databases/metagenomics --check-only

./metagenomics_pipeline.sh --hpc --prepare-databases \
    --db-root /shared/databases/metagenomics \
    --only human_reference,bowtie2_index,checkm2,gunc
```

The preparer writes:

- `/shared/databases/metagenomics/database_manifest.tsv`, containing version,
  status, path, provenance, check time, and validation details;
- `/shared/databases/metagenomics/metagenomics_databases.config`, containing
  ready resource paths for Nextflow.

The preparer invokes `bowtie2-build`, `checkm2 database --download`, and
`gunc download_db` where those official tools are required. It downloads the
three eggNOG 5.0.2 archives directly from the current official eggNOG download
host because the URL embedded in eggNOG-mapper 2.1.13 is obsolete. An optional
patched downloader can be supplied with `--eggnog-downloader`. Command overrides
are listed by `bin/prepare_databases.sh --help`.

GeneMarkS-2 is restricted software. The preparer never downloads or copies the
installation or key; it only validates the user-provided paths and writes them
to the private generated configuration. Keep the database configuration and
key outside Git with restrictive permissions.

## Execution model

Environment and software runtime are independent profile dimensions:

```text
--<environment> --<runtime> --<mode>
```

- Environment: `--local` or `--hpc`.
- Runtime: `--docker`, `--conda`, `--apptainer`, or `--singularity`.
- Mode: `--run`.
- Database preparation: `--prepare-databases`; it is separate and takes no
  runtime flag.

| Environment | Runtime | `--run` |
| --- | --- | :---: |
| local | Docker | yes* |
| local | Conda + Mamba | yes |
| local | Apptainer | yes* |
| local | Singularity | yes* |
| HPC/SLURM | Docker | no |
| HPC/SLURM | Conda + Mamba | yes |
| HPC/SLURM | Apptainer | yes* |
| HPC/SLURM | Singularity | yes* |

Docker requires a running daemon. The Conda profile requires both Conda and
Mamba and enables Mamba solving. Its environment cache defaults to the ignored
repository path `.conda/`; set `--conda_cache_dir` to use a shared or external
cache.
Apptainer is the preferred HPC container runtime; the Singularity profile is a
compatibility option for sites that still expose that executable. Nextflow uses
the same OCI image declarations for Docker, Apptainer, and Singularity.

The PhyloPhlAn Conda environment pins both PhyloPhlAn 3.1.1 and IQ-TREE 3.0.1.
The stock PhyloPhlAn container does not guarantee that exact IQ-TREE version;
the starred container production cells above are supported after the repository
combined image is built and, for distributed execution, published to a registry
available from every compute node:

```bash
docker build \
    --file containers/phylophlan-iqtree/Dockerfile \
    --tag metagenomics/phylophlan-iqtree:3.1.1-3.0.1 \
    .
```

The local tag is the Docker default. For Apptainer, Singularity, or distributed
Docker execution, retag and publish it to an accessible OCI registry, then pass
that image to the pipeline:

```bash
./metagenomics_pipeline.sh --local --docker --run \
    --database-config /data/db/metagenomics_databases.config \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --input /data/project/samplesheet.csv
```

See
[containers/phylophlan-iqtree/README.md](containers/phylophlan-iqtree/README.md)
for the exact build contract.

Container-backed BioProject execution also requires the repository SRA image,
which pins SRA Toolkit 3.4.1, pigz 2.8, and Python 3.12.11:

```bash
docker build -t metagenomics/sra-tools:3.4.1-pigz-2.8 containers/sra-tools
```

For Docker, the launcher runs every SRA lifecycle container with the invoking
host UID:GID and a writable scratch-backed `HOME`. This keeps work,
checkpoints, ownership records, and cleanup permissions usable by the host
launcher across process boundaries.

Build or publish it before production; the pipeline does not build or pull
images as part of its tests. See
[containers/sra-tools/README.md](containers/sra-tools/README.md).
For distributed Docker/Apptainer/Singularity execution, publish it to an OCI
registry visible to all nodes and pass that reference as `--sraContainer`.

### Launcher

Make the launcher executable once if needed:

```bash
chmod +x metagenomics_pipeline.sh bin/prepare_databases.sh
```

Production examples:

```bash
./metagenomics_pipeline.sh --local --docker --run \
    --database-config /data/db/metagenomics_databases.config \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --input /data/project/samplesheet.csv \
    --outdir /data/project/results

./metagenomics_pipeline.sh --local --conda --run \
    --database-config /data/db/metagenomics_databases.config \
    --input /data/project/samplesheet.csv

./metagenomics_pipeline.sh --local --apptainer --run \
    --database-config /data/db/metagenomics_databases.config \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --input /data/project/samplesheet.csv

./metagenomics_pipeline.sh --hpc --apptainer --run \
    --database-config /shared/db/metagenomics_databases.config \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --input /shared/project/samplesheet.csv \
    --slurm_account ACCOUNT --slurm_queue PARTITION --slurm_qos QOS
```

BioProject mode requires explicit, external checkpoint and scratch roots. It
automatically serializes queued tasks to reduce overlapping disk demand:

```bash
./metagenomics_pipeline.sh --hpc --apptainer --run \
    --database-config /shared/db/metagenomics_databases.config \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --sra-project PRJNA123456 \
    --sra-checkpoint-dir /shared/checkpoints/PRJNA123456 \
    --sra-scratch-dir /scratch/project/PRJNA123456 \
    --outdir /shared/project/PRJNA123456/results \
    --resource-database-root /shared/db/metagenomics \
    --slurm_account ACCOUNT --slurm_queue PARTITION --slurm_qos QOS
```

`--sra-cache-dir` defaults to `<sra-scratch-dir>/sra-cache` and
`--sra-temp-dir` to `<sra-scratch-dir>/fasterq-temp`. `--sra-email` adds an
NCBI contact address, `--sra-platforms` changes the discovery allowlist, and
`--sra-max-size` is passed to `prefetch` (default `u`, unlimited). All staged
SRA/results paths reject single or double quotes, backticks, dollar signs,
backslashes, and line breaks. Container-backed SRA storage paths additionally
reject whitespace, comma, and colon so bind arguments remain unambiguous.

The launcher always selects the `disk_efficient` profile for BioProject mode;
`--storage-constrained` selects the same profile for local-FASTQ mode. It sets
the executor queue size to one without altering scientific dependencies,
inputs, or parameters, and prevents independent large tasks such as MEGAHIT
and metaSPAdes from occupying disk concurrently. It is a scheduling control,
not a 500 GB quota or a guarantee: actual peak storage remains dataset- and
filesystem-dependent and must be checked in the generated telemetry. Enforce
the site's filesystem quota and first confirm that the largest single task can
fit, because serialization cannot make an individual assembly smaller.

Use `--resume` to reuse `-resume <UUID>` only when the same stable invocation
key already has a recorded Nextflow session; a first execution of that key
starts without `-resume`. Use `--dry-run` to inspect the translated commands.
The `--` passthrough delimiter is intentionally rejected so internal staged
parameters cannot be overridden; pass supported pipeline parameters directly.
Run `./metagenomics_pipeline.sh --help` for the complete interface.

### Optional GPU execution

CPU execution is the default and remains the reference path. `--enable-gpu`
enables only the three pinned tools with a verified upstream GPU interface:

- COMEBin 1.0.4 with CUDA-visible PyTorch 2.1.2/CUDA 11.8;
- SemiBin2 1.5.0 with `--engine gpu` and PyTorch 2.1.2/CUDA 11.8;
- Vamb 5.0.4 with `--cuda` and PyTorch 2.6.0/CUDA 12.4.

Each enabled task requests exactly one GPU; the launcher rejects any other
`--gpu-accelerators` value. Docker receives `--gpus 1`, while Apptainer and
Singularity receive `--nv`. For a local Conda, Apptainer, or Singularity run,
select exactly one device in the launch environment so the task cannot see all
GPUs on a multi-GPU host:

```bash
CUDA_VISIBLE_DEVICES=0 ./metagenomics_pipeline.sh --local --conda --run \
    --input /data/project/samplesheet.csv --enable-gpu
```

The launcher rejects an absent or comma-separated selector for those three
local runtime modes. An HPC GPU run must additionally provide the
site-specific GRES string:

```bash
./metagenomics_pipeline.sh --hpc --apptainer --run \
    --database-config /shared/db/metagenomics_databases.config \
    --input /shared/project/samplesheet.csv \
    --enable-gpu --gpu-accelerators 1 \
    --slurm-gpu-gres gpu:a100:1
```

Build and publish the pinned GPU images before using a container profile:

```bash
docker build -t metagenomics/comebin-gpu:1.0.4-cuda11.8 containers/comebin-gpu
docker build -t metagenomics/semibin-gpu:1.5.0-cuda11.8 containers/semibin-gpu
docker build -t metagenomics/vamb-gpu:5.0.4-cuda12.4 containers/vamb-gpu
```

Distributed runtimes need registry-visible images; override
`--comebinGpuContainer`, `--semibinGpuContainer`, and `--vambGpuContainer` when
the default local tags are not visible to compute nodes. The host must provide
a compatible NVIDIA driver, and Docker additionally needs the NVIDIA Container
Toolkit. The corresponding exact-version Conda environments are available for
the Conda profile.

Before each accelerated scientific command, the wrapper verifies the exact
PyTorch and CUDA runtime versions shown above, `torch.cuda.is_available()`, and
exactly one visible CUDA device. A mismatch fails that task instead of silently
running with a different runtime or falling back to CPU.

GPU training can differ numerically or stochastically from CPU execution, and
Vamb upstream does not guarantee deterministic training even with a seed. The
pipeline does not claim bitwise CPU/GPU equivalence or acceleration for any
other tool; scientific thresholds and non-device options remain unchanged. See
[docs/gpu_capability_audit.tsv](docs/gpu_capability_audit.tsv) and
[containers/gpu/README.md](containers/gpu/README.md) for the evidence and build
contract.

### Direct Nextflow execution

Direct Nextflow execution is supported for the existing-FASTQ samplesheet mode.
BioProject production must use the launcher because its separate invocations,
external checkpoints, per-sample work cleanup, and final success gate are part
of the storage and recovery contract. With a generated database configuration:

```bash
nextflow -c /data/db/metagenomics_databases.config run . \
    -profile local,docker \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --input /data/project/samplesheet.csv \
    --outdir /data/project/results
```

Without that file, provide every required database and license parameter:

```bash
nextflow run . -profile hpc,apptainer \
    --input /shared/project/samplesheet.csv \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --host_bowtie2_index /shared/db/human/bowtie2/GRCh38_p14 \
    --checkm2_db /shared/db/checkm2/1.1.0 \
    --gunc_db /shared/db/gunc/1.0.6/progenomes_2.1/gunc_db_progenomes2.1.dmnd \
    --gtdbtk_db /shared/db/gtdbtk/release226 \
    --phylophlan_db /shared/db/phylophlan/phylophlan \
    --genemark_home /opt/genemarks2-1.15 \
    --genemark_key /secure/licenses/gm_key \
    --eggnog_db /shared/db/eggnog/5.0.2 \
    --interproscan_data /shared/db/interproscan/5.59-91.0/data \
    --outdir /shared/project/results
```

### Important parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--outdir` | `results` | User-facing output root |
| `--work-dir` | `work` locally; `<sra-scratch-dir>/nextflow-work` for SRA | Nextflow work root |
| `--input` | unset | Existing FASTQ samplesheet; mutually exclusive with `--sra-project` |
| `--sra-project` | unset | BioProject accession; mutually exclusive with `--input` |
| `--sra-checkpoint-dir` | required in SRA mode | External durable non-host read checkpoint root |
| `--sra-scratch-dir` | required in SRA mode | External disposable acquisition and work root |
| `--sra-cache-dir` | `<sra-scratch-dir>/sra-cache` | SRA `prefetch` cache |
| `--sra-temp-dir` | `<sra-scratch-dir>/fasterq-temp` | Uncompressed `fasterq-dump` temporary root |
| `--sra-email` | unset | Optional NCBI E-utilities contact email |
| `--sra-platforms` | `ILLUMINA,BGISEQ` | SRA discovery platform allowlist |
| `--sra-max-size` | `u` | Maximum size value passed to `prefetch` |
| `--keep-sra-checkpoints` | false | Retain checkpoint FASTQ pairs after validated global success |
| `--storage-constrained` | false (automatic for SRA) | Select queue-size-one disk scheduling for local FASTQ mode; SRA always selects it |
| `--resource-sample-interval` | `60` | Storage sampling interval in seconds |
| `--resource-database-root` | unset | Optional database tree included in storage accounting |
| `--enable-gpu` | false | Enable only COMEBin, SemiBin2, and Vamb GPU paths |
| `--gpu-accelerators` | `1` | Per-task GPU request; implemented mode requires exactly one |
| `--gpu-telemetry-interval` | `10` | Best-effort `nvidia-smi` sampling interval in seconds |
| `--slurm-gpu-gres` | required for HPC GPU mode | Site-specific SLURM GRES value |
| `--resume` | false | Reuse the last recorded Nextflow UUID for the same stable invocation key |
| `--dry-run` | false | Print translated command(s) without execution |
| `--publish_dir_mode` | `copy` | Nextflow publication mode |
| `--conda_cache_dir` | unset (`.conda/`) | Conda environment cache; the shown repository path is the effective fallback |
| `--phylophlan_container` | `metagenomics/phylophlan-iqtree:3.1.1-3.0.1` | Combined PhyloPhlAn 3.1.1 and IQ-TREE 3.0.1 OCI image for container production |
| `--min_contig_length` | `1500` | Minimum contig length entering binning |
| `--hq_completeness` | `90.0` | Strict lower bound; MAG must be greater than this |
| `--hq_contamination` | `5.0` | Strict upper bound; MAG must be less than this |
| `--derep_ani` | `0.99` | Strain-level dRep secondary ANI |
| `--species_ani` | `0.95` | Species-level dRep secondary ANI |
| `--derep_coverage` | `0.30` | Minimum comparable alignment fraction |
| `--save_clean_reads` | `false` | Publish fastp-cleaned reads |
| `--save_host_removed_reads` | `true` | Publish final non-host reads |
| `--save_bam` | `false` | Publish binning BAMs |
| `--save_intermediates` | `false` | Publish selected large native intermediate directories |
| `--local_max_jobs` | `2` | Local executor queue size |
| `--max_cpus` | `32` | Per-task profile resource ceiling |
| `--max_memory` | `1800 GB` | Per-task profile resource ceiling |
| `--max_time` | `14d` | Per-task profile resource ceiling |
| `--slurm_account` | unset | SLURM account |
| `--slurm_queue` | unset | SLURM partition/queue |
| `--slurm_qos` | unset | SLURM QoS |
| `--slurm_cluster_options` | unset | Additional site-specific SLURM options |
| `--slurm_queue_size` | `100` | Maximum queued SLURM tasks |

In BioProject mode the launcher forces `save_clean_reads=false` and
`save_host_removed_reads=false` for normal publication because the validated
external checkpoint is the single durable large-read copy. This does not alter
the local samplesheet defaults.

## Resource telemetry and accounting

The launcher starts one storage monitor across the complete project lifecycle
and registers every Nextflow invocation, including the separate discovery,
checkpoint, per-sample, and global invocations used by BioProject mode. Each
invocation has its own log, trace, report, timeline, and DAG under
`pipeline_info/invocations/`; the registry and final merger prevent one stage
from overwriting another.

Accounting is emitted at four aggregation levels:

| Level | Primary output | Contents |
| --- | --- | --- |
| Task | `pipeline_info/resources/resource_usage_by_task.tsv` | Requested resources, runtime, CPU, RSS/VMEM, I/O, sampled task-work peak, accelerator request, GPU observations, status, invocation, and raw trace fields |
| Process | `pipeline_info/resources/resource_usage_by_process.tsv` | Executed/cached counts, requested-versus-observed resources, CPU-hours, memory efficiency, maximum individual task-work peak, sampled concurrent process-work peak, I/O, and GPU aggregates |
| Subworkflow | `pipeline_info/resources/resource_usage_by_subworkflow.tsv` | The same aggregates plus process counts, sampled concurrent work, and the sampled whole-stage dynamic peak available for SRA preprocessing |
| Project | `pipeline_info/resources/resource_usage_summary.{tsv,json,html}` | Outcome, wall time, total consumption, dynamic and per-category storage peaks, largest consumers, coverage, GPU, SLURM, warnings, and limitations |

`pipeline_info/execution_trace.tsv` is the merged provenance-preserving trace,
while `pipeline_info/resources/trace_registry.tsv` binds each source trace to
its invocation/session/stage and completion status. Cached rows remain visible
but contribute zero resource consumption; duplicate registered trace rows are
excluded. This makes resumed and multi-invocation runs additive without
counting a cached task as newly executed work.

The principal formulas use task realtime in seconds:

- allocated CPU-hours = `requested_cpus * realtime / 3600`;
- observed CPU-hours = `(%cpu / 100) * realtime / 3600`;
- observed CPU efficiency = `observed_cpu_hours / allocated_cpu_hours` for
  tasks with a CPU measurement;
- requested GPU-hours = `requested_accelerators * realtime / 3600`.

Cumulative task time is a sum and may exceed wall-clock time when tasks overlap.
Project wall time is the earliest registered invocation start through the latest
finish, with trace task bounds used only when registry timestamps are
incomplete.

`storage_usage_timeseries.tsv` samples allocated filesystem bytes without
following symbolic links for work, durable checkpoints, SRA cache, acquisition
scratch, fasterq temporary files, and results. Every row carries the current
launcher invocation and stage. It records both per-category peaks and
`total_dynamic_bytes`, which excludes the immutable database tree.
The optional database root is measured once in the background and reported
separately when that measurement completes, so a large database scan does not
delay dynamic sampling or inflate the working-storage peak.
`task_workdir_peaks.tsv` preserves the largest observed size of each Nextflow
task directory across samples. `task_workdir_timeseries.tsv` retains each
same-timestamp task-directory sample; the summarizer joins those paths to trace
tasks and sums concurrent members of a process or outer scientific scope.

`--resource-sample-interval` controls periodic sampling; the default is 60
seconds. The launcher additionally uses a bounded request/acknowledgement
sample before deleting completed SRA sample work, validated checkpoint reads,
or successful global work, and takes a final sample during shutdown.

In GPU mode each enabled process writes `gpu_tasks/*.gpu_metrics.tsv` from
best-effort `nvidia-smi` sampling at the configurable
`--gpu-telemetry-interval` (10 seconds by default). Filenames and rows carry the Nextflow session
ID and task attempt; the merger joins on session, process, sample/tag, and
attempt so retries from separate invocations are not misattributed. It reports
device models, mean and maximum utilization, peak observed device memory, and
sample coverage. On HPC, the launcher makes one globally bounded 60-second
attempt, shared by all query batches, to collect `sacct` records by native job
ID into `slurm_accounting.tsv`; collection
state and errors are persisted separately in
`slurm_accounting.status.json`. Timeout, unavailable, empty, and failed states
remain explicit. This enrichment never replaces the portable Nextflow trace.

Telemetry is observational, not a scheduler guarantee. Sampled concurrent
process/subworkflow work excludes external SRA cache/scratch/temp categories,
while the SRA-preprocessing stage peak includes all configured dynamic roots at
that stage. Sampling can miss a short disk or GPU spike; GPU utilization is
device-wide within the task's
visible device set rather than a per-PID measurement; trace RSS, CPU, I/O,
queue, host, and accelerator fields depend on executor/Nextflow availability;
`nvidia-smi` and `sacct` may be unavailable; and a forced process or machine
termination can prevent the final sample. Input FASTQ files outside the
configured roots are not included. Published output size cannot be attributed
reliably to a process or subworkflow and is therefore reported at storage/project
scope rather than fabricated at those levels. Missing measurements are reported
as limitations rather than invented as zero or used to change scientific
success.

Historical SLURM CPU/RAM requests and their current selectors are mapped in
[docs/original_resource_mapping.tsv](docs/original_resource_mapping.tsv). In
particular, the recovered metaSPAdes request is 32 CPUs and 1,800 GB RAM.
Historical scripts did not provide elapsed-time requests and requested no GPUs;
the audit records those facts explicitly as `not recorded` and `none requested`.
Current time limits are identified as policy fallbacks, not benchmarks.

## Results

Within each invocation, processes exchange files through Nextflow channels and
work directories. Separate BioProject invocations meet only at the validated
external checkpoint boundary. No downstream stage reads its inputs back from
the published results directory.

```text
results/
|-- 01_quality_control_and_filtering/
|   |-- fastqc/raw/<sample>/
|   |-- fastp/<sample>/
|   |-- fastqc/clean/<sample>/
|   `-- host_removal/<sample>/
|-- 02_mag_construction/
|   |-- megahit/
|   |   |-- assembly/metaquast/
|   |   |-- binning/{coverage,comebin,metabat2,semibin2,vamb,dastool}/
|   |   `-- refinement/{raw,high_quality,drep_99,species_95,clean}/
|   |-- spades/
|   |   |-- assembly/metaquast/
|   |   |-- binning/{coverage,comebin,metabat2,semibin2,vamb,dastool}/
|   |   `-- refinement/{raw,high_quality,drep_99,species_95,clean}/
|   `-- final_catalog/
|       |-- final_catalog/*.fa
|       |-- final_catalog.provenance.tsv
|       |-- final_catalog.quality.tsv
|       |-- drep_99/
|       |-- species_95/
|       `-- evaluation/{checkm2,gunc}/
|-- 03_taxonomic_classification_and_phylogenomics/
|   |-- gtdbtk/
|   `-- phylophlan/
|-- 04_gene_prediction_and_functional_annotation/
|   |-- genemarks2/<MAG>/
|   |-- eggnog_mapper/<MAG>/
|   |-- interproscan/<MAG>/
|   `-- integrated/<MAG>/
|-- 05_mag_abundance_estimation/
|-- 06_global_processing_evaluation/
`-- pipeline_info/
    |-- execution_trace.tsv
    |-- software_versions.tsv
    |-- invocations/<run_token>_<stage>/
    |-- resources/
    |   |-- trace_registry.tsv
    |   |-- resume_sessions.tsv
    |   |-- storage_usage_timeseries.tsv
    |   |-- task_workdir_peaks.tsv
    |   |-- task_workdir_timeseries.tsv
    |   |-- resource_usage_by_{task,process,subworkflow}.tsv
    |   |-- gpu_tasks/*.gpu_metrics.tsv          # GPU mode
    |   |-- slurm_accounting.tsv                 # best-effort HPC enrichment
    |   |-- slurm_accounting.status.json         # collection state/error
    |   `-- resource_usage_summary.{tsv,json,html}
    `-- sra/                         # BioProject mode only
        |-- sra_project_{runinfo.csv,manifest.tsv,exclusions.tsv,summary.json}
        |-- sra_sample_manifest.tsv
        |-- sra_{checkpoint_manifest.tsv,pending_samples.tsv,checkpoint_status.json}
        `-- sra_global_success.json       # full scientific-output hash inventory
```

The external SRA checkpoint root has a deliberately small, stable interface:

```text
<sra-checkpoint-dir>/
|-- sra_checkpoint_owner.json        # immutable project/cohort ownership
|-- reads/<sample>_host_removed_R{1,2}.fastq.gz
|-- reports/<sample>/...
|-- records/<sample>.checkpoint.json
`-- sra_checkpoint_cleanup.json       # after validated default cleanup
```

After normal global success, the `reads/` files are absent unless
`--keep-sra-checkpoints` was used; reports, completion records, frozen state,
and cleanup provenance remain available for audit.

Key deliverables include:

- final MAG FASTA files plus provenance and quality tables;
- final CheckM2 and GUNC summaries;
- 99% and 95% dRep cluster tables and representative sets;
- GTDB-Tk bacterial and archaeal classifications;
- PhyloPhlAn concatenated alignment and final Newick tree;
- per-MAG GeneMarkS-2 predictions, eggNOG and InterProScan native outputs, and
  integrated functional tables;
- CoverM wide abundance output and normalized long-form abundance table, with
  one validated row per completed sample and final-catalog MAG;
- `global_processing_evaluation.multiqc.html` and its MultiQC data directory;
- Nextflow report, timeline, trace, DAG, local-mode normalized samplesheet, and
  `software_versions.tsv` under `pipeline_info/`;
- the four-level resource summaries and raw storage/GPU/SLURM telemetry described
  above;
- in BioProject mode, the frozen cohort, checkpoint reconciliation, and
  global-success marker containing the complete scientific-output hash
  inventory.

## Repository policy

- `modules/core/` contains project-maintained reusable bioinformatics modules.
- `modules/local/` and `bin/` contain transparent pipeline-specific helpers.
- nf-core is an implementation reference only; no `modules/nf-core` code is
  vendored.
- Generated reads, databases, containers, work directories, and results are not
  committed.
- GeneMark licenses, private site configuration, and machine-specific
  samplesheets must remain outside Git.
