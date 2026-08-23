# Shotgun metagenomics pipeline

This repository implements a production-oriented shotgun metagenomics pipeline
in Nextflow DSL2. It starts from paired-end FASTQ files, removes human reads,
reconstructs metagenome-assembled genomes (MAGs) independently with MEGAHIT and
metaSPAdes, refines and combines both catalogs, and produces taxonomy,
phylogenomics, functional annotation, abundance, and global processing reports.

The normal production workflow always runs both assemblers, all four binners,
DAS Tool, both assembler-specific refinement paths, and final cross-assembler
catalog selection. The 15 selectable assembler/binner combinations are isolated
under `tests/` and never alter production behavior.

## Workflow

```text
paired-end FASTQ
  -> FastQC (raw)
  -> fastp
  -> FastQC (clean)
  -> Bowtie2 host removal
  -> paired non-host reads
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
   relative abundance, mean coverage, covered fraction, and genome length.
9. **Global evaluation.** MultiQC collects genuine native parser inputs for
   FastQC, fastp, Bowtie2, MEGAHIT, QUAST/MetaQUAST, CheckM2, and GTDB-Tk.
   Outputs from tools without a dependable parser remain available in their
   native result directories and are not claimed as automatically parsed.

## Pinned software

| Component | Version | Role or compatibility note |
| --- | ---: | --- |
| Nextflow | 26.04.6 or newer | DSL2 execution engine |
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

- Linux or WSL2 with Bash, Python 3, Java 17 or newer, and Nextflow 26.04.6 or
  newer.
- Read access to input FASTQ files and production databases.
- Sufficient storage for Nextflow work files, container/Conda caches, databases,
  coassemblies, alignments, and results.
- One supported software runtime: Docker, Conda with Mamba installed, Apptainer,
  or compatible Singularity. The `conda` profile requires the `mamba`
  executable; a Conda installation without Mamba is insufficient.

Local execution uses the Nextflow local executor. Default production resource
requests preserve the historical workflow, including 500 GB for MEGAHIT and
1,800 GB for SPAdes. A workstation that cannot satisfy those requests should
run stubs or a deliberately scoped synthetic workflow; do not lower production
resources without assessing the dataset.

HPC execution additionally requires a SLURM cluster, a shared filesystem, a
site-approved Conda or container runtime, and configured account/queue/QoS
values where the site requires them. Obsolete hostnames are not embedded in the
pipeline.

## Input samplesheet

Production starts from already available paired-end FASTQ files. SRA download
is outside the workflow.

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
- Mode: `--run`, `--stub`, `--test-local`, or `--test-hpc`.
- Database preparation: `--prepare-databases`; it is separate and takes no
  runtime flag.

| Environment | Runtime | `--run` | `--stub` | `--test-local` | `--test-hpc` |
| --- | --- | :---: | :---: | :---: | :---: |
| local | Docker | yes* | yes | yes | no |
| local | Conda + Mamba | yes | yes | yes | no |
| local | Apptainer | yes* | yes | yes | no |
| local | Singularity | yes* | yes | yes | no |
| HPC/SLURM | Docker | no | no | no | no |
| HPC/SLURM | Conda + Mamba | yes | yes | no | yes |
| HPC/SLURM | Apptainer | yes* | yes | no | yes* |
| HPC/SLURM | Singularity | yes* | yes | no | yes* |

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

Stub examples traverse the production channel graph without real biological
databases. They generate ignored synthetic fixtures and structured module
outputs. Stubs validate process contracts and orchestration, not the real tool
commands or databases.

```bash
./metagenomics_pipeline.sh --local --docker --stub
./metagenomics_pipeline.sh --local --conda --stub
./metagenomics_pipeline.sh --local --apptainer --stub
```

`--test-local` generates deterministic synthetic reads, builds a real tiny
Bowtie2 index, and selects the isolated `tests/workflows/synthetic_real.nf`
entrypoint. That workflow runs real FastQC, fastp, Bowtie2, MEGAHIT, SPAdes, and
MetaQUAST processes under conservative limits. It deliberately stops before
database-heavy and licensed stages, so it does not require or accept
`--database-config`. Its results are written below
`tests/results/synthetic_real/`.

```bash
./metagenomics_pipeline.sh --local --conda --test-local
```

`--test-hpc` is a real production run on SLURM. It requires the generated
database configuration and an explicit samplesheet; it never prepares or
redownloads databases.

```bash
./metagenomics_pipeline.sh --hpc --apptainer --test-hpc \
    --database-config /shared/db/metagenomics_databases.config \
    --phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1 \
    --input /shared/project/samplesheet.csv \
    --outdir /shared/project/test_hpc_results
```

Use `--resume` to add Nextflow `-resume`, `--dry-run` to inspect the translated
command, and `--` to stop launcher option parsing. Run
`./metagenomics_pipeline.sh --help` for the complete interface.

### Direct Nextflow execution

The launcher is optional. With a generated database configuration:

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

Direct production stub traversal uses the same profile order as the launcher:

```bash
python3 tests/scripts/generate_synthetic_data.py
nextflow run . -profile local,docker,test,stub \
    -stub-run --stub_run true
```

### Important parameters

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `--outdir` | `results` | User-facing output root |
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

## Results

Processes exchange files through Nextflow channels and work directories. No
downstream stage reads its inputs back from the published results directory.

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
```

Key deliverables include:

- final MAG FASTA files plus provenance and quality tables;
- final CheckM2 and GUNC summaries;
- 99% and 95% dRep cluster tables and representative sets;
- GTDB-Tk bacterial and archaeal classifications;
- PhyloPhlAn concatenated alignment and final Newick tree;
- per-MAG GeneMarkS-2 predictions, eggNOG and InterProScan native outputs, and
  integrated functional tables;
- CoverM wide abundance output and normalized long-form abundance table;
- `global_processing_evaluation.multiqc.html` and its MultiQC data directory;
- Nextflow report, timeline, trace, DAG, normalized samplesheet, and
  `software_versions.tsv` under `pipeline_info/`.

## Testing and validation

Generated FASTQ, reference indexes, database placeholders, test results, work
directories, and runtime caches are ignored by Git.

Static checks:

```bash
nextflow lint
python3 -m compileall -q bin tests/scripts
bash -n metagenomics_pipeline.sh bin/prepare_databases.sh
python3 -m unittest discover -s tests/scripts -p 'test_*.py'
```

Production stubs exercise structured downstream contracts without restricted
software or large databases:

```bash
./metagenomics_pipeline.sh --local --docker --stub
./metagenomics_pipeline.sh --local --conda --stub
./metagenomics_pipeline.sh --local --apptainer --stub
```

Backend availability in the configuration is not itself proof of a successful
runtime test. In the current development environment, a running Docker daemon,
Conda with Mamba, Apptainer/Singularity, all production databases, a
GeneMarkS-2 license, and SLURM/HPC access have not all been simultaneously
available. Container production also requires building or publishing the
combined PhyloPhlAn/IQ-TREE image described above.
Therefore complete real-data execution, database-heavy real-tool execution,
licensed GeneMarkS-2 validation, and HPC validation remain external pending
checks. Stubs must not be reported as real tool success.

## 15-configuration benchmark

The benchmark harness is under `tests/` and reuses production assembly, binning,
refinement, and final-catalog components. It evaluates:

- assembler: `megahit`, `spades`, or `both`;
- binning strategy: `comebin`, `metabat2`, `semibin2`, `vamb`, or `all`.

This produces 3 × 5 = 15 variants. An individual-binner variant sends that
binner's bins directly to the normal refinement path and does not run DAS Tool.
The `all` strategy runs all four binners and DAS Tool. The `both` assembler
variants retain independent assembler branches and then apply the normal final
combined-catalog logic.

Generate ignored fixtures and traverse all 15 graphs with stubs:

```bash
python3 tests/scripts/generate_synthetic_data.py
python3 tests/scripts/run_variants.py --local --stub
```

Validate a particular stub backend or one variant:

```bash
python3 tests/scripts/run_variants.py --local --docker --stub
python3 tests/scripts/run_variants.py --local --conda --stub \
    --variant both_all
python3 tests/scripts/run_variants.py --local --apptainer --stub \
    --jobs 2
```

Run the real matrix later on SLURM using filtered non-host reads and prepared
CheckM2/GUNC databases:

```bash
python3 tests/scripts/run_variants.py \
    --hpc --apptainer --run \
    --input /shared/project/filtered_reads.csv \
    --checkm2-db /shared/db/checkm2/1.1.0/uniref100.KO.1.dmnd \
    --gunc-db /shared/db/gunc/1.0.6/progenomes_2.1/gunc_db_progenomes2.1.dmnd \
    --results-root /shared/project/benchmark_results \
    --work-root /scratch/project/benchmark_work \
    --jobs 3 --resume
```

Each run writes `tests/results/<assembler>_<binner>/` by default. The runner
then calls `tests/scripts/summarize_variants.py`, which creates:

- `variant_comparison.tsv`;
- `variant_ranking.tsv`;
- `variant_summary.md`;
- `variant_comparison.json`.

The summarizer extracts native assembly, binning, CheckM2, GUNC, dRep, GTDB,
MultiQC-link, and Nextflow trace metrics when they exist. Missing values remain
`NA` or JSON `null`. It ranks variants lexicographically, with no opaque score:

1. more strict high-quality MAGs;
2. higher median completeness;
3. lower median contamination;
4. fewer GUNC failures;
5. more final non-redundant MAGs;
6. lower wall-clock runtime as a tie-breaker.

The ranking describes measured outputs only and does not imply unmeasured
biological superiority. See [tests/README.md](tests/README.md) for the focused
benchmark interface and [docs/architecture.md](docs/architecture.md) for code
organization and data contracts.

## Repository policy

- `modules/core/` contains project-maintained reusable bioinformatics modules.
- `modules/local/` and `bin/` contain transparent pipeline-specific helpers.
- nf-core is an implementation reference only; no `modules/nf-core` code is
  vendored.
- Generated reads, databases, containers, work directories, and results are not
  committed.
- GeneMark licenses, private site configuration, and machine-specific
  samplesheets must remain outside Git.
