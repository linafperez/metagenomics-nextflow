# SRA Toolkit image

This image provides NCBI SRA Toolkit 3.4.1 and pigz 2.8 for the per-run
`SRATOOLS_ACQUIRE` and per-BioSample `MERGE_SRA_RUNS` processes. Acquisition uses
Bash/AWK and SRA Toolkit directly; it does not use a Python sample controller.
The environment explicitly supplies gawk. Existing images remain compatible if
they provide the pinned toolkit/pigz versions and an `awk` executable, even if
they also contain Python. Check this before using a previously cached image;
rebuild the recipe if it lacks AWK. No image was built or pulled during refactoring.

```bash
docker build -t ghcr.io/linafperez/metagenomics-sra-tools:3.4.1-pigz-2.8 containers/sra-tools
```

Build/pull operations are separate from static validation. The metadata-only
`FETCH_SRA_RUNINFO` task has its own Python environment and container.
