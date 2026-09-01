# SRA preprocessing image

This image pins NCBI SRA Toolkit 3.4.1, pigz 2.8, and Python 3.12.11 for
`SRA_ACQUIRE`. Build it before production use and publish it under the value of
`params.sraContainer` (default: `metagenomics/sra-tools:3.4.1-pigz-2.8`).

```bash
docker build -t metagenomics/sra-tools:3.4.1-pigz-2.8 containers/sra-tools
```

Building or pulling this image is intentionally not part of repository tests.
