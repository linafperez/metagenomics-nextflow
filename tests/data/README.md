# Test sequencing data

No biological data belong in this tracked directory. The deterministic fixture
generator writes paired FASTQ files and `samplesheet.csv` to the ignored
`tests/generated_data/` directory:

```bash
python3 tests/scripts/generate_synthetic_data.py
```

The launcher invokes the generator automatically for production stubs and the
scoped `--test-local` real-tool workflow. The fixtures are synthetic and are not
a substitute for production metagenomic reads.
