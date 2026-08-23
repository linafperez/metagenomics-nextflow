# Test Bowtie2 index

No reference data or indexes belong in this tracked directory. The deterministic
fixture generator writes `tests/generated_reference/GRCh38.p14.fa` and
structured stub index fixtures below the ignored `tests/generated_reference/`
directory:

```bash
python3 tests/scripts/generate_synthetic_data.py
```

For `--test-local`, `tests/workflows/synthetic_real.nf` builds a genuine tiny
Bowtie2 index from that FASTA before host-removal testing. Production stubs use
the generated structural index fixtures. Neither mode downloads GRCh38 or a
production Bowtie2 index.
