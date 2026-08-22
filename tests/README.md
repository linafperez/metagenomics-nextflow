# Progressive testing plan

Phase 1 uses the real pipeline modules and user-supplied small inputs. No test
data, reference genome, or Bowtie2 index is downloaded automatically.

The intended progression is:

1. Validate samplesheet parsing and normalized channel inputs.
2. Run the reusable FastQC module on raw paired reads.
3. Run raw FastQC followed by fastp.
4. Reuse FastQC on the cleaned reads after fastp.
5. Run the complete `QUALITY_CONTROL_AND_FILTERING` subworkflow after an
   existing human Bowtie2 index is supplied.

The current `test` profile targets the complete fifth test. The earlier tests
will be exercised while diagnosing the vertical slice with the same real
modules; no simplified test-only bioinformatics processes will be introduced.

See `tests/data/README.md` and `tests/reference/README.md` for the required
user-provided files.
