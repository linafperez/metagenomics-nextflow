from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
AWK = os.environ.get("AWK") or shutil.which("awk")
PREFIXES = ("megahit_coassembly", "spades_coassembly")


def abundance_command(prefix: str, sample_ids: list[str]) -> list[str]:
    """Extract the real command and render its isolated Groovy string escapes.

    This models only the escapes used in this fragment; it never runs Nextflow,
    Groovy or CoverM. Rendering before shell tokenization catches a literal
    newline introduced into an AWK string by a single Groovy backslash.
    """
    source = (REPOSITORY / "modules/core/coverm/main.nf").read_text(encoding="utf-8")
    contig_process = source.split("process COVERM_GENOME", 1)[0]
    match = re.search(
        r"(?ms)^\s*(awk -v expected=.*?^\s*' .*?vamb_abundance\.tsv\")",
        contig_process,
    )
    if match is None:
        raise AssertionError("COVERM_CONTIG abundance AWK command not found")
    escapes = {"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"', "$": "$"}
    rendered = re.sub(r'\\([ntr\\"$])', lambda item: escapes[item[1]], match[1])
    rendered = rendered.replace("${expected_columns}", str(len(sample_ids) + 1))
    rendered = rendered.replace("${vamb_header}", r"\t".join(["contigname", *sample_ids]))
    rendered = rendered.replace("${prefix}", prefix)
    return shlex.split(rendered)


@unittest.skipUnless(AWK, "awk is required (or set AWK to its executable path)")
class CovermContigAwkTests(unittest.TestCase):
    def run_fragment(
        self, prefix: str, sample_ids: list[str], content: str
    ) -> tuple[subprocess.CompletedProcess[str], str]:
        command = abundance_command(prefix, sample_ids)
        self.assertEqual(command[0], "awk")
        self.assertEqual(command[-2:], [">", f"{prefix}.vamb_abundance.tsv"])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / f"{prefix}.mean_depth.raw.tsv").write_text(content, encoding="utf-8")
            output = root / command[-1]
            with output.open("w", encoding="utf-8") as handle:
                result = subprocess.run(
                    [str(AWK), *command[1:-2]],
                    cwd=root,
                    stdout=handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
            return result, output.read_text(encoding="utf-8")

    def test_normalizes_header_and_preserves_mean_depth_rows(self) -> None:
        for prefix in PREFIXES:
            for count in (1, 2, 8):
                with self.subTest(prefix=prefix, samples=count):
                    samples = [f"sample_{index}" for index in range(1, count + 1)]
                    raw_header = "\t".join(["Contig", *(f"{sample}.bam Mean" for sample in samples)])
                    rows = "".join(
                        "\t".join([contig, *([depth] * count)]) + "\n"
                        for contig, depth in (("contig_1", "0.000"), ("contig_2", "12.345"))
                    )
                    result, output = self.run_fragment(prefix, samples, raw_header + "\n" + rows)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(result.stderr, "")
                    self.assertEqual(output, "\t".join(["contigname", *samples]) + "\n" + rows)

    def test_rejects_too_few_or_too_many_header_columns(self) -> None:
        samples = [f"sample_{index}" for index in range(1, 9)]
        for prefix in PREFIXES:
            for actual in (8, 10):
                with self.subTest(prefix=prefix, columns=actual):
                    content = "\t".join(["Contig", *(["sample.bam Mean"] * (actual - 1))]) + "\n"
                    result, output = self.run_fragment(prefix, samples, content)
                    self.assertEqual(result.returncode, 1)
                    self.assertEqual(
                        result.stderr,
                        f"Unexpected CoverM mean-depth column count: {actual} (expected 9)\n",
                    )
                    self.assertEqual(output, "")


if __name__ == "__main__":
    unittest.main()
