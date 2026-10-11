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
BASH = os.environ.get("BASH") or (shutil.which("bash") if os.name != "nt" else None)
GROOVY = shutil.which("groovy")
GROOVY_CLASSPATH = os.environ.get("GROOVY_CLASSPATH")
JAVA = shutil.which("java")
VALUES = {
    "mag": "MAG_TEST_0001.fa",
    "task.cpus": 4,
    "prefix": "MAG_TEST_0001",
    "mag_id": "MAG_TEST_0001",
    "genemark_home": "genemark_home",
    "genemark_key": "genemark_key",
    "args": "",
}


def module_source() -> str:
    return (REPOSITORY / "modules/core/genemarks2/main.nf").read_text(encoding="utf-8")


def script_source() -> str:
    return module_source().split('    """', 2)[1]


def render_script(args: str = "") -> str:
    """Model this GString's escapes in one pass, without reinterpreting values.

    Unknown escapes and unbound interpolations fail explicitly. This portable
    check is supplemented by real Groovy and Bash tests when installed.
    """
    values = {**VALUES, "args": args}
    escapes = {"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"', "$": "$", "\n": ""}

    def replace(match: re.Match[str]) -> str:
        token = match[0]
        if token.startswith("\\"):
            if token[1:] not in escapes:
                raise AssertionError(f"Unsupported Groovy escape: {token!r}")
            return escapes[token[1:]]
        if token.startswith("${"):
            return str(values[token[2:-1]])
        raise AssertionError(f"Unescaped Bash dollar sign: {token!r}")

    return re.sub(r"\\[\s\S]|\$\{[^}]+\}|\$", replace, script_source())


def prediction_command(script: str) -> str:
    return script.split('    cd "${GENEMARK_USER_DIR}"\n', 1)[1].split(
        '    cd "${output_root}"', 1
    )[0]


def version_check(script: str) -> str:
    return "    gms2_usage=" + script.split("    gms2_usage=", 1)[1].split("    mag_path=", 1)[0]


class GeneMarkS2StaticTests(unittest.TestCase):
    def test_version_comes_from_usage_without_regex_or_version_option(self) -> None:
        source = script_source()
        self.assertNotIn("--version", source)
        self.assertNotIn(r"1\.15", source)
        self.assertIn(r'gms2_usage=\$(perl "\${gms2_script}" 2>&1 || true)', source)
        self.assertIn('grep \'^Version:\' <<< "\\${gms2_usage}"', source)
        self.assertIn(r'if [[ "\${detected_version}" != 1.15* ]]; then', source)

    def test_nextflow_values_are_bound_once_and_bash_expansions_are_escaped(self) -> None:
        source = script_source()
        bindings = (
            ("NF_MAG", "mag"),
            ("NF_CPUS", "task.cpus"),
            ("NF_PREFIX", "prefix"),
            ("NF_MAG_ID", "mag_id"),
            ("NF_GENEMARK_HOME", "genemark_home"),
            ("NF_GENEMARK_KEY", "genemark_key"),
        )
        for name, expression in bindings:
            self.assertIn(f'{name}="${{{expression}}}"', source)
            self.assertEqual(source.count("${" + expression + "}"), 1)
        self.assertEqual(source.count("${args}"), 1)
        render_script()  # Also rejects unescaped Bash expansions anywhere in the block.

    def test_prediction_arguments_and_outputs_are_preserved(self) -> None:
        source = module_source()
        self.assertIn('def args   = task.ext.args ?: \'\'', source)
        self.assertIn(
            "tuple val(meta), path('*.proteins.faa'), path('*.cds.fna'), path('*.genes.gff'), emit: predictions",
            source,
        )
        self.assertIn("tuple val(meta), path('*.genemarks2.log'), emit: log", source)
        command = prediction_command(render_script())
        self.assertEqual(
            shlex.split(command),
            [
                "perl", "${gms2_script}", "--seq", "${mag_path}",
                "--genome-type", "bacteria", "--threads", "${NF_CPUS}",
                "--format", "gff", "--output", "${output_root}/${NF_PREFIX}.raw.gff",
                "--fnn", "${output_root}/${NF_PREFIX}.raw.fnn",
                "--faa", "${output_root}/${NF_PREFIX}.raw.faa",
                ">", "${output_root}/${NF_PREFIX}.genemarks2.log", "2>&1",
            ],
        )
        self.assertIn('mag_path=$(readlink -f "${NF_MAG}")', render_script())
        self.assertEqual(script_source().count(r'awk -v MAG="\${NF_MAG_ID}"'), 3)

    def test_extra_arguments_are_inserted_once_without_reinterpolation(self) -> None:
        args = "--test-option 'two words' --test-path '$PATH'"
        tokens = shlex.split(prediction_command(render_script(args)))
        self.assertEqual(tokens[-7:-3], ["--test-option", "two words", "--test-path", "$PATH"])


@unittest.skipUnless(GROOVY or (JAVA and GROOVY_CLASSPATH), "Groovy or JAVA + GROOVY_CLASSPATH is required")
class GeneMarkS2GroovyRenderingTests(unittest.TestCase):
    def test_real_gstring_rendering_matches_portable_rendering(self) -> None:
        for args in ("", "--test-option 'two words' --test-path '$PATH'"):
            with self.subTest(args=args):
                values = {**VALUES, "args": args}
                bindings = []
                for name, value in values.items():
                    if name == "task.cpus":
                        bindings.append(f"def task = [cpus: {value}]")
                    else:
                        literal = str(value).replace("\\", "\\\\").replace("'", "\\'")
                        bindings.append(f"def {name} = '{literal}'")
                harness = "\n".join(bindings) + '\nprint """' + script_source() + '"""\n'
                command = [GROOVY] if GROOVY else [JAVA, "-cp", GROOVY_CLASSPATH, "groovy.ui.GroovyMain"]
                result = subprocess.run([*command, "-e", harness], capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, render_script(args))


@unittest.skipUnless(BASH, "Bash is required (or set BASH to its executable path)")
class GeneMarkS2BashTests(unittest.TestCase):
    def run_bash(self, script: str, **kwargs) -> subprocess.CompletedProcess[str]:
        return subprocess.run([BASH, "-c", script], capture_output=True, text=True, timeout=15, **kwargs)

    def test_full_rendered_script_has_valid_bash_syntax(self) -> None:
        result = subprocess.run([BASH, "-n"], input=render_script(), capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_usage_version_acceptance_and_rejection(self) -> None:
        for version, accepted in (("1.15", True), ("1.15_1.25_lic", True), ("1.14", False), ("1.16", False), ("", False)):
            with self.subTest(version=version):
                usage = f"Usage: gms2.pl --seq SEQ --genome-type TYPE\nVersion: {version}\n" if version else "Usage: gms2.pl\n"
                stub = (
                    "set -euo pipefail\ngms2_script=mock_gms2.pl\n"
                    "perl() { [[ $# -eq 1 ]] || return 97; printf '%s' \"$MOCK_USAGE\" >&2; return 1; }\n"
                )
                result = self.run_bash(stub + version_check(render_script()), env={**os.environ, "MOCK_USAGE": usage})
                self.assertEqual(result.returncode, 0 if accepted else 1, result.stderr)
                if not accepted:
                    self.assertIn(f"detected: {version or 'unknown'}", result.stderr)

    def test_synthetic_mag_threads_filenames_and_empty_or_extra_args(self) -> None:
        for args, extra in (("", []), ("--test-option 'two words'", ["--test-option", "two words"])):
            with self.subTest(args=args), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve()
                (root / str(VALUES["mag"])).write_text(">synthetic_contig\nATG\n", encoding="utf-8")
                script = render_script(args)
                bindings = script.split('    if [[ ! -d "${NF_GENEMARK_HOME}"', 1)[0]
                stub = (
                    'gms2_script=mock_gms2.pl\noutput_root="$PWD"\n'
                    'mag_path=$(readlink -f "${NF_MAG}")\n'
                    "perl() { printf '%s\\0' \"$@\" > captured_args; }\n"
                )
                result = self.run_bash(bindings + stub + prediction_command(script), cwd=root)
                self.assertEqual(result.returncode, 0, result.stderr)
                actual = (root / "captured_args").read_bytes().decode().rstrip("\0").split("\0")
                self.assertEqual(
                    actual,
                    [
                        "mock_gms2.pl", "--seq", str(root / "MAG_TEST_0001.fa"),
                        "--genome-type", "bacteria", "--threads", "4", "--format", "gff",
                        "--output", str(root / "MAG_TEST_0001.raw.gff"),
                        "--fnn", str(root / "MAG_TEST_0001.raw.fnn"),
                        "--faa", str(root / "MAG_TEST_0001.raw.faa"), *extra,
                    ],
                )
                self.assertTrue((root / "MAG_TEST_0001.genemarks2.log").is_file())


if __name__ == "__main__":
    unittest.main()
