from __future__ import annotations

import json
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
MAG_FILES = ("MAG_TEST_0001.fa", "MAG_TEST_0002.fna")
VALUES = {
    "task.cpus": 16,
    "prefix": "final_catalog",
    "phylophlan_db": "phylophlan",
    "database_name": "phylophlan",
    "phylophlan_config": "phylophlan_iqtree.cfg",
    "stage_commands": '\n'.join(
        f'ln -s "$PWD/{filename}" "genomes/{Path(filename).stem}.fna"'
        for filename in MAG_FILES
    ),
    "args": "",
}


def module_source() -> str:
    return (REPOSITORY / "modules/core/phylophlan/main.nf").read_text(encoding="utf-8")


def script_source() -> str:
    return module_source().split('    """', 1)[1].rsplit('"""', 1)[0]


def render_script(args: str = "") -> str:
    """Model the isolated GString in one pass; actual Groovy is checked below."""
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


def phylophlan_command(script: str) -> str:
    return "phylophlan " + script.split("    phylophlan ", 1)[1].split("    mapfile ", 1)[0]


class PhylophlanStaticTests(unittest.TestCase):
    def test_nextflow_bindings_and_all_bash_expansions_are_separated(self) -> None:
        source = script_source()
        for name, expression in (
            ("NF_CPUS", "task.cpus"),
            ("NF_PREFIX", "prefix"),
            ("NF_DATABASE", "phylophlan_db"),
            ("NF_DATABASE_NAME", "database_name"),
            ("NF_CONFIG", "phylophlan_config"),
        ):
            self.assertIn(f'{name}="${{{expression}}}"', source)
            self.assertEqual(source.count("${" + expression + "}"), 1)
        self.assertEqual(source.count("${args}"), 1)
        self.assertEqual(source.count("${stage_commands}"), 1)
        render_script()  # Rejects unsupported escapes and unbound Bash dollars anywhere.

    def test_iqtree_model_matches_forced_nucleotide_phylogeny(self) -> None:
        config = (REPOSITORY / "assets/phylophlan_iqtree.cfg").read_text(encoding="utf-8")
        self.assertIn("params = --quiet -m GTR+ASC", config)
        self.assertNotIn("-m LG", config)

    def test_scientific_command_and_empty_args_are_preserved(self) -> None:
        self.assertEqual(
            shlex.split(phylophlan_command(render_script())),
            [
                "phylophlan", "-i", "genomes", "-d", "${NF_DATABASE_NAME}",
                "-t", "a", "-f", "${NF_CONFIG}", "--databases_folder", "databases",
                "--genome_extension", ".fna", "--force_nucleotides", "--diversity", "high",
                "--accurate", "--nproc", "${NF_CPUS}", "--output_folder", "${NF_PREFIX}.phylophlan",
                "--verbose", "2>&1", "|", "tee", "${NF_PREFIX}.phylophlan.log",
            ],
        )

    def test_staging_and_output_declarations_are_preserved(self) -> None:
        source = module_source()
        self.assertIn('def database_name = phylophlan_db.name', source)
        self.assertIn(r'ln -s \"\$PWD/${entry.file}\" \"genomes/${entry.id}.fna\"', source)
        self.assertIn(r'gzip -cd \"${entry.file}\" > \"genomes/${entry.id}.fna\"', source)
        outputs = source.split("    output:\n", 1)[1].split("    when:", 1)[0]
        self.assertEqual(
            [line.strip() for line in outputs.strip().splitlines()],
            [
                "tuple val(meta), path('*.phylophlan.tree.nwk'), emit: tree",
                "tuple val(meta), path('*.phylophlan.alignment.fasta'), emit: alignment",
                "tuple val(meta), path('*.phylophlan'), emit: results",
                "tuple val(meta), path('*.phylophlan.log'), emit: log",
                'tuple val("${task.process}"), val(\'phylophlan\'), val(\'3.1.1\'), emit: versions',
                'tuple val("${task.process}"), val(\'iqtree\'), val(\'3.0.1\'), emit: versions_iqtree',
            ],
        )
        script = render_script()
        for suffix in ("tree.nwk", "alignment.fasta"):
            self.assertIn(f'test -s "${{NF_PREFIX}}.phylophlan.{suffix}"', script)
            self.assertNotIn(f'test -s ".phylophlan.{suffix}"', script)
        self.assertIn('if [[ "${#tree_files[@]}" -ne 1 ]]', script)
        self.assertIn('if [[ "${#alignments[@]}" -ne 1 ]]', script)

    def test_iqtree_version_check_uses_literal_dots_without_invalid_escapes(self) -> None:
        source = script_source()
        self.assertNotIn(r"3\.0\.1", source)
        self.assertIn(r'iqtree_version=\$(iqtree --version 2>&1)', source)
        self.assertIn("3[.]0[.]1", source)
        render_script()

    def test_extra_args_are_inserted_without_reinterpreting_their_values(self) -> None:
        tokens = shlex.split(phylophlan_command(render_script("--test-option 'two words' --test-path '$PATH'")))
        self.assertEqual(
            tokens[tokens.index("--verbose") + 1 : tokens.index("2>&1")],
            ["--test-option", "two words", "--test-path", "$PATH"],
        )


@unittest.skipUnless(GROOVY or (JAVA and GROOVY_CLASSPATH), "Groovy or JAVA + GROOVY_CLASSPATH is required")
class PhylophlanGroovyRenderingTests(unittest.TestCase):
    def test_full_script_block_and_mag_staging_render_with_real_groovy(self) -> None:
        args_values = ("", "--test-option 'two words' --test-path '$PATH'")

        def literal(value: str) -> str:
            return "'" + value.replace("\\", "\\\\").replace("'", "\\'") + "'"

        body = module_source().split("    script:\n", 1)[1].rsplit("\n}", 1)[0]
        harness = (
            "def scripts = [" + ", ".join(map(literal, args_values)) + "].collect { suppliedArgs ->\n"
            "def task = [cpus: 16, ext: [args: suppliedArgs]]\n"
            "def meta = [id: 'final_catalog']\n"
            "def mags = [new File('MAG_TEST_0001.fa'), new File('MAG_TEST_0002.fna')]\n"
            "def phylophlan_db = new File('phylophlan')\n"
            "def phylophlan_config = new File('phylophlan_iqtree.cfg')\n"
            "def error = { message -> throw new IllegalArgumentException(message) }\n"
            "def result = {\n" + body + "\n}.call()\nresult.toString()\n}\n"
            "print groovy.json.JsonOutput.toJson(scripts)\n"
        )
        command = [GROOVY] if GROOVY else [JAVA, "-cp", GROOVY_CLASSPATH, "groovy.ui.GroovyMain"]
        result = subprocess.run(
            [*command, "-e", harness], capture_output=True, text=True,
            timeout=int(os.environ.get("GROOVY_TEST_TIMEOUT", "60")),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), [render_script(args) for args in args_values])


MOCK_COMMANDS = r'''
iqtree() {
    printf 'IQ-TREE multicore version %s for Linux 64-bit\n' "$MOCK_IQTREE_VERSION"
}
phylophlan() {
    printf '%s\0' "$@" > captured_args
    [[ -L genomes/MAG_TEST_0001.fna && -L genomes/MAG_TEST_0002.fna ]] || return 90
    [[ -L databases/phylophlan ]] || return 91
    printf '%s\n' genomes/*.fna > staged_mags
    mkdir -p "${NF_PREFIX}.phylophlan/genomes_phylophlan"
    for ((mock_idx=0; mock_idx<MOCK_TREE_COUNT; mock_idx++)); do
        if [[ "$MOCK_EMPTY_TREE" == 1 ]]; then
            : > "${NF_PREFIX}.phylophlan/genomes_phylophlan/tree_${mock_idx}.treefile"
        else
            printf '(MAG_TEST_0001,MAG_TEST_0002);\n' > "${NF_PREFIX}.phylophlan/genomes_phylophlan/tree_${mock_idx}.treefile"
        fi
    done
    for ((mock_idx=0; mock_idx<MOCK_ALIGNMENT_COUNT; mock_idx++)); do
        if [[ "$MOCK_EMPTY_ALIGNMENT" == 1 ]]; then
            : > "${NF_PREFIX}.phylophlan/genomes_phylophlan/alignment_${mock_idx}_concatenated.aln"
        else
            printf '>MAG_TEST_0001\nATG\n>MAG_TEST_0002\nATG\n' > "${NF_PREFIX}.phylophlan/genomes_phylophlan/alignment_${mock_idx}_concatenated.aln"
        fi
    done
    printf 'mock PhyloPhlAn completed\n'
}
'''


@unittest.skipUnless(BASH, "Bash is required (or set BASH to its executable path)")
class PhylophlanBashTests(unittest.TestCase):
    def run_mock_script(self, args: str = "", **overrides):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "phylophlan").mkdir()
            (root / "phylophlan_iqtree.cfg").write_text("# synthetic config\n", encoding="utf-8")
            for filename in MAG_FILES:
                (root / filename).write_text(">synthetic_contig\nATG\n", encoding="utf-8")
            env = {
                **os.environ,
                "MOCK_IQTREE_VERSION": "3.0.1",
                "MOCK_TREE_COUNT": "1",
                "MOCK_ALIGNMENT_COUNT": "1",
                "MOCK_EMPTY_TREE": "0",
                "MOCK_EMPTY_ALIGNMENT": "0",
                **overrides,
            }
            result = subprocess.run(
                [BASH, "-c", MOCK_COMMANDS + render_script(args)],
                cwd=root, env=env, capture_output=True, text=True, timeout=15,
            )
            captured = root / "captured_args"
            argv = captured.read_bytes().decode().rstrip("\0").split("\0") if captured.exists() else []
            outputs = {
                suffix: (root / f"final_catalog.phylophlan.{suffix}").read_text(encoding="utf-8")
                for suffix in ("tree.nwk", "alignment.fasta", "log")
                if (root / f"final_catalog.phylophlan.{suffix}").is_file()
            }
            staged = root / "staged_mags"
            mag_links = staged.read_text().splitlines() if staged.exists() else []
            return result, argv, outputs, mag_links

    def test_full_mock_execution_has_correct_args_symlinks_and_prefixed_outputs(self) -> None:
        for args, extra in (("", []), ("--test-option 'two words'", ["--test-option", "two words"])):
            with self.subTest(args=args):
                result, argv, outputs, mag_links = self.run_mock_script(args)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(argv, [
                    "-i", "genomes", "-d", "phylophlan", "-t", "a", "-f", "phylophlan_iqtree.cfg",
                    "--databases_folder", "databases", "--genome_extension", ".fna", "--force_nucleotides",
                    "--diversity", "high", "--accurate", "--nproc", "16",
                    "--output_folder", "final_catalog.phylophlan", "--verbose", *extra,
                ])
                self.assertEqual(mag_links, ["genomes/MAG_TEST_0001.fna", "genomes/MAG_TEST_0002.fna"])
                self.assertEqual(outputs, {
                    "tree.nwk": "(MAG_TEST_0001,MAG_TEST_0002);\n",
                    "alignment.fasta": ">MAG_TEST_0001\nATG\n>MAG_TEST_0002\nATG\n",
                    "log": "mock PhyloPhlAn completed\n",
                })

    def test_iqtree_version_gate_rejects_wrong_or_missing_versions(self) -> None:
        for version in ("3.0.0", "3.0.10", "13.0.1", ""):
            with self.subTest(version=version):
                result, argv, _, _ = self.run_mock_script(MOCK_IQTREE_VERSION=version)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn("requires IQ-TREE 3.0.1", result.stderr)
                self.assertEqual(argv, [])

    def test_output_count_and_nonempty_validation_are_preserved(self) -> None:
        for overrides, error in (
            ({"MOCK_TREE_COUNT": "0"}, "Expected one IQ-TREE treefile, found 0"),
            ({"MOCK_TREE_COUNT": "2"}, "Expected one IQ-TREE treefile, found 2"),
            ({"MOCK_ALIGNMENT_COUNT": "0"}, "Expected one PhyloPhlAn concatenated alignment, found 0"),
            ({"MOCK_ALIGNMENT_COUNT": "2"}, "Expected one PhyloPhlAn concatenated alignment, found 2"),
            ({"MOCK_EMPTY_TREE": "1"}, ""),
            ({"MOCK_EMPTY_ALIGNMENT": "1"}, ""),
        ):
            with self.subTest(overrides=overrides):
                result, _, _, _ = self.run_mock_script(**overrides)
                self.assertEqual(result.returncode, 1, result.stderr)
                if error:
                    self.assertIn(error, result.stderr)


if __name__ == "__main__":
    unittest.main()
