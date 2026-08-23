#!/usr/bin/env python3
"""Generate deterministic, ignored test reads and placeholder resources."""

from __future__ import annotations

import argparse
import csv
import gzip
import random
from pathlib import Path


DNA = "ACGT"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, default=Path(__file__).parents[2])
    parser.add_argument("--pairs", type=int, default=400)
    return parser.parse_args()


def random_sequence(length: int, seed: int) -> str:
    generator = random.Random(seed)
    return "".join(generator.choice(DNA) for _ in range(length))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_fasta(path: Path, name: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{name}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset : offset + 80] + "\n")


def write_reads(
    first_path: Path,
    second_path: Path,
    sample: str,
    host: str,
    microbial: str,
    pairs: int,
) -> None:
    read_length = 150
    insert_length = 280
    quality = "I" * read_length
    generator = random.Random(f"{sample}-synthetic")

    first_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(first_path, "wt", encoding="ascii", newline="\n") as first, gzip.open(
        second_path, "wt", encoding="ascii", newline="\n"
    ) as second:
        for index in range(pairs):
            source = host if index % 5 == 0 else microbial
            start = generator.randrange(0, len(source) - insert_length)
            fragment = source[start : start + insert_length]
            read_1 = fragment[:read_length]
            read_2 = reverse_complement(fragment[-read_length:])
            name = f"{sample}_{index + 1}"
            first.write(f"@{name}/1\n{read_1}\n+\n{quality}\n")
            second.write(f"@{name}/2\n{read_2}\n+\n{quality}\n")


def create_stub_resources(project_dir: Path) -> tuple[str, list[str]]:
    data_dir = project_dir / "tests" / "generated_data"
    reference_dir = project_dir / "tests" / "generated_reference"

    host = random_sequence(12_000, 44)
    microbes = [random_sequence(9_000, 101), random_sequence(8_000, 202)]
    write_fasta(reference_dir / "GRCh38.p14.fa", "synthetic_GRCh38_p14", host)
    write_fasta(data_dir / "microbial_fixture.fa", "synthetic_microbe", microbes[0])
    for assembler, sequence in (("megahit", microbes[0]), ("spades", microbes[1])):
        for index in (1, 2):
            write_fasta(
                data_dir / "bins" / assembler / f"bin_{index:03d}.fa",
                f"{assembler}_bin_{index:03d}",
                sequence[:2400],
            )

    index_prefix = reference_dir / "GRCh38_p14"
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        index_path = Path(f"{index_prefix}{suffix}")
        if not index_path.exists():
            index_path.write_text("stub index\n", encoding="utf-8")

    directories = (
        data_dir / "databases" / "checkm2",
        data_dir / "databases" / "gtdbtk",
        data_dir / "databases" / "phylophlan",
        data_dir / "databases" / "eggnog",
        data_dir / "databases" / "interproscan",
        data_dir / "databases" / "gunc",
        data_dir / "licenses",
        data_dir / "software" / "genemark",
    )
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)

    (data_dir / "databases" / "gunc" / "gunc_db.dmnd").write_text(
        "stub database\n", encoding="utf-8"
    )
    (data_dir / "databases" / "checkm2" / "uniref100.KO.1.dmnd").write_text(
        "stub database\n", encoding="utf-8"
    )
    (data_dir / "licenses" / "gm_key").write_text("STUB_KEY\n", encoding="utf-8")

    return host, microbes


def main() -> int:
    args = parse_args()
    project_dir = args.project_dir.resolve()
    data_dir = project_dir / "tests" / "generated_data"
    host, microbes = create_stub_resources(project_dir)

    rows: list[dict[str, str]] = []
    for index, sample in enumerate(("sample_A", "sample_B")):
        first = data_dir / f"{sample}_R1.fastq.gz"
        second = data_dir / f"{sample}_R2.fastq.gz"
        write_reads(first, second, sample, host, microbes[index], args.pairs)
        rows.append(
            {
                "sample": sample,
                "fastq_1": str(first.resolve()),
                "fastq_2": str(second.resolve()),
            }
        )

    with (data_dir / "samplesheet.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=("sample", "fastq_1", "fastq_2"))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Generated synthetic resources under {data_dir.parent}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
