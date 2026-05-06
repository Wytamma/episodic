from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

Range = tuple[int, int]
FastaRecords = dict[str, str]

# Coordinates embedded from tests/data/genome_annotations.json
BASE_PARTITIONS: dict[str, list[Range]] = {
    "ORF1a": [(266, 13468)],
    "ORF1b": [(13468, 21555)],
    "S": [(21563, 25384)],
    "N": [(28274, 29533)],
}

ACCESSORY_PARTITIONS: dict[str, list[Range]] = {
    "ORF3a": [(25393, 26220)],
    "E": [(26245, 26472)],
    "M": [(26523, 27191)],
    "ORF6": [(27202, 27387)],
    "ORF7a": [(27394, 27759)],
    "ORF7b": [(27756, 27887)],
    "ORF8": [(27894, 28259)],
}


def build_partitions(*, concat_orf3_orf8: bool) -> dict[str, list[Range]]:
    partitions = dict(BASE_PARTITIONS)
    if concat_orf3_orf8:
        partitions["ORF3-ORF8"] = [
            ACCESSORY_PARTITIONS["ORF3a"][0],
            ACCESSORY_PARTITIONS["E"][0],
            ACCESSORY_PARTITIONS["M"][0],
            ACCESSORY_PARTITIONS["ORF6"][0],
            ACCESSORY_PARTITIONS["ORF7a"][0],
            ACCESSORY_PARTITIONS["ORF7b"][0],
            ACCESSORY_PARTITIONS["ORF8"][0],
        ]
    else:
        partitions.update(ACCESSORY_PARTITIONS)
    return partitions


def read_fasta(path: str | Path) -> FastaRecords:
    records: FastaRecords = {}
    header: str | None = None
    seq_lines: list[str] = []

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    if not seq_lines:
                        msg = f"Empty sequence for {header}."
                        raise ValueError(msg)
                    records[header] = "".join(seq_lines)

                header = line[1:].strip()
                seq_lines = []
                continue

            if header is None:
                msg = f"Found sequence content before first FASTA header at line {line_number}."
                raise ValueError(msg)
            seq_lines.append(line)

    if header is not None:
        if not seq_lines:
            msg = f"Empty sequence for {header}."
            raise ValueError(msg)
        records[header] = "".join(seq_lines)

    if not records:
        msg = "No FASTA records found."
        raise ValueError(msg)

    return records


def slice_sequence_by_ranges(sequence: str, ranges: Iterable[Range]) -> str:
    fragments: list[str] = []
    for start, end in ranges:
        if start < 1 or end < start:
            msg = f"Invalid range {start}-{end}."
            raise ValueError(msg)
        if end > len(sequence):
            msg = f"Range {start}-{end} exceeds sequence length {len(sequence)}."
            raise ValueError(msg)
        fragments.append(sequence[start - 1 : end])
    return "".join(fragments)


def write_partition_fastas(
    input_fasta: str | Path,
    output_dir: str | Path,
    partitions: dict[str, list[Range]],
) -> list[Path]:
    records = read_fasta(input_fasta)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    written_files: list[Path] = []
    for partition_name, ranges in partitions.items():
        output_file = output_path / f"{partition_name}.fasta"
        with output_file.open("w", encoding="utf-8") as handle:
            for taxon, full_sequence in records.items():
                partition_sequence = slice_sequence_by_ranges(full_sequence, ranges)
                handle.write(f">{taxon}\n{partition_sequence}\n")
        written_files.append(output_file)

    return written_files


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract SARS-CoV-2 partitions from an aligned FASTA using embedded "
            "coordinates (ORF1a, ORF1b, ORF3-ORF8 plus E and M, S, N)."
        )
    )
    parser.add_argument("input_fasta", type=Path, help="Path to input aligned FASTA")
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for partition FASTA outputs (default: current directory)",
    )
    parser.add_argument(
        "--concat-ORF3-E-M-ORF8",
        dest="concat_orf3_e_m_orf8",
        action="store_true",
        help="Concatenate ORF3a/E/M/ORF6/ORF7a/ORF7b/ORF8 into a single ORF3-ORF8 plus E and M partition.",
    )

    args = parser.parse_args()
    partitions = build_partitions(concat_orf3_orf8=args.concat_orf3_e_m_orf8)
    written = write_partition_fastas(args.input_fasta, args.output_dir, partitions)

    for path in written:
        sys.stdout.write(f"{path}\n")


if __name__ == "__main__":
    main()
