import argparse
import logging
import re

from pathlib import Path
from typing import Dict, Iterable, List, Tuple


LOGGER = logging.getLogger(__name__)
EXPECTED_MATRIX_PARTS = 2

Range = Tuple[int, int]
Partitions = Dict[str, List[Range]]
FastaRecords = Dict[str, str]


def _remove_nexus_comments(nexus_text: str) -> str:
    """Remove square-bracket NEXUS comments."""
    return re.sub(r"\[.*?\]", "", nexus_text, flags=re.DOTALL)


def _parse_range_token(token: str) -> Range:
    """
    Parse a NEXUS charset range token.

    Supported examples:
      - '151-920'
      - '1-150'
      - '921-1695'
      - '42'   -> interpreted as (42, 42)
      - '1-300\\3' or '1-300/3' -> interpreted as (1, 300); step is ignored
    """
    token = token.strip()
    if not token:
        msg = "Empty range token"
        raise ValueError(msg)

    token = re.split(r"[\\/]", token, maxsplit=1)[0].strip()

    if "-" in token:
        start_str, end_str = token.split("-", 1)
        start = int(start_str)
        end = int(end_str)
    else:
        start = end = int(token)

    if start < 1:
        msg = f"Invalid range start: {token}"
        raise ValueError(msg)
    if start > end:
        msg = f"Invalid range: {token}"
        raise ValueError(msg)

    return start, end


def extract_nexus_partitions(nexus_text: str) -> Partitions:
    """
    Extract charset partitions from a NEXUS string.

    Example matched lines:
        charset HA1 = 151-920;
        charset HA2 = 1-150,921-1695;

    Returns:
        {
            "HA1": [(151, 920)],
            "HA2": [(1, 150), (921, 1695)]
        }
    """
    nexus_text = _remove_nexus_comments(nexus_text)

    pattern = re.compile(
        r"charset\s+([A-Za-z0-9_.-]+)\s*=\s*([^;]+);",
        flags=re.IGNORECASE | re.DOTALL,
    )

    partitions: Partitions = {}

    for match in pattern.finditer(nexus_text):
        name = match.group(1)
        coords_text = match.group(2)

        ranges: List[Range] = []
        for range_token in coords_text.split(","):
            cleaned_token = range_token.strip()
            if cleaned_token:
                ranges.append(_parse_range_token(cleaned_token))

        partitions[name] = ranges

    return partitions


def _extract_matrix_block(nexus_text: str) -> str:
    """Extract the MATRIX block from a NEXUS alignment."""
    nexus_text = _remove_nexus_comments(nexus_text)
    match = re.search(r"\bmatrix\b(.*?)\s*;", nexus_text, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        msg = "Could not find a MATRIX block in the NEXUS file."
        raise ValueError(msg)
    return match.group(1)


def extract_alignment_from_nexus(nexus_text: str) -> FastaRecords:
    """
    Extract taxon sequences from a NEXUS MATRIX block.

    Supports repeated taxon rows by concatenating sequence fragments.
    """
    matrix_text = _extract_matrix_block(nexus_text)
    records: FastaRecords = {}

    for raw_line in matrix_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        parts = re.split(r"\s+", line, maxsplit=1)
        if len(parts) != EXPECTED_MATRIX_PARTS:
            continue

        taxon, sequence = parts
        sequence = re.sub(r"\s+", "", sequence)
        if not sequence:
            continue

        records[taxon] = records.get(taxon, "") + sequence

    if not records:
        msg = "No sequences were found in the NEXUS MATRIX block."
        raise ValueError(msg)

    return records


def extract_nexus_partitions_from_file(path: str | Path) -> Partitions:
    """Read a NEXUS file and extract charset partitions."""
    text = Path(path).read_text(encoding="utf-8")
    return extract_nexus_partitions(text)


def extract_alignment_from_nexus_file(path: str | Path) -> FastaRecords:
    """Read a NEXUS file and extract the alignment matrix."""
    text = Path(path).read_text(encoding="utf-8")
    return extract_alignment_from_nexus(text)


def slice_sequence_by_ranges(sequence: str, ranges: Iterable[Range]) -> str:
    """Slice a sequence using 1-based inclusive NEXUS ranges."""
    fragments: List[str] = []
    for start, end in ranges:
        if end > len(sequence):
            msg = f"Partition range {start}-{end} exceeds sequence length {len(sequence)}."
            raise ValueError(msg)
        fragments.append(sequence[start - 1 : end])
    return "".join(fragments)


def write_partition_fastas(nexus_file: str | Path, output_dir: str | Path | None = None) -> List[Path]:
    """
    Extract charset partitions and write one FASTA file per partition.

    Returns the written FASTA paths.
    """
    nexus_path = Path(nexus_file)
    output_path = Path(output_dir) if output_dir is not None else nexus_path.parent
    output_path.mkdir(parents=True, exist_ok=True)

    nexus_text = nexus_path.read_text(encoding="utf-8")
    partitions = extract_nexus_partitions(nexus_text)
    alignment = extract_alignment_from_nexus(nexus_text)

    if not partitions:
        msg = "No charset partitions were found in the NEXUS file."
        raise ValueError(msg)

    written_files: List[Path] = []
    for partition_name, ranges in partitions.items():
        fasta_path = output_path / f"{partition_name}.fasta"
        with fasta_path.open("w", encoding="utf-8") as handle:
            for taxon, sequence in alignment.items():
                partition_sequence = slice_sequence_by_ranges(sequence, ranges)
                handle.write(f">{taxon}\n{partition_sequence}\n")
        written_files.append(fasta_path)

    return written_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract charset partitions from a NEXUS file.")
    parser.add_argument("nexus_file", type=Path, help="Path to the input NEXUS file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write partition FASTA files into. Defaults to the input file directory.",
    )
    args = parser.parse_args()

    partitions = extract_nexus_partitions_from_file(args.nexus_file)
    for name, ranges in partitions.items():
        LOGGER.info("%s: %s", name, ranges)

    written_files = write_partition_fastas(args.nexus_file, args.output_dir)
    for fasta_path in written_files:
        LOGGER.info("Wrote %s", fasta_path)
