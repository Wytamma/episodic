import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable, Tuple


def extract_fasta_from_beast2_xml(xml_path: Path) -> Iterable[Tuple[str, str]]:
    """Extract (taxon_id, sequence) pairs from a BEAST2 XML alignment."""
    root = ET.parse(xml_path).getroot()

    for sequence_elem in root.findall(".//{*}sequence"):
        taxon_elem = sequence_elem.find("{*}taxon")
        if taxon_elem is None:
            continue

        taxon_id = taxon_elem.attrib.get("idref")
        sequence = (taxon_elem.tail or "").strip()
        if not taxon_id or not sequence:
            continue

        yield taxon_id, sequence


def write_fasta(records: Iterable[Tuple[str, str]], output_path: Path) -> int:
    """Write sequence records to a FASTA file and return record count."""
    count = 0
    with output_path.open("w") as handle:
        for taxon_id, sequence in records:
            handle.write(f">{taxon_id}\n{sequence}\n")
            count += 1
    return count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract FASTA records from BEAST2 XML sequences."
    )
    parser.add_argument("xml", type=Path, help="Input BEAST2 XML file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output FASTA path (default: <xml_stem>.fasta)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output = args.output or args.xml.with_suffix(".fasta")
    count = write_fasta(extract_fasta_from_beast2_xml(args.xml), output)
    print(f"Wrote {count} sequences to {output}")


if __name__ == "__main__":
    main()
