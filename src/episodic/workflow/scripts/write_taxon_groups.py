import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import typer


def read_fasta_headers(alignment_path: Path) -> List[str]:
    """Return FASTA headers without the leading '>' character."""
    headers: List[str] = []
    with alignment_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line.strip()[1:])
    return headers


def assign_group(taxon: str, groups: Iterable[str]) -> str:
    """Assign the longest matching configured group to a taxon label."""
    matches = [group for group in groups if group in taxon]
    if not matches:
        return ""
    return max(matches, key=len)


def build_group_members(taxa: Iterable[str], groups: Iterable[str]) -> Dict[str, List[str]]:
    """Map each configured group to the taxa assigned to it."""
    group_members: Dict[str, List[str]] = {group: [] for group in groups}
    for taxon in taxa:
        group = assign_group(taxon, groups)
        if group:
            group_members[group].append(taxon)
    return group_members


def read_group_members(groups_path: Path) -> Dict[str, List[str]]:
    """Read a TSV mapping taxa to groups and return members grouped by label."""
    group_members: Dict[str, List[str]] = {}

    with groups_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxon = row["taxon"].strip()
            group = row["group"].strip()
            if not group:
                continue
            group_members.setdefault(group, []).append(taxon)

    return group_members


def write_taxon_groups(alignment_path: Path, output_path: Path, groups: List[str]) -> None:
    """Write a TSV mapping taxa to configured groups."""
    headers = read_fasta_headers(alignment_path)
    group_members = build_group_members(headers, groups)
    assigned_groups = {taxon: group for group, taxa in group_members.items() for taxon in taxa}
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as handle:
        handle.write("taxon\tgroup\n")
        for taxon in headers:
            handle.write(f"{taxon}\t{assigned_groups.get(taxon, '')}\n")


def main(
    alignment_path: Path = typer.Argument(..., help="Alignment FASTA used in the analysis."),
    output_path: Path = typer.Argument(..., help="Output TSV path for taxon group assignments."),
    groups: Optional[List[str]] = typer.Option(None, "--group", help="Configured group labels."),
) -> None:
    """Create a taxon-to-group mapping file for downstream tree annotation."""
    write_taxon_groups(alignment_path=alignment_path, output_path=output_path, groups=groups or [])


if __name__ == "__main__":
    typer.run(main)
