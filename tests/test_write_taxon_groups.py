from pathlib import Path

from episodic.workflow.scripts.write_taxon_groups import assign_group, read_group_members, write_taxon_groups


def test_assign_group_prefers_longest_match():
    assert assign_group("sample@B.1.1.7@2021.1", ["B.1", "B.1.1.7"]) == "B.1.1.7"


def test_write_taxon_groups_creates_mapping(tmp_path):
    output_path = tmp_path / "taxon_groups.tsv"

    write_taxon_groups(
        alignment_path=Path("tests/data/Tay2023.afa"),
        output_path=output_path,
        groups=["B.1.1.7", "B.1.351", "P.1", "B.1.617.2"],
    )

    contents = output_path.read_text()

    assert contents.startswith("taxon\tgroup\n")
    assert "EPI_ISL_3085251@B.1.1.7@Italy@2021.11\tB.1.1.7" in contents


def test_read_group_members_ignores_unassigned_taxa(tmp_path):
    groups_path = tmp_path / "taxon_groups.tsv"
    groups_path.write_text("taxon\tgroup\nfoo\tA\nbar\t\n")

    assert read_group_members(groups_path) == {"A": ["foo"]}