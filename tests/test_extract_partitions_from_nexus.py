from episodic.workflow.scripts.extract_partitions_from_nexus import extract_alignment_from_nexus
from episodic.workflow.scripts.extract_partitions_from_nexus import extract_nexus_partitions
from episodic.workflow.scripts.extract_partitions_from_nexus import slice_sequence_by_ranges
from episodic.workflow.scripts.extract_partitions_from_nexus import write_partition_fastas


NEXUS_TEXT = """#NEXUS
begin data;
    dimensions ntax=2 nchar=8;
    format datatype=dna missing=? gap=-;
    matrix
    taxon1  ACGTACGT
    taxon2  TGCATGCA
    ;
end;

begin sets;
    charset part1 = 1-4;
    charset part2 = 5-8;
end;
"""


def test_extract_nexus_partitions():
    partitions = extract_nexus_partitions(NEXUS_TEXT)

    assert partitions == {
        "part1": [(1, 4)],
        "part2": [(5, 8)],
    }


def test_extract_alignment_from_nexus():
    alignment = extract_alignment_from_nexus(NEXUS_TEXT)

    assert alignment == {
        "taxon1": "ACGTACGT",
        "taxon2": "TGCATGCA",
    }


def test_slice_sequence_by_ranges():
    assert slice_sequence_by_ranges("ACGTACGT", [(1, 2), (5, 6)]) == "ACAC"


def test_write_partition_fastas(tmp_path):
    nexus_path = tmp_path / "example.nexus"
    nexus_path.write_text(NEXUS_TEXT, encoding="utf-8")

    written = write_partition_fastas(nexus_path, tmp_path)

    assert [path.name for path in written] == ["part1.fasta", "part2.fasta"]
    assert (tmp_path / "part1.fasta").read_text(encoding="utf-8") == (
        ">taxon1\nACGT\n>taxon2\nTGCA\n"
    )
    assert (tmp_path / "part2.fasta").read_text(encoding="utf-8") == (
        ">taxon1\nACGT\n>taxon2\nTGCA\n"
    )
