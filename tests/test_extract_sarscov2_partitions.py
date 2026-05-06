from episodic.workflow.scripts.extract_sarscov2_partitions import (
    ACCESSORY_PARTITIONS,
    build_partitions,
    slice_sequence_by_ranges,
)


def test_build_partitions_includes_e_and_m_when_not_concatenated():
    partitions = build_partitions(concat_orf3_orf8=False)

    assert partitions["E"] == ACCESSORY_PARTITIONS["E"]
    assert partitions["M"] == ACCESSORY_PARTITIONS["M"]
    assert "ORF3-ORF8" not in partitions


def test_build_partitions_concatenates_orf3_to_orf8_with_e_and_m():
    partitions = build_partitions(concat_orf3_orf8=True)

    assert partitions["ORF3-ORF8"] == [
        ACCESSORY_PARTITIONS["ORF3a"][0],
        ACCESSORY_PARTITIONS["E"][0],
        ACCESSORY_PARTITIONS["M"][0],
        ACCESSORY_PARTITIONS["ORF6"][0],
        ACCESSORY_PARTITIONS["ORF7a"][0],
        ACCESSORY_PARTITIONS["ORF7b"][0],
        ACCESSORY_PARTITIONS["ORF8"][0],
    ]


def test_concatenated_partition_slices_all_accessory_regions():
    partitions = build_partitions(concat_orf3_orf8=True)
    sequence = "A" * 29903

    partition_sequence = slice_sequence_by_ranges(sequence, partitions["ORF3-ORF8"])

    expected_length = sum(end - start + 1 for start, end in partitions["ORF3-ORF8"])
    assert len(partition_sequence) == expected_length
