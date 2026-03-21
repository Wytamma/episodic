from pathlib import Path

from episodic.workflow.scripts.populate_beast_template import populate_beast_template


def test_populate_beast_template_supports_multiple_partitions(tmp_path):
    template_path = Path("src/episodic/workflow/templates/beast_xml_template.jinja")
    alignment_path = Path("tests/data/BA.2.86.afa")

    xml = populate_beast_template(
        work_dir=tmp_path,
        name="partitioned",
        template_path=template_path,
        alignment_paths=[alignment_path, alignment_path],
        groups=["BA.2.86"],
        clock="flc-stem",
        date_delimiter="@",
        date_index=-1,
    )

    assert xml.count("<alignment id=") == 2
    assert "alignment.BA.2.86\"" in xml
    assert "alignment.BA.2.86_2\"" in xml
    assert "BA.2.86.clock.rate" in xml
    assert "BA.2.86_2.clock.rate" in xml
    assert "BA.2.86.BA.2.86.stem.rate" in xml
    assert "BA.2.86_2.BA.2.86.stem.rate" in xml


def test_populate_beast_template_supports_stem_and_clade_models(tmp_path):
    template_path = Path("src/episodic/workflow/templates/beast_xml_template.jinja")
    alignment_path = Path("tests/data/BA.2.86.afa")

    xml = populate_beast_template(
        work_dir=tmp_path,
        name="partitioned",
        template_path=template_path,
        alignment_paths=[alignment_path, alignment_path],
        groups=["BA.2.86"],
        clock="flc-stem-and-clade",
        date_delimiter="@",
        date_index=-1,
    )

    assert "BA.2.86.BA.2.86.stem_and_clade.rate" in xml
    assert "BA.2.86_2.BA.2.86.stem_and_clade.rate" in xml
    assert "BA.2.86.BA.2.86.stem.rate" not in xml
    assert "BA.2.86.BA.2.86.clade.rate" not in xml


def test_populate_beast_template_supports_shared_stem_and_clade_models(tmp_path):
    template_path = Path("src/episodic/workflow/templates/beast_xml_template.jinja")
    alignment_path = Path("tests/data/BA.2.86.afa")

    xml = populate_beast_template(
        work_dir=tmp_path,
        name="partitioned",
        template_path=template_path,
        alignment_paths=[alignment_path, alignment_path],
        groups=["BA.2.86"],
        clock="flc-shared-stem-and-clade",
        date_delimiter="@",
        date_index=-1,
    )

    assert "BA.2.86.stem_and_clade.rate" in xml
    assert "BA.2.86_2.stem_and_clade.rate" in xml
    assert "BA.2.86.stem.rate" not in xml
    assert "BA.2.86.clade.rate" not in xml
