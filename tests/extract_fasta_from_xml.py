import xml.etree.ElementTree as ET
from pathlib import Path


def extract_fasta_from_beast2_xml(xml: Path):
    """Extract fasta from BEAST2 xml file."""
    with xml.open() as f:
        xml_str = f.read()
    xml_str = xml_str.replace('xmlns="http://www.beast2.org/2.0"', "")
    root = ET.fromstring(xml_str)
    for elem in root.iter():
        if elem.tag.endswith("sequence"):
            taxon = elem.find("taxon")
            yield taxon.attrib["idref"], taxon.tail

if __name__ == "__main__":
    xml = Path("/Users/wwirth/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/MDU/programming/stem-rate/BA.2.86/FLC_MLE_only/variants1208Aligned_FLC_getMLE.xml")
    for seq_id, seq in extract_fasta_from_beast2_xml(xml):
        print(f">{seq_id}\n{seq.strip()}")
