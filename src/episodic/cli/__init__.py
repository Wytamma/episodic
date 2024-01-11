# SPDX-FileCopyrightText: 2023-present Wytamma Wirth <wytamma.wirth@me.com>
#
# SPDX-License-Identifier: MIT
from pathlib import Path

from snk.cli import CLI

episodic = CLI(pipeline_dir_path = Path(__file__).parent.parent)

@episodic.app.command()
def template():
    """Display the BEAST XML template."""
    template_path = episodic.pipeline.path / "workflow/templates/beast_xml_template.jinja"
    template = template_path.read_text()
    print(template)
