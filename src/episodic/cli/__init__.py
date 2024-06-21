# SPDX-FileCopyrightText: 2023-present Wytamma Wirth <wytamma.wirth@me.com>
#
# SPDX-License-Identifier: MIT
from pathlib import Path

from snk_cli import CLI

episodic = CLI(Path(__file__).parent.parent)

@episodic.app.command()
def template():
    """Show the BEAST XML template."""
    template_path = episodic.workflow.path / "workflow/templates/beast_xml_template.jinja"
    template = template_path.read_text()
    episodic.echo(template)
