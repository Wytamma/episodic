# SPDX-FileCopyrightText: 2023-present Wytamma Wirth <wytamma.wirth@me.com>
#
# SPDX-License-Identifier: MIT
from pathlib import Path

from snk.cli import CLI

episodic = CLI(pipeline_dir_path = Path(__file__).parent.parent)
