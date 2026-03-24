# convert config to yaml
from pathlib import Path

import yaml


def _to_plain_types(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _to_plain_types(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_to_plain_types(item) for item in value]
    if isinstance(value, tuple):
        return [_to_plain_types(item) for item in value]
    return value

rule write_config_to_file:
    output: OUT_DIR / "config.yaml"
    params: yaml_config=yaml.safe_dump(_to_plain_types(config), sort_keys=False)
    shell:
        """
        echo "{params.yaml_config}" > {output}
        """