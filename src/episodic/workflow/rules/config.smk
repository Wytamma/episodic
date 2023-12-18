# convert config to yaml
import yaml

rule write_config_to_file:
    output: OUT_DIR / "config.yaml"
    params: yaml_config=yaml.dump(config)
    resources: **config["resources"].get("default", {}),
    shell:
        """
        echo "{params.yaml_config}" > {output}
        """