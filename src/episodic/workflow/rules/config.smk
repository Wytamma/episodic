# convert config to yaml
import yaml

rule write_config_to_file:
    output: OUT_DIR / "config.yaml"
    params: yaml_config=yaml.dump(config)
    shell:
        """
        echo "{params.yaml_config}" > {output}
        """