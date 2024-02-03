
rule create_beast_xml:
    input:
        alignment = config["alignment"],
    output:
        beast_XML_file = CLOCK_DIR / "{clock}" / "{name}" / "{name}.xml",
    params:
        template = beast_xml_template,
        date_delimiter = "\|" if config.get("date_delimiter") == "|" else config.get("date_delimiter"),
        date_index = config.get("date_index", -1),
        groups = " ".join(config.get("group")),
        clock = lambda wildcards: wildcards.clock.split("_")[0],
        rate_gamma_prior_shape = config.get("rate_gamma_prior_shape"),
        rate_gamma_prior_scale = config.get("rate_gamma_prior_scale"),
        chain_length = config["beast"].get("chain_length"),
        samples = config["beast"].get("samples"),
        mle = lambda wildcards: "--mle" if "mle" in wildcards.name else "",
        mle_chain_length = f"--mle-chain-length {config['marginal_likelihood'].get('chain_length')}",
        mle_path_steps = f"--mle-path-steps {config['marginal_likelihood'].get('paths')}",
        mle_log_every = f"--mle-log-every {config['marginal_likelihood'].get('log_every')}",
        no_trace = lambda wildcards: "--no-trace" if "mle" in wildcards.name else "",
        no_trees = lambda wildcards: "--no-trees" if "mle" in wildcards.name else "",
        fixed_tree = f'--fixed-tree {config.get("newick")}'  if config.get("newick") else "",
    shell:
        """
        python {SCRIPT_DIR}/populate_beast_template.py \
            {params.template} \
            --output {output.beast_XML_file} \
            --alignment {input.alignment} \
            --date-delimiter {params.date_delimiter} \
            --date-index {params.date_index} \
            --groups {params.groups} \
            --clock {params.clock} \
            --rate-gamma-prior-shape {params.rate_gamma_prior_shape} \
            --rate-gamma-prior-scale {params.rate_gamma_prior_scale} \
            --chain-length {params.chain_length} \
            --samples {params.samples} \
            {params.mle} \
            {params.mle_chain_length} \
            {params.mle_path_steps} \
            {params.mle_log_every} \
            {params.no_trace} \
            {params.no_trees} \
            {params.fixed_tree}
        """

TREES = {"beast_trees_file": CLOCK_DIR / "{clock}" / "{name}" / "{name}.trees"} if config.get("trees") else {}

rule run_beast:
    input:
        beast_XML_file = rules.create_beast_xml.output.beast_XML_file,
    output:
        beast_stdout_file = CLOCK_DIR / "{clock}" / "{name}" / "{name}.stdout",
        beast_log_file = CLOCK_DIR / "{clock}" / "{name}" / "{name}.log",
        beast_trees_file =  CLOCK_DIR / "{clock}" / "{name}" / "{name}.trees",
    threads: config["beast"].get("threads")
    envmodules:
        *config["beast"].get("envmodules", []),
    conda:
        "../envs/beast.yml"
    shell:
        """
        beast -working -overwrite -beagle_GPU -threads {threads} {input.beast_XML_file} > {output.beast_stdout_file}
        """

MLE_OUT_DIR = OUT_DIR / "mle" / "{clock}"

use rule create_beast_xml as create_mle_xml with:
    output:
        beast_XML_file = MLE_OUT_DIR / "{name}.xml",

use rule run_beast as run_mle_beast with:
    input:
        beast_XML_file = rules.create_mle_xml.output.beast_XML_file,
    output:
        beast_stdout_file = MLE_OUT_DIR / "{name}.stdout",