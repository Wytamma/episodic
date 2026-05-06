
import hashlib

MAX_BEAST_SEED = 2**31 - 1

def beast_seed(wildcards):
    """Generate a deterministic seed for BEAST based on the base seed, clock, and name.
    This ensures that the same seed is used for the same clock and name across runs, while
    allowing for different seeds for different clocks and names.
    """
    base_seed = config["beast"].get("seed", "episodic")
    seed_key = f"{base_seed}:{wildcards.clock}:{wildcards.name}"
    seed_hash = hashlib.sha256(seed_key.encode()).hexdigest()
    return int(seed_hash[:8], 16) % MAX_BEAST_SEED + 1

rule create_beast_xml:
    input:
        alignments = config["alignment"],
        groups_file = OUT_DIR / "taxon_groups.tsv",
    output:
        beast_XML_file = CLOCK_DIR / "{clock}" / "{name}" / "{name}.xml",
    params:
        template = beast_xml_template,
        alignments = lambda wildcards, input: " ".join(f"--alignment {alignment}" for alignment in input.alignments),
        date_delimiter = "\|" if config.get("date_delimiter") == "|" else config.get("date_delimiter"),
        date_index = config.get("date_index", -1),
        clock = lambda wildcards: wildcards.clock.split("_")[0],
        rate_gamma_prior_shape = config.get("rate_gamma_prior_shape"),
        rate_gamma_prior_scale = config.get("rate_gamma_prior_scale"),
        chain_length = config["beast"].get("chain_length"),
        samples = config["beast"].get("samples"),
        mle = lambda wildcards: "--mle" if "mle" in wildcards.name else "",
        mle_chain_length = f"--mle-chain-length {config['marginal_likelihood'].get('chain_length')}",
        mle_path_steps = f"--mle-path-steps {config['marginal_likelihood'].get('path_steps')}",
        mle_log_every = f"--mle-log-every {config['marginal_likelihood'].get('log_every')}",
        no_trace = lambda wildcards: "--no-trace" if "mle" in wildcards.name else "",
        no_trees = lambda wildcards: "--no-trees" if "mle" in wildcards.name else "",
        fixed_tree = f'--fixed-tree {config.get("newick")}'  if config.get("newick") else "",
        foreground_label = f'--foreground-label {config.get("foreground_label")}' if config.get("foreground_label") else "",
        background_label = f'--background-label {config.get("background_label")}' if config.get("background_label") else "",
    shell:
        """
        python {SCRIPT_DIR}/populate_beast_template.py \
            {params.template} \
            --output {output.beast_XML_file} \
            {params.alignments} \
            --date-delimiter {params.date_delimiter} \
            --date-index {params.date_index} \
            --groups-file {input.groups_file} \
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
            {params.fixed_tree} \
            {params.foreground_label} \
            {params.background_label}
        """

TREES = {"beast_trees_file": CLOCK_DIR / "{clock}" / "{name}" / "{name}.trees"} if config.get("trees") else {}

rule beast:
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
    params:
        beast_args = lambda wildcards, resources: (
            config["beast"].get("args", "").replace("-beagle_CPU", "-beagle_GPU")
            if "gpu" in str(getattr(resources, "gres", ""))
            else config["beast"].get("args", "")
        ),
        seed = beast_seed,
    shell:
        """
        beast -seed {params.seed} -working -overwrite {params.beast_args} -threads {threads} {input.beast_XML_file} > {output.beast_stdout_file}
        """

MLE_OUT_DIR = OUT_DIR / "mle" / "{clock}"

use rule create_beast_xml as create_mle_xml with:
    output:
        beast_XML_file = MLE_OUT_DIR / "{name}.xml",

use rule beast as mle with:
    input:
        beast_XML_file = rules.create_mle_xml.output.beast_XML_file,
    output:
        beast_stdout_file = MLE_OUT_DIR / "{name}.stdout",
