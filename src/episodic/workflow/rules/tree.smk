
rule max_clade_credibility_tree:
    """
    Makes trace plots from the beast log file.
    """
    input:
        rules.run_beast.output.beast_trees_file,
    output:
        CLOCK_DIR / "{clock}" / "{name}" / "{name}.mcc.{heights}.nexus",
    params:
        burnin = int(int(config['beast']['samples']) * 0.1),
    conda:
        "../envs/beast.yml"
    shell:
        """
        treeannotator -burninTrees {params.burnin} -heights {wildcards.heights} {input} {output}
        """


rule max_clade_credibility_tree_newick:
    """
    Makes trace plots from the beast log file.
    """
    input:
        rules.max_clade_credibility_tree.output,
    output:
        CLOCK_DIR / "{clock}" / "{name}" / "{name}.mcc.{heights}.nwk",
    conda:
        "../envs/phylo.yml"
    shell:
        "${{CONDA_PREFIX}}/bin/python {SCRIPT_DIR}/tree_converter.py {input} {output} --node-label posterior"


rule max_clade_credibility_tree_render:
    """
    Renders the MCC tree in SVG format.
    """
    input:
        rules.max_clade_credibility_tree.output,
    output:
        CLOCK_DIR / "{clock}" / "{name}" / "{name}.mcc.{heights}.svg",
        CLOCK_DIR / "{clock}" / "{name}" / "{name}.mcc.{heights}.height_0.95_HPD.svg",
        CLOCK_DIR / "{clock}" / "{name}" / "{name}.mcc.{heights}.posterior.svg"
    params:
        mrsd = most_recent_sampling_date,
        prefix = lambda wildcards: f"{CLOCK_DIR}/{wildcards.clock}/{wildcards.name}/{wildcards.name}.mcc.{wildcards.heights}",
    conda:
        "../envs/ggtree.yml"
    shell:
        "${{CONDA_PREFIX}}/bin/Rscript {SCRIPT_DIR}/plot_mcc_tree.R --input {input} --output-prefix {params.prefix} --mrsd {params.mrsd}"


rule rate_quantile_analysis:
    input:
        rules.run_beast.output.beast_trees_file,
    output:
        csv = CLOCK_DIR / "{clock}" / "{name}" / "{name}.stem.rate_quantiles.csv",
        svg = CLOCK_DIR / "{clock}" / "{name}" / "{name}.stem.rate_quantiles.svg",
    params:
        groups = " ".join(f"-g {group}" for group in config['group']),
    conda:
        "../envs/phylo.yml"
    shell:
        """
        python {SCRIPT_DIR}/phylo_rate_quantile_analysis.py \
          {input} \
          --output-csv {output.csv} \
          --output-plot {output.svg} \
          {params.groups}
        """

