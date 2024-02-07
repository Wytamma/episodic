MLE_DIR = OUT_DIR / "mle"

rule extract_mle:
    input:
        expand(OUT_DIR / "mle" / "{clock}" / "{clock}_mle_{duplicate}.stdout", clock=clocks, duplicate=mle_duplicates),
    output:
        MLE_DIR / "mle.svg",
        MLE_DIR / "mle.csv",
        MLE_DIR / "mle.grouped.csv",
    conda:
        "../envs/python.yml"
    shell:
        """
        python {SCRIPT_DIR}/extract_mle.py {MLE_DIR}
        """

rule plot_traces:
    """
    Makes trace plots from the beast log file.
    """
    input:
        rules.run_beast.output.beast_log_file,
    output:
        directory(CLOCK_DIR / "{clock}" / "{name}" / "{name}_trace_plots/"),
    conda:
        "../envs/plot_traces.yml"
    shell:
        """
        ${{CONDA_PREFIX}}/bin/python {SCRIPT_DIR}/plot_traces.py {input} {output}
        """

rule summary:
    """
    Makes combined summarys from the beast log file.
    """
    input:
        PER_CLOCK_LOG_FILES,
    output:
        posterior_svg=CLOCK_DIR / "{clock}" / "{clock}-summary.csv",
    params:
        output=lambda wildcards: CLOCK_DIR / f"{wildcards.clock}" / f"{wildcards.clock}-summary.csv",
    conda:
        "../envs/arviz.yml"
    shell:
        """
        ${{CONDA_PREFIX}}/bin/python {SCRIPT_DIR}/arviz_output.py summary {input} {params.output}
        """

rule plot_flc_rates:
    """
    Makes plots from the flc clock files.
    """
    input:
        lambda wildcards: [CLOCK_DIR / wildcards.clock / f"{wildcards.clock}_{duplicate}" / f"{wildcards.clock}_{duplicate}.log" for duplicate in duplicates if "flc" in wildcards.clock],
    output:
        rate_svg=CLOCK_DIR / "{clock}" / "{clock}-violin.svg" 
    params:
        output_prefix=lambda wildcards: CLOCK_DIR / f"{wildcards.clock}" / f"{wildcards.clock}",
        gamma_shape=rate_gamma_prior_shape,
        gamma_scale=rate_gamma_prior_scale,
    conda:
        "../envs/arviz.yml"
    shell:
        """
        ${{CONDA_PREFIX}}/bin/python {SCRIPT_DIR}/arviz_output.py rates {input} --output-prefix {params.output_prefix} --gamma-shape {params.gamma_shape} --gamma-scale {params.gamma_scale}
        """

use rule plot_flc_rates as plot_rates with:
    input:
        expand(CLOCK_DIR / "{clock}" / "{clock}_{duplicate}" / "{clock}_{duplicate}.log", clock=clocks, duplicate=duplicates),
    output:
        clocks_violin=CLOCK_DIR / "clocks_{rate_gamma_prior_shape}_{rate_gamma_prior_scale}-violin.svg",
        clocks_trace=CLOCK_DIR / "clocks_{rate_gamma_prior_shape}_{rate_gamma_prior_scale}-trace.svg",
    params:
        output_prefix=CLOCK_DIR / f"clocks_{rate_gamma_prior_shape}_{rate_gamma_prior_scale}",
        gamma_shape=rate_gamma_prior_shape,
        gamma_scale=rate_gamma_prior_scale,

rule calculate_odds:
    input:
        lambda wildcards: [CLOCK_DIR / wildcards.clock / f"{wildcards.clock}_{duplicate}" / f"{wildcards.clock}_{duplicate}.log" for duplicate in duplicates if "flc" in wildcards.clock],
    output:
        rate_svg=CLOCK_DIR / "{clock}" / "{clock}-odds.csv" 
    params:
        directory=lambda wildcards: CLOCK_DIR / f"{wildcards.clock}",
        gamma_shape=rate_gamma_prior_shape,
        gamma_scale=rate_gamma_prior_scale,
    conda:
        "../envs/python.yml"
    shell:
        """
        ${{CONDA_PREFIX}}/bin/python {SCRIPT_DIR}/calculate_odds.py {input} {output} --gamma-shape {params.gamma_shape} --gamma-scale {params.gamma_scale}
        for file in {input}
        do
            ${{CONDA_PREFIX}}/bin/python {SCRIPT_DIR}/calculate_odds.py $file ${{file%.log}}-odds.csv --gamma-shape {params.gamma_shape} --gamma-scale {params.gamma_scale}
        done
        """