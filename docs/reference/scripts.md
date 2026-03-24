## Run scripts outside the pipeline

You can run workflow helper scripts directly via the CLI without running the full Snakemake pipeline.

General form:

```console
episodic script run <script_name> --env <env_name> -- <script_args>
```

Example (show script help):

```console
episodic script run plot_partition_local_rate_posteriors --env python -- --help
```

Notes:

- `--env` selects the workflow conda environment declared for that script (for example `python`, `phylo`, or `ggtree`).
- The extra `--` separates `episodic script run` options from the script's own arguments.

::: src.episodic.workflow.utils

::: src.episodic.workflow.scripts.extract_mle

::: src.episodic.workflow.scripts.phylo_rate_quantile_analysis

::: src.episodic.workflow.scripts.arviz_output

::: src.episodic.workflow.scripts.calculate_odds

::: src.episodic.workflow.scripts.plot_partition_local_rate_posteriors

::: src.episodic.workflow.scripts.plot_traces

::: src.episodic.workflow.scripts.tree_converter

::: src.episodic.workflow.scripts.populate_beast_template

