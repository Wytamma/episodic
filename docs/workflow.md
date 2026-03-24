# Workflow Outputs

This page documents all files produced by the Episodic workflow, including:

- Snakemake target outputs (requested by `rule all`)
- Additional side-effect files produced by plotting/report scripts

## Output root

All outputs are written under:

- `output.dir` (from config), or
- `output.dir/YYYY-MM-DDTHH:MM:SS` when `output.dated: true`

Examples below use `OUT_DIR` to represent that resolved output root.

## Naming tokens

- `{clock}`: clock model name including gamma prior suffix, e.g. `flc-stem_0.5_0.1`
- `{duplicate}`: BEAST run replicate index (1..N)
- `{heights}`: MCC summary option, e.g. `mean`, `median`, `keep`, `ca`

For BEAST run artifacts, `{name}` is typically `{clock}_{duplicate}`.

## Always produced

| File | Description |
|---|---|
| `OUT_DIR/config.yaml` | Rendered workflow configuration used for the run |
| `OUT_DIR/taxon_groups.tsv` | Taxon-to-group assignments used by local clock models |

## Core per-clock BEAST outputs (`beast.fit_clocks: true`)

For each `{clock}` and `{duplicate}`:

| File pattern | Description |
|---|---|
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.xml` | BEAST XML analysis file |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.stdout` | BEAST stdout log |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.log` | BEAST posterior trace log |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.trees` | Posterior tree samples |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}_trace_plots/` | Trace PNGs per variable |

## Per-clock summaries and rate plots

| File pattern | Description |
|---|---|
| `OUT_DIR/clocks/{clock}/{clock}-summary.csv` | ArviZ posterior summary table |
| `OUT_DIR/clocks/{clock}/{clock}-violin.svg` | Clock-rate violin plot |
| `OUT_DIR/clocks/clocks_{shape}_{scale}-violin.svg` | Combined violin plot across clocks |
| `OUT_DIR/clocks/clocks_{shape}_{scale}-trace.svg` | Combined trace plot across clocks |

### Additional side-effect files from rate plotting scripts

The plotting scripts also emit extra artifacts sharing the same output prefix:

- `*-violin-rug.svg`
- `*-violin-trimmed.svg`
- `*-violin-rug-trimmed.svg`
- `*-trace.svg`
- `*-trace.html`

These are generated for both per-clock and combined clock prefixes.

## MCC tree outputs

For each `{clock}`, `{duplicate}`, and configured `{heights}`:

| File pattern | Description |
|---|---|
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.mcc.{heights}.nexus` | MCC tree in NEXUS format |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.mcc.{heights}.nwk` | MCC tree in Newick format |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.mcc.{heights}.svg` | MCC tree render |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.mcc.{heights}.height_0.95_HPD.svg` | MCC tree with 95% HPD node ranges |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.mcc.{heights}.posterior.svg` | MCC tree with posterior node labels |

## Clock-type-specific outputs

### FLC odds / effect-size Bayes factors (`flc*` clocks)

| File pattern | Description |
|---|---|
| `OUT_DIR/clocks/{clock}/{clock}-odds.csv` | Per-clock odds/BF summary |

Additional side-effect per-duplicate files:

- `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}-odds.csv`

### Relaxed-clock stem quantile analysis (`relaxed*` clocks)

| File pattern | Description |
|---|---|
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.stem.rate_quantiles.csv` | Stem branch rank/quantile table |
| `OUT_DIR/clocks/{clock}/{clock}_{duplicate}/{clock}_{duplicate}.stem.rate_quantiles.svg` | Stem branch rank/quantile plot |

### FLC partition local-rate posterior comparison (`flc*` clocks)

| File pattern | Description |
|---|---|
| `OUT_DIR/clocks/{clock}/{clock}-partition_local_rate_posteriors.csv` | Long-format posterior rates used for plotting |
| `OUT_DIR/clocks/{clock}/{clock}-partition_local_rate_posteriors.svg` | Overlaid background vs local posterior densities by partition |

## Marginal likelihood outputs (`marginal_likelihood.estimate: true`)

Per MLE run files:

| File pattern | Description |
|---|---|
| `OUT_DIR/mle/{clock}/{clock}_mle_{duplicate}.xml` | MLE XML |
| `OUT_DIR/mle/{clock}/{clock}_mle_{duplicate}.stdout` | MLE stdout |

Aggregated MLE summaries:

| File | Description |
|---|---|
| `OUT_DIR/mle/mle.svg` | Marginal likelihood comparison plot |
| `OUT_DIR/mle/mle.csv` | Raw extracted MLE values |
| `OUT_DIR/mle/mle.grouped.csv` | Grouped/aggregated MLE values |

## Minimal directory example

```text
OUT_DIR/
	config.yaml
	taxon_groups.tsv
	clocks/
		flc-shared-stem-and-clade_0.5_0.1/
			flc-shared-stem-and-clade_0.5_0.1-summary.csv
			flc-shared-stem-and-clade_0.5_0.1-odds.csv
			flc-shared-stem-and-clade_0.5_0.1-partition_local_rate_posteriors.csv
			flc-shared-stem-and-clade_0.5_0.1-partition_local_rate_posteriors.svg
			flc-shared-stem-and-clade_0.5_0.1_1/
				flc-shared-stem-and-clade_0.5_0.1_1.log
				flc-shared-stem-and-clade_0.5_0.1_1.trees
				flc-shared-stem-and-clade_0.5_0.1_1.mcc.mean.nexus
				flc-shared-stem-and-clade_0.5_0.1_1.mcc.mean.nwk
				flc-shared-stem-and-clade_0.5_0.1_1.mcc.mean.svg
				flc-shared-stem-and-clade_0.5_0.1_1_trace_plots/
	mle/
		mle.svg
		mle.csv
		mle.grouped.csv
```

## Which outputs are present?

- `beast.fit_clocks: false` → suppresses most `clocks/` BEAST-derived outputs
- `marginal_likelihood.estimate: false` → no `mle/` outputs
- No `relaxed` clocks configured → no `*.stem.rate_quantiles.*`
- No `flc` clocks configured → no `*-partition_local_rate_posteriors.*`

