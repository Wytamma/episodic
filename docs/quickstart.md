# Quickstart

Before running Episodic, we recommend reading Tay et al. ([2022](https://pubmed.ncbi.nlm.nih.gov/35038741/) and [2023](https://academic.oup.com/mbe/article/40/10/msad212/7280106)) and becoming familiar with Fixed Local Clock (FLC) models.

## Installation

Install Episodic with pip:

```console
pip install episodic
```

Episodic uses Snakemake and conda environments for the workflow tools. After installing the package, create the workflow environments:

```console
episodic env create
```

See [Installation](installation.md) for more detail.

## CLI overview

The main entry point is:

```console
episodic run [OPTIONS]
```

For a minimal FLC analysis, provide an alignment and one or more foreground groups:

```console
episodic run \
  --alignment data.fasta \
  --group group1 \
  --group group2
```

List-valued options can be repeated. For example, multiple alignment partitions are supplied by repeating `--alignment`, and multiple clock models are supplied by repeating `--clock`:

```console
episodic run \
  --alignment HA1.fasta \
  --alignment HA2.fasta \
  --group bird \
  --clock strict \
  --clock flc-stem
```

Episodic accepts both workflow options and Snakemake options. Any unrecognised arguments are passed through to Snakemake:

```console
episodic run \
  --config config.yaml \
  --profile slurm \
  --set-resources beast:mem_mb=16G \
  --set-resources beast:runtime=04:00:00
```

Useful general options:

| CLI option | Description |
| --- | --- |
| `--config FILE` | Load a YAML workflow config. CLI values override matching config values. |
| `--profile`, `-p` | Use a Snakemake profile, for example `slurm` or `slurm-gpu`. |
| `--dry`, `-n` | Dry-run the workflow and print what would run. |
| `--dag`, `-d` | Write the workflow DAG to `.pdf`, `.png`, or `.svg`. |
| `--cores`, `-c` | Number of local cores Snakemake may use. |
| `--no-conda` | Disable workflow conda environments. |
| `--verbose`, `-v` | Run with verbose output. |
| `--help-snakemake`, `-hs` | Show Snakemake help. |

Other top-level commands include:

| Command | Description |
| --- | --- |
| `episodic info` | Show workflow information. |
| `episodic config` | Print the workflow configuration schema/defaults. |
| `episodic template` | Print the bundled BEAST XML template. |
| `episodic env` | Manage workflow conda environments. |
| `episodic script` | Run workflow helper scripts directly. |
| `episodic profile` | Show bundled Snakemake profiles. |

## Configuration files

You can run the workflow entirely from command-line options, but a YAML config is usually easier to reproduce:

```yaml
alignment:
  - data.fasta
group:
  - group1
  - group2
clock:
  - flc-stem
  - relaxed
rate_gamma_prior_shape: 0.5
rate_gamma_prior_scale: 0.1
date_delimiter: "@"
date_index: -1
output:
  dir: results
  dated: false
beast:
  chain_length: 10000000
  samples: 10000
  seed: 12345
  duplicates: 2
  threads: 4
  args: "-beagle -beagle_CPU"
marginal_likelihood:
  estimate: false
```

Run it with:

```console
episodic run --config config.yaml
```

CLI options override values from the config file:

```console
episodic run --config config.yaml --beast-chain-length 1000000 --beast-samples 1000
```

## Input naming

Episodic extracts sampling dates from sequence names. By default, headers are split on `@`, and the last field is interpreted as the date:

```text
EPI_ISL_12345678@BA.2.86@2023.75
```

Use `--date-delimiter` and `--date-index` if your headers use a different format.

For FLC models, `--group` values are matched against sequence headers to define foreground lineages. At least one group is required when any selected clock contains `flc`.

## All Workflow Parameters

The table below lists the workflow config keys, their CLI flags, and defaults.

| Config key | CLI flag | Default | Description |
| --- | --- | --- | --- |
| `alignment` | `--alignment`, `-a` | `null` | Path to a FASTA alignment partition. Repeat for multiple partitions. Required for real runs. |
| `group` | `--group`, `-g` | `null` | Header-matching group string used to define FLC foreground lineages. Repeat for multiple groups. Required for FLC clocks. |
| `clock` | `--clock` | `flc-stem` | Clock model to run. Repeat for multiple models. Options: `strict`, `relaxed`, `flc-stem`, `flc-shared-stem`, `flc-clade`, `flc-shared-clade`, `flc-stem-and-clade`, `flc-shared-stem-and-clade`. |
| `rate_gamma_prior_shape` | `--rate-gamma-prior-shape`, `-shape` | `0.5` | Shape parameter for the gamma prior on clock rates. |
| `rate_gamma_prior_scale` | `--rate-gamma-prior-scale`, `-scale` | `0.1` | Scale parameter for the gamma prior on clock rates. |
| `mcc_tree.heights` | `--mcc-tree-heights` | `mean` | Node-height summary for MCC trees. Repeat for multiple values. Options include `mean`, `median`, `keep`, and `ca`. |
| `date_delimiter` | `--date-delimiter` | `@` | Delimiter used to split dates from sequence headers. |
| `date_index` | `--date-index` | `-1` | Zero-based field index containing the sampling date after splitting the header. `-1` means the last field. |
| `newick` | `--newick` | `null` | Optional Newick tree. If provided, the topology is fixed. |
| `output.dir` | `--output-dir`, `-o` | `.` | Output directory. |
| `output.dated` | `--output-dated`, `--no-output-dated` | `false` | Create a timestamped subdirectory under `output.dir`. |
| `beast.chain_length` | `--beast-chain-length`, `-cl` | `10000000` | BEAST MCMC chain length. |
| `beast.samples` | `--beast-samples`, `-s` | `10000` | Number of samples to draw from the BEAST chain. |
| `beast.seed` | `--beast-seed` | `null` | Optional base seed used to derive stable, job-specific BEAST seeds. |
| `beast.duplicates` | `--beast-duplicates` | `2` | Number of duplicate BEAST runs for convergence checks. |
| `beast.template` | `--beast-template` | bundled template | Optional Episodic BEAST XML template. |
| `beast.fit_clocks` | `--beast-fit-clocks`, `--no-beast-fit-clocks` | `true` | Run BEAST clock fitting. Disable to run only other requested workflow branches, such as MLE. |
| `beast.threads` | `--beast-threads` | `4` | Threads passed to BEAST. |
| `beast.args` | `--beast-args` | `-beagle -beagle_CPU` | Extra command-line arguments passed to BEAST. |
| `beast.envmodules` | `--beast-envmodules` | `GCC/11.3.0`, `beagle-lib/4.0.1-CUDA-12.2.0` | Environment modules to load for BEAST when Snakemake module loading is enabled. Repeat for multiple modules. |
| `marginal_likelihood.estimate` | `--marginal-likelihood-estimate`, `-mle` | `false` | Run path-sampling/stepping-stone marginal likelihood estimation. |
| `marginal_likelihood.path_steps` | `--marginal-likelihood-path-steps` | `100` | Number of path steps for marginal likelihood estimation. |
| `marginal_likelihood.chain_length` | `--marginal-likelihood-chain-length` | `1000000` | Chain length for marginal likelihood estimation. |
| `marginal_likelihood.log_every` | `--marginal-likelihood-log-every` | `10000` | Logging interval for marginal likelihood estimation. |
| `marginal_likelihood.duplicates` | `--marginal-likelihood-duplicates` | `3` | Number of duplicate marginal likelihood runs. |

## BEAST seeds

Episodic always passes an explicit `-seed` to BEAST. This avoids BEAST assigning the same time-based seed to jobs that start at the same time on a cluster.

Seeds are derived deterministically from:

- `beast.seed`, if provided, otherwise an internal default base value.
- The clock wildcard.
- The run name wildcard, including duplicate number.

This means duplicate jobs submitted at the same time receive different seeds, while rerunning the same analysis with the same config and output naming uses the same seeds. To start a new independent set of chains, change `beast.seed`:

```yaml
beast:
  seed: 20260506
```

## Outputs

By default, results are written to `output.dir`. If `output.dated: true`, a timestamped subdirectory is created under that directory.

For a full list of output files and naming patterns, see [Workflow Outputs](workflow.md).
