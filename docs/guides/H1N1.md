# Partitioned H1N1 strict-clock analysis

This guide documents an end-to-end workflow for running a partitioned H1N1 HA analysis with `episodic` using a strict molecular clock and no FLC foreground groups.

![](/images/clocks_0.5_0.1-violin-rug-trimmed.svg)


## Study context

This example is aligned with the BEAST partitioning workflow pattern:

- https://beast.community/constructing_models#importing-partitions-from-multiple-files

## Scope

This guide is intended as a reproducible template for standard strict-clock analyses where local-clock foreground groups are not required. We use two HA partitions (HA1 and HA2), run duplicate MCMC chains, and produce summary, trace, and MCC tree outputs.

We will do the following:

- run a strict clock analysis with partitioned alignments
- generate per-chain BEAST logs and trees
- generate MCC trees and render plots
- summarize posterior and trace diagnostics

## Phylogenetic analyses

We analyze H1N1 HA sequences in a Bayesian phylogenetic framework using BEAST 1.10 through the `episodic` workflow. The molecular clock model is strict clock (SC), with one rate parameter per partition in this partitioned setup.

The analysis uses partitioned alignments of HA1 and HA2, which are treated as separate partitions in the template. This yields partition-specific substitution/clock components while sharing the underlying tree process.

Because no FLC model is selected, no foreground group definitions are required. This is useful for baseline analyses and benchmarking before introducing local-clock hypotheses.

Posterior inference is performed by MCMC with user-defined chain length and sampling frequency (`beast.chain_length`, `beast.samples`). Duplicate runs are used to assess stability across independent chains. The workflow then generates trace plots, posterior summaries, and MCC trees for interpretation.

## Input data

This example starts from the NEXUS alignment:

- `tests/data/H1N1_HA.nex`

Taxon labels encode sampling time using `@`, and can include date uncertainty.
For example, `hSCar@1918/1` is interpreted as:

- date = `1918`
- uncertainty = `1` (for example, year-level uncertainty window)

The workflow parser supports this format (`DATE/UNCERTAINTY`) and reads the date field using `date_delimiter: '@'` and `date_index: -1`.

### Split NEXUS partitions into FASTA files

The NEXUS file contains partition definitions (`charset HA1`, `charset HA2`).
Use the partition-extraction script to produce one FASTA per partition:

```console
episodic script run extract_partitions_from_nexus --env python -- tests/data/H1N1_HA.nex --output-dir tests/data
```

This writes:

- `tests/data/HA1.fasta`
- `tests/data/HA2.fasta`


## Run

![](/images/H1N1_dag.png)

First run a dry run to validate the workflow configuration:

```console
episodic run -a tests/data/HA1.fasta -a tests/data/HA2.fasta -o H1N1 -cl 1000000 -s 1000 --clock strict --dry
```

Generate a DAG to visualize workflow structure:

```console
episodic run -a tests/data/HA1.fasta -a tests/data/HA2.fasta -o H1N1 -cl 1000000 -s 1000 --clock strict --dag H1N1_dag.pdf
```

Run the full analysis:

```console
episodic run -a tests/data/HA1.fasta -a tests/data/HA2.fasta -o H1N1 -cl 1000000 -s 1000 --clock strict
```

## Expected outputs

Main output categories:

- rendered run config and taxon grouping table
- per-clock BEAST logs/trees/XML/stdout
- MCC trees and MCC SVG renders
- per-clock summary CSV files
- trace plots and diagnostic figures

For the full output file matrix and naming patterns, see:

- [Workflow Outputs](../workflow.md)

### Compare traces from the two independent runs:

H1N1/clocks/strict_0.5_0.1/strict_0.5_0.1_1/strict_0.5_0.1_1_trace_plots/age(root).png vs H1N1/clocks/strict_0.5_0.1/strict_0.5_0.1_2/strict_0.5_0.1_2_trace_plots/age(root).png

![](/images/age(root)_1.png)
![](/images/age(root)_2.png)