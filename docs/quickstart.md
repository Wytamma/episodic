# QuickStart

Before running the pipeline it is sujjested that your read Tay et al. ([2022](https://pubmed.ncbi.nlm.nih.gov/35038741/) and [2023](https://academic.oup.com/mbe/article/40/10/msad212/7280106)) and familiarise yourself with the concepts of Fixed Local Clock (FLC) models.

## Installation

The episodic package can be installed via pip.

```console
pip install episodic
```

## Running the pipeline

The pipeline can be run using the `episodic` command line interface. 

```console
episodic run --alignment data.fasta --group group1 --group group2
```

## Where outputs are written

By default, results are written to `output.dir` from your config.
If `output.dated: true`, a timestamped subdirectory is created under that directory.

For a full list of output files and naming patterns, see:

- [Workflow Outputs](workflow.md)


## Confuguration

The pipeline requires a configuration file to be provided. The configuration file is a YAML file that specifies the input data and the groups to test for episodic evolution.