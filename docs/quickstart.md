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
episodic run --fasta tests/data/.yaml
```


## Confuguration

The pipeline requires a configuration file to be provided. The configuration file is a YAML file that specifies the input data and the groups to test for episodic evolution.