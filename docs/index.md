---
title: Home
---
# Episodic <img src='docs/images/logo.png' align="right" height="210" />

A complete pipeline for fitting and testing Fixed Local Clock (FLC) molecular clock models for episodic evolution.

[![PyPI - Version](https://img.shields.io/pypi/v/episodic.svg)](https://pypi.org/project/episodic)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/episodic.svg)](https://pypi.org/project/episodic)

-----

# About

Episodic is a tool for fitting and testing Fixed Local Clock (FLC) molecular clock models for episodic evolution. The package is built on top of [SNK](https://snk.wytamma.com/), and provides a complete pipeline for fitting and testing models of episodic evolution using [BEAST](https://beast.community/).

Episodic implements the ideas of Tay et al. ([2022](https://pubmed.ncbi.nlm.nih.gov/35038741/) and [2023](https://academic.oup.com/mbe/article/40/10/msad212/7280106)) and detects episodic evolution through Bayesian inference of molecular clock models. 

Given a multiple sequence alignment and a list of groups to test for episodic evolution, episodic will:
- Configure BEAST analyses for strict, relaxed (UCGD) and stem fixed local clock models. 
- Configure marginal likelihood analyses for each clock model.
- Run all the BEAST and marginal likelihood analyses.
- Plot and summarise the results.
- Compute and plot Bayes factors for the marginal likelihood analyses.
- Produce maximum clade credibility (MCC) trees for each clock model.
- Compute bayes factor on effect size for the FLC models (foreground vs background).
- Run rank and quantile tests on the all the models.
- Handel the execution of the pipeline on a HPC cluster via snakemake profiles.
- Produce a report of the results (TBD).

## Features

- **Complete pipeline** - `episodic` provides a complete pipeline for fitting and testing FLC models of episodic evolution.
- **Flexible** - `episodic` is built on top of [SNK](https://snk.wytamma.com/), and provides a flexible framework for fitting and testing FLC models of episodic evolution.
- **Easy to use** - `episodic` is easy to use, and provides a simple interface for fitting and testing FLC models of episodic evolution.
- **robust** - `episodic` is robust, and provides a robust framework for fitting and testing FLC models of episodic evolution. 