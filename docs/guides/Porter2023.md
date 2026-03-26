# Reproducing Porter et al. (2023)

## Overview

This guide reproduces the mink host-jump analysis using fixed local clock (FLC) models, where branch groups are defined *a priori* and tested for episodic rate shifts.

![](/images/Porter2023_dag.png)

## Reference

- Porter et al. (2023)

## Inputs

- Alignment: `tests/data/Porter2023.afa`
- Starting tree: `tests/data/Porter2023.newick`
- Foreground groups: `mink_Netherlands`, `mink_Denmark`

## Analysis setup

We considered six FLC model configurations (following Tay et al. 2022):

- **FLC (stem):** elevated rate only on the stem branches leading into mink clades.
- **FLC (clade):** elevated rate only within mink clades.
- **FLC (clade and stem):** elevated rates on both stem branches and within mink clades.
- **FLC (shared, stem):** one shared stem-rate parameter across all mink clades.
- **FLC (shared, clade):** one shared clade-rate parameter across all mink clades.
- **FLC (shared, clade and stem):** shared rate parameters for stem and clade components across all mink clades.

The biological motivation is that adaptation to a new host may cause a short, episodic rate increase on stem branches (the transition period), with possible extension into the clade depending on the scenario. Across all models, the tree topology is fixed and the root age is constrained using prior knowledge that the SARS-CoV-2 tMRCA is in the second half of 2019 (Duchene et al. 2020; Ghafari et al. 2022). We used a gamma distribution with shape = 1 and scale = 0.001 such that the 95 per cent percentile range was 2.5 × 10–5 to 3.7 × 10–3.

## Template customization

Here we add a uniform prior on the root age to the default template, which is necessary to constrain the root age in this analysis. The template can be generated using `episodic template` and modified using `sed` or a text editor.

```bash
episodic template \
  | sed 's#</prior>#<uniformPrior lower="2019.5" upper="2020"><parameter idref="age(root)"/></uniformPrior>\
</prior>#' \
  > template_porter2023.xml
```

## Configure and Execute workflow

This run command executes the same dataset under six FLC parameterizations so their posterior rate estimates and model fit can be compared directly. The same alignment, fixed starting tree, and foreground groups are used across all model variants to isolate the effect of clock structure.

Key settings in this command are:

- `--beast-template template_porter2023.xml`: applies the customized root-age prior.
- `-a` and `--newick`: provide the sequence alignment and fixed topology.
- `-s 10000` and `-cl 10000000`: collects 10,000 samples from a 10,000,000-step chain.
- `--clock ...` (repeated): runs all six FLC model variants in one workflow.
- `-shape 1 -scale 0.001`: sets the gamma prior used for clock-rate parameters.
- `-g mink_Netherlands -g mink_Denmark`: defines mink foreground clades for local-clock assignment.


```console
episodic run \
  --beast-template template_porter2023.xml \
  -a tests/data/Porter2023.afa \
  --newick tests/data/Porter2023.newick \
  -s 10000 \
  -cl 10000000 \
  -o Porter2023 \
  --clock flc-stem \
  --clock flc-clade \
  --clock flc-stem-and-clade \
  --clock flc-shared-stem \
  --clock flc-shared-clade \
  --clock flc-shared-stem-and-clade \
  -shape 1 -scale 0.001 \
  -g mink_Netherlands \
  -g mink_Denmark
```

## Key outputs

In the clocks_1.0_0.001-violin-trimmed.png plot, the FLC stem models show a clear rate increase on the stem branches leading into mink clades, while the clade models show no increase within the mink clades. This suggests that the rate shift is concentrated on the transition branches rather than being sustained within the mink lineages, consistent with an episodic burst of adaptation during host-jump. The Netherlands mink clade shows a stronger rate increase than the Denmark clade, which is averaged out in the shared-parameter models.

![](/images/clocks_1.0_0.001-violin-trimmed.png)

The odds ratios and Bayes factors comparing the posterior probability of elevated rates in the foreground groups relative to the background are consistent with this interpretation: there is weak evidence against an elevated rate for the Denmark stem branch (BF = 0.724) and decisive evidence for an elevated rate in the Netherlands stem branch (BF = inf).

| Rate Column | Background Column | p_p | p_odds | pos_p | pos_odds | bf |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Porter2023.mink_Denmark.stem.rate | Porter2023.clock.rate | 0.48961226530385515 | 0.9592947322594689 | 0.6643150761026553 | 1.9789839483700151 | 0.7241404740143056 |
| Porter2023.mink_Netherlands.stem.rate | Porter2023.clock.rate | 0.48961226530385515 | 0.9592947322594689 | 1.0 | inf | inf |
