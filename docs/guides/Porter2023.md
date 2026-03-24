# Reproducing Porter et al. (2023)

In this guide, we focus on fixed local clock (FLC) models, which require *a priori* biological hypotheses about which branches should share a rate. Here, mink-associated branches are used to test whether host-jump events are associated with transient increases in evolutionary rate.

![](/images/Porter2023_dag.png)

We considered six FLC model configurations (following Tay et al. 2022):

- **FLC (stem):** elevated rate only on the stem branches leading into mink clades.
- **FLC (clade):** elevated rate only within mink clades.
- **FLC (clade and stem):** elevated rates on both stem branches and within mink clades.
- **FLC (shared, stem):** one shared stem-rate parameter across all mink clades.
- **FLC (shared, clade):** one shared clade-rate parameter across all mink clades.
- **FLC (shared, clade and stem):** shared rate parameters for stem and clade components across all mink clades.

The biological motivation is that adaptation to a new host may cause a short, episodic rate increase on stem branches (the transition period), with possible extension into the clade depending on the scenario. Across all models, the tree topology is fixed and the root age is constrained using prior knowledge that the SARS-CoV-2 tMRCA is in the second half of 2019 (Duchene et al. 2020; Ghafari et al. 2022). We used a gamma distribution with shape = 1 and scale = 0.01 such that the 95 per cent percentile range was 2.5 × 10–5 to 3.7 × 10–3. 

## Format template 

Here we add a uniform prior on the root age to the default template, which is necessary to constrain the root age in this analysis. The template can be generated using `episodic template` and modified using `sed` or a text editor.

```bash
episodic template \
  | sed 's#</prior>#<uniformPrior lower="2019.5" upper="2020"><parameter idref="age(root)"/></uniformPrior>\
</prior>#' \
  > template_porter2023.xml
```

## Example run

```bash
episodic run \
  --beast-template template_porter2023.xml \
  -a tests/data/Porter2023.afa \
  --newick tests/data/Porter2023.newick \
  -s 1000 \
  -cl 100000 \
  -o Porter2023 \
  --clock flc-stem \
  --clock flc-clade \
  --clock flc-stem-and-clade \
  --clock flc-shared-stem \
  --clock flc-shared-clade \
  --clock flc-shared-stem-and-clade \
  -shape 1 -scale 0.01 \
  -g mink_Netherlands \
  -g mink_Denmark --dag dag.png
```

## Outputs

![](/images/flc-stem-and-clade_0.5_0.1_2.mcc.mean.svg)