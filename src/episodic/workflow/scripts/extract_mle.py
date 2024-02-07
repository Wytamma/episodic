# CLI for plotting the marginal likelihoods from the mle files
# loop over the files in the directory and plot the mle
# plot using a jitter_stripplot from seaborn
# the files are named as follows: results/mle/flc-stem/flc-stem_mle_1.stdout
# the mle is on a line like this: log marginal likelihood (using stepping stone sampling) from pathLikelihood.delta = -143654.15496497147

import argparse
from pathlib import Path
import re
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description='Plot the marginal likelihoods from the mle files.')
parser.add_argument('directory', type=Path, help='Path to the directory containing the mle files.')

# Parse the command line arguments
args = parser.parse_args()

output = args.directory / 'mle'

# Set up the data frame
df = pd.DataFrame(columns=['Clock', 'log(MLE)'])

# Loop over the files in the directory

for file in args.directory.glob('**/*.stdout'):
    # Get the clock name
    clock = file.parent.name

    # Get the mle
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('log marginal likelihood'):
                mle = float(re.search(r'-?\d+\.\d+', line).group())

    # Add the data to the data frame
    df = pd.concat(
        [
            df,
            pd.DataFrame(
                {
                    'Clock': clock,
                    'log(MLE)': mle,
                },
                index=[0],
            ),
        ],
        ignore_index=True,
    )

sns.set_theme(style="whitegrid")

# Initialize the figure
fig, ax = plt.subplots()
sns.despine(bottom=True, left=True)

# Show each observation with a scatterplot
sns.stripplot(
    data=df, x='log(MLE)',
    y='Clock', hue="Clock",
    dodge=False, alpha=.7, zorder=1, legend=False,
)


# Show the conditional means, aligning each pointplot in the
# center of the strips by adjusting the width allotted to each
# category (.8 by default) by the number of hue levels
sns.pointplot(
    data=df, x='log(MLE)',
    y='Clock', hue="Clock",
    dodge=False, palette="dark", errorbar="sd",
    markers="d", markersize=4, linestyle="--", legend=False,
)

# Save the plot
fig.savefig(output.with_suffix('.svg'), dpi=300, bbox_inches='tight')

# save the data frame
df.to_csv(output.with_suffix('.csv'), index=False)

# group the data frame by clock and save the mean and standard deviation
df = df.groupby('Clock').agg(['mean', 'std'])
# reset the headers
df.columns = df.columns.droplevel()
# save the grouped data frame
df.to_csv(output.with_suffix('.grouped.csv'))