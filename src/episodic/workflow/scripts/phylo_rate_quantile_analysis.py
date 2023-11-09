from typing import List, Dict
import csv

import dendropy
import matplotlib.pyplot as plt
import numpy as np
import typer
from scipy.stats import percentileofscore

app = typer.Typer()

@app.command()
def analyze_rates(
    trees_path: str = typer.Argument(..., help="Path to the BEAST output trees file"),
    groups: List[str] = typer.Option(..., "--group", "-g", help="Group labels to analyze"),
    output_plot_path: str = typer.Option(..., "--output-plot", help="Output path for the plot file"),
    output_csv_path: str = typer.Option(..., "--output-csv", help="Output path for the CSV file"),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of trees to discard as burn-in")
):
    # Load the BEAST trees file
    trees = dendropy.TreeList.get(path=trees_path, schema="nexus", preserve_underscores=True)

    # discard burn-in
    trees = trees[int(len(trees) * burnin):]

    print(f"Loaded {len(trees)} trees from {trees_path}")

    # Initialize a dictionary to hold ranks and quantiles for each group
    group_stats: Dict[str, Dict[str, List]] = {g: {"ranks": [], "quantiles": []} for g in groups}

    # Iterate over each tree in the trees file
    for tree in trees:
        # Extract rates from tree annotations)
        rates = [float(node.annotations.get_value("rate")) for node in tree if node.annotations.get_value("rate")]
        sorted_rates = sorted(rates)
        # For each group, find the rate and determine its rank and quantile
        for group in groups:
            nodes_in_clade = [node for node in tree if node.taxon is not None and group in node.taxon.label]
            # get MRCA of all nodes in clade
            mrca = tree.mrca(taxon_labels=[node.taxon.label for node in nodes_in_clade])
            group_rate = float(mrca.annotations.get_value("rate"))
            rank = sorted_rates.index(group_rate) + 1  # ranks are 1-based
            group_stats[group]["ranks"].append(rank)
            quantile = percentileofscore(rates, group_rate, "rank")
            group_stats[group]["quantiles"].append(quantile)

    # Prepare the CSV data
    csv_data = [["Group", "Mean Rank", "Rank Credible Interval", "Mean Quantile", "Quantile Credible Interval"]]

    # Plot settings
    plt.figure(figsize=(15, 5 * len(groups)))

    # Calculate and plot the mean and 95% credible interval for the ranks and quantiles for each group
    for i, group in enumerate(groups, start=1):
        ranks = group_stats[group]["ranks"]
        quantiles = group_stats[group]["quantiles"]

        # Rank statistics
        mean_rank = np.mean(ranks)
        rank_credible_interval = (np.percentile(ranks, 2.5), np.percentile(ranks, 97.5))

        # Quantile statistics
        mean_quantile = np.mean(quantiles)
        quantile_credible_interval = (np.percentile(quantiles, 2.5), np.percentile(quantiles, 97.5))

        # Update the CSV data
        csv_data.append([
            group,
            f"{mean_rank:.4f}",
            f"[{rank_credible_interval[0]:.0f}, {rank_credible_interval[1]:.0f}]",
            f"{mean_quantile:.4f}",
            f"[{quantile_credible_interval[0]:.2f}, {quantile_credible_interval[1]:.2f}]"
        ])

        # Plot the rank distribution
        plt.subplot(len(groups), 2, 2 * i - 1)
        plt.hist(ranks, bins=30, alpha=0.7, color="blue")
        rank_credible_interval_str = f"[{rank_credible_interval[0]:.0f}, {rank_credible_interval[1]:.0f}]"
        plt.title(f"Ranks for {group} - Mean: {mean_rank:.2f}, 95% CI: {rank_credible_interval_str}")
        plt.xlabel("Rank")
        plt.ylabel("Frequency")

        # Plot the quantile distribution
        plt.subplot(len(groups), 2, 2 * i)
        plt.hist(quantiles, bins=30, alpha=0.7, color="green")
        quantile_credible_interval_str = f"[{quantile_credible_interval[0]:.2f}, {quantile_credible_interval[1]:.2f}]"
        plt.title(f"Quantiles for {group} - Mean: {mean_quantile:.2f}, 95% CI: {quantile_credible_interval_str}")
        plt.xlabel("Quantile")
        plt.ylabel("Frequency")

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the plot to a file instead of showing it
    plt.savefig(output_plot_path)

    # Save the results to a CSV file
    with open(output_csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerows(csv_data)

if __name__ == "__main__":
    app()
