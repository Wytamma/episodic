import csv
from bisect import bisect_left
from datetime import datetime
from typing import Dict, List

import dendropy
import matplotlib.pyplot as plt
import numpy as np
import typer

app = typer.Typer()

def extract_and_sort_rates(tree):
    """
    Extracts and sorts rates from a given tree.

    Args:
      tree (dendropy.Tree): The tree to extract rates from.

    Returns:
      List[float]: A list of sorted rates.

    Examples:
      >>> extract_and_sort_rates(tree)
      [0.1, 0.2, 0.3]
    """
    rates = [
        float(node.annotations.get_value("rate"))
        for node in tree
        if node.annotations.get_value("rate")
    ]
    sorted_rates = sorted(rates)
    return sorted_rates

def analyze_tree(tree, groups, group_stats):
    """
    Analyzes a given tree and updates group statistics.

    Args:
      tree (dendropy.Tree): The tree to analyze.
      groups (List[str]): The group labels to analyze.
      group_stats (Dict[str, Dict[str, List]]): A dictionary containing group statistics.

    Returns:
      None

    Examples:
      >>> analyze_tree(tree, ['A', 'B'], {'A': {'ranks': [], 'quantiles': []}, 'B': {'ranks': [], 'quantiles': []}})
    """
    # Assuming sorted_rates is generated here for each tree passed to this function
    sorted_rates = extract_and_sort_rates(tree)
    rates = np.array(sorted_rates)  # For efficient operations with numpy
    for group in groups:
        nodes_in_clade = [node for node in tree if node.taxon is not None and group in node.taxon.label]
        mrca = tree.mrca(taxon_labels=[node.taxon.label for node in nodes_in_clade])
        group_rate = float(mrca.annotations.get_value("rate"))

        # Use bisect_left for efficient rank finding in a sorted list
        rank = bisect_left(sorted_rates, group_rate) + 1
        group_stats[group]["ranks"].append(rank)

        # Use numpy for efficient percentile calculation
        quantile = np.percentile(rates, group_rate, method="lower")
        group_stats[group]["quantiles"].append(quantile)


@app.command()
def analyze_rates(
    trees_path: str = typer.Argument(..., help="Path to the BEAST output trees file"),
    groups: List[str] = typer.Option(..., "--group", "-g", help="Group labels to analyze"),
    output_plot_path: str = typer.Option(..., "--output-plot", help="Output path for the plot file"),
    output_csv_path: str = typer.Option(..., "--output-csv", help="Output path for the CSV file"),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of trees to discard as burn-in"),
):
    """
    Analyzes rates from a given BEAST output trees file and generates a plot and CSV file.

    Args:
      trees_path (str): The path to the BEAST output trees file.
      groups (List[str]): The group labels to analyze.
      output_plot_path (str): The output path for the plot file.
      output_csv_path (str): The output path for the CSV file.
      burnin (float): The fraction of trees to discard as burn-in.

    Returns:
      None

    Examples:
      >>> analyze_rates('trees.nexus', ['A', 'B'], 'plot.png', 'stats.csv', 0.1)
    """
    # time ow long it takes to run
    now = datetime.now()
    tree_yielder = dendropy.Tree.yield_from_files(files=[trees_path], schema="nexus", preserve_underscores=True)
    total_trees = sum(1 for _ in tree_yielder)  # Count total trees
    burnin_count = int(total_trees * burnin)
    total_time = datetime.now() - now
    print(f"Total time to count trees: {total_time}")

    tree_yielder = dendropy.Tree.yield_from_files(  # Reinitialize generator
        files=[trees_path], schema="nexus", preserve_underscores=True
    )
    group_stats: Dict[str, Dict[str, List]] = {g: {"ranks": [], "quantiles": []} for g in groups}

    with typer.progressbar(
            tree_yielder,
            length=total_trees,
            label="Processing trees",
            show_pos=True,
            show_percent=True
        ) as progress:
        for tree_idx, tree in enumerate(progress):
            if tree_idx < burnin_count:
                continue
            analyze_tree(tree, groups, group_stats)

    csv_data = [
        [
            "Group",
            "Mean Rank",
            "Rank Credible Interval",
            "Mean Quantile",
            "Quantile Credible Interval",
        ]
    ]

    plt.figure(figsize=(15, 5 * len(groups)))

    for i, group in enumerate(groups, start=1):
        ranks = group_stats[group]["ranks"]
        quantiles = group_stats[group]["quantiles"]

        mean_rank = np.mean(ranks)
        rank_credible_interval = (np.percentile(ranks, 2.5), np.percentile(ranks, 97.5))

        mean_quantile = np.mean(quantiles)
        quantile_credible_interval = (
            np.percentile(quantiles, 2.5),
            np.percentile(quantiles, 97.5),
        )

        csv_data.append(
            [
                group,
                f"{mean_rank:.4f}",
                f"[{rank_credible_interval[0]:.0f}, {rank_credible_interval[1]:.0f}]",
                f"{mean_quantile:.4f}",
                f"[{quantile_credible_interval[0]:.2f}, {quantile_credible_interval[1]:.2f}]",
            ]
        )

        plt.subplot(len(groups), 2, 2 * i - 1)
        plt.hist(ranks, bins=30, alpha=0.7, color="blue")
        rank_credible_interval_str = f"[{rank_credible_interval[0]:.0f}, {rank_credible_interval[1]:.0f}]"
        plt.title(f"Ranks for {group} - Mean: {mean_rank:.2f}, 95% CI: {rank_credible_interval_str}")
        plt.xlabel("Rank")
        plt.ylabel("Frequency")

        plt.subplot(len(groups), 2, 2 * i)
        plt.hist(quantiles, bins=30, alpha=0.7, color="green")
        quantile_credible_interval_str = f"[{quantile_credible_interval[0]:.2f}, {quantile_credible_interval[1]:.2f}]"
        plt.title(f"Quantiles for {group} - Mean: {mean_quantile:.2f}, 95% CI: {quantile_credible_interval_str}")
        plt.xlabel("Quantile")
        plt.ylabel("Frequency")

    plt.tight_layout()

    plt.savefig(output_plot_path)

    with open(output_csv_path, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerows(csv_data)


if __name__ == "__main__":
    app()
