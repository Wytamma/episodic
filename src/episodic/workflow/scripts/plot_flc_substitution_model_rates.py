"""Plot posterior GTR substitution-model relative rates from FLC analyses.

The current Episodic BEAST template estimates one GTR relative-rate vector per
alignment partition. FLC models change branch clock rates, not branch-specific
substitution-model relative rates, so this plot shows the logged GTR relative
rates that are actually estimated for each partition.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import typer

app = typer.Typer()

SUBSTITUTION_ORDER = ["AC", "AG", "AT", "CG", "CT", "GT"]
TRANSVERSIONS = {"AC", "AT", "CG", "GT"}


def _read_log(path: Path, burnin: float) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", comment="#")
    return raw.iloc[int(len(raw) * burnin) :].reset_index(drop=True)


def _extract_relative_rates(df: pd.DataFrame) -> pd.DataFrame:
    pattern = re.compile(r"^(?P<partition>.+)\.gtr\.rates\.rate(?P<substitution>[ACGT]{2})$")
    rows: list[pd.DataFrame] = []

    for column in df.columns:
        match = pattern.match(column)
        if not match:
            continue

        values = df[column].astype(float)
        rows.append(
            pd.DataFrame(
                {
                    "partition": match.group("partition"),
                    "substitution": match.group("substitution"),
                    "log_relative_rate": np.log(values),
                }
            )
        )

    if not rows:
        raise ValueError("Could not find GTR relative-rate columns like '<partition>.gtr.rates.rateAC'.")

    rates = pd.concat(rows, ignore_index=True)
    rates = rates.replace([np.inf, -np.inf], np.nan).dropna(subset=["log_relative_rate"])
    if rates.empty:
        raise ValueError("No finite GTR relative-rate samples available for plotting.")
    return rates


def _display_substitution(substitution: str) -> str:
    if substitution in TRANSVERSIONS:
        return f"{substitution}*"
    return substitution


def _plot_partition(ax, partition_df: pd.DataFrame, partition: str) -> None:
    sns.boxplot(
        data=partition_df,
        x="substitution",
        y="log_relative_rate",
        order=SUBSTITUTION_ORDER,
        color="#d9d9d9",
        fliersize=0,
        linewidth=1,
        ax=ax,
    )
    ax.axhline(0, color="#777777", linestyle=(0, (2, 3)), linewidth=1)
    ax.set_title(partition, fontsize=11, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("log(relative rate)")
    ax.set_xticks(range(len(SUBSTITUTION_ORDER)))
    ax.set_xticklabels([_display_substitution(s) for s in SUBSTITUTION_ORDER], rotation=90)
    ax.grid(axis="y", alpha=0.2)


@app.command()
def plot_flc_substitution_model_rates(
    logs: List[Path] = typer.Argument(..., help="BEAST log files from a partitioned FLC analysis."),
    output_plot: Path = typer.Argument(..., help="Output SVG/PDF/PNG path."),
    output_table: Path = typer.Argument(..., help="Output CSV path for long-format logged relative rates."),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of each chain to discard as burn-in."),
):
    """Plot posterior log GTR relative rates by partition and substitution type."""
    if burnin < 0 or burnin >= 1:
        raise typer.BadParameter("--burnin must be in [0, 1).")

    df = pd.concat([_read_log(path, burnin=burnin) for path in logs], ignore_index=True)
    rates = _extract_relative_rates(df)

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    output_table.parent.mkdir(parents=True, exist_ok=True)
    rates.to_csv(output_table, index=False)

    sns.set_style("whitegrid")
    partitions = sorted(rates["partition"].unique())
    ncols = min(3, len(partitions))
    nrows = math.ceil(len(partitions) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(max(4, ncols * 3.4), max(3.4, nrows * 3.2)),
        squeeze=False,
    )

    for ax in axes.flatten():
        ax.set_visible(False)

    for ax, partition in zip(axes.flatten(), partitions):
        ax.set_visible(True)
        _plot_partition(ax, rates[rates["partition"] == partition], partition)

    fig.suptitle("GTR substitution-model relative rates", fontsize=13)
    fig.tight_layout()
    fig.savefig(output_plot)
    plt.close(fig)


if __name__ == "__main__":
    app()
