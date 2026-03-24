"""Plot partition-level background vs local posterior substitution-rate densities.

This script consumes one or more posterior sample tables (typically BEAST `.log` files)
and writes:

- an overlaid density plot of per-partition background and local rates
- a long-format CSV table used to generate the plot

Supported input formats:

1. Wide BEAST-style columns, e.g.:
    - `<partition>.clock.rate` (background)
    - `<partition>.clade.rate`, `<partition>.stem.rate`, `<partition>.stem_and_clade.rate`, or equivalent

2. Long columns:
    - `partition`, `background_rate`, `clade_rate`

For current FLC logs used by Episodic, local-rate columns are treated as absolute rates.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import typer

app = typer.Typer()

BACKGROUND_PATTERNS = [
    re.compile(r"^(?P<partition>.+)\.background_rate$"),
    re.compile(r"^(?P<partition>.+)\.clock\.rate$"),
]

CLADE_PATTERNS = [
    re.compile(r"^(?P<partition>.+)\.clade_rate$"),
    re.compile(r"^(?P<partition>.+)\.clade\.rate$"),
    re.compile(r"^(?P<partition>.+)\.stem\.rate$"),
    re.compile(r"^(?P<partition>.+)\.stem_and_clade\.rate$"),
    re.compile(r"^(?P<partition>.+)\.[^.]+\.clade\.rate$"),
    re.compile(r"^(?P<partition>.+)\.[^.]+\.stem\.rate$"),
    re.compile(r"^(?P<partition>.+)\.[^.]+\.stem_and_clade\.rate$"),
]


def _read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".log", ".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, comment="#")


def _drop_burnin(df: pd.DataFrame, burnin: float) -> pd.DataFrame:
    start = int(len(df) * burnin)
    return df.iloc[start:].reset_index(drop=True)


def _extract_wide(df: pd.DataFrame) -> pd.DataFrame:
    background_cols: dict[str, str] = {}
    clade_cols: dict[str, str] = {}

    for column in df.columns:
        for pattern in BACKGROUND_PATTERNS:
            match = pattern.match(column)
            if match:
                background_cols[match.group("partition")] = column
                break

    for partition in sorted(background_cols):
        candidate_columns = [
            c
            for c in df.columns
            if c.startswith(f"{partition}.")
            and (
                c.endswith(".clade.rate")
                or c.endswith(".stem.rate")
                or c.endswith(".stem_and_clade.rate")
                or c.endswith(".clade_rate")
            )
        ]

        if not candidate_columns:
            continue

        priority = [
            f"{partition}.clade.rate",
            f"{partition}.stem.rate",
            f"{partition}.stem_and_clade.rate",
            f"{partition}.clade_rate",
        ]

        selected = None
        for name in priority:
            if name in candidate_columns:
                selected = name
                break
        if selected is None:
            selected = sorted(candidate_columns)[0]

        clade_cols[partition] = selected

    partitions = sorted(set(background_cols) & set(clade_cols))
    if not partitions:
        msg = (
            "Could not find matching partition background/clade columns in wide format. "
            "Expected columns like '<partition>.clock.rate' with one of "
            "'<partition>.clade.rate', '<partition>.stem.rate', or "
            "'<partition>.stem_and_clade.rate'."
        )
        raise ValueError(msg)

    rows: list[pd.DataFrame] = []
    for partition in partitions:
        background = df[background_cols[partition]].astype(float)
        local = df[clade_cols[partition]].astype(float)

        rows.append(
            pd.DataFrame(
                {
                    "partition": partition,
                    "clock_state": "background",
                    "rate": background,
                }
            )
        )
        rows.append(
            pd.DataFrame(
                {
                    "partition": partition,
                    "clock_state": "local",
                    "rate": local,
                }
            )
        )

    return pd.concat(rows, ignore_index=True)


def _extract_long(df: pd.DataFrame) -> pd.DataFrame:
    required = {"partition", "background_rate", "clade_rate"}
    if not required.issubset(df.columns):
        missing = sorted(required.difference(df.columns))
        raise ValueError(
            f"Long format missing required columns: {', '.join(missing)}. "
            "Expected: partition, background_rate, clade_rate"
        )

    long_df = df[["partition", "background_rate", "clade_rate"]].copy()
    long_df["background_rate"] = long_df["background_rate"].astype(float)
    long_df["local_rate"] = long_df["clade_rate"].astype(float)

    background = long_df[["partition", "background_rate"]].rename(columns={"background_rate": "rate"})
    background["clock_state"] = "background"

    local = long_df[["partition", "local_rate"]].rename(columns={"local_rate": "rate"})
    local["clock_state"] = "local"

    return pd.concat([background, local], ignore_index=True)


def _build_palette(partitions: Iterable[str]) -> dict[tuple[str, str], tuple[float, float, float]]:
    partition_list = list(partitions)

    if len(partition_list) == 1:
        background_palette = sns.color_palette("YlOrBr", 3)[1:2]
        local_palette = sns.color_palette("GnBu", 3)[1:2]
    else:
        background_palette = sns.color_palette("YlOrBr", len(partition_list) + 2)[1:-1]
        local_palette = sns.color_palette("GnBu", len(partition_list) + 2)[1:-1]

    colors: dict[tuple[str, str], tuple[float, float, float]] = {}

    for partition, background_color, local_color in zip(
        partition_list, background_palette, local_palette
    ):
        colors[(partition, "background")] = background_color
        colors[(partition, "local")] = local_color

    return colors


def _estimate_density_peak(samples: np.ndarray) -> tuple[float, float]:
    if len(samples) < 2:
        sample_val = float(samples[0]) if len(samples) == 1 else 0.0
        return sample_val, 0.0

    sigma = float(np.std(samples, ddof=1))
    if not np.isfinite(sigma) or sigma <= 0:
        median_val = float(np.median(samples))
        return median_val, 0.0

    bandwidth = 1.06 * sigma * (len(samples) ** (-1.0 / 5.0))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        median_val = float(np.median(samples))
        return median_val, 0.0

    x_min = float(np.min(samples))
    x_max = float(np.max(samples))
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
        return float(np.median(samples)), 0.0

    grid = np.linspace(x_min, x_max, 256)
    z_scores = (grid[:, None] - samples[None, :]) / bandwidth
    kernel_vals = np.exp(-0.5 * z_scores * z_scores) / (np.sqrt(2 * np.pi) * bandwidth)
    densities = np.mean(kernel_vals, axis=1)

    peak_index = int(np.argmax(densities))
    return float(grid[peak_index]), float(densities[peak_index])


@app.command()
def plot_partition_local_rates(
    posterior_samples: List[Path] = typer.Argument(..., help="Posterior samples files (BEAST .log or CSV/TSV)."),
    output_plot: Path = typer.Argument(..., help="Output SVG/PDF/PNG path for the density plot."),
    output_table: Path = typer.Argument(..., help="Output CSV for long-format rates used in plotting."),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of each chain to discard as burn-in."),
):
    """Generate partitioned posterior density plot for background vs local rates.

    Args:
        posterior_samples: One or more posterior sample files (`.log`, `.tsv`, `.csv`).
        output_plot: Output path for overlaid density plot (`.svg`, `.png`, `.pdf`).
        output_table: Output path for long-format CSV with columns
            `partition`, `clock_state`, and `rate`.
        burnin: Fraction of each input chain to discard from the start.

    Raises:
        typer.BadParameter: If `burnin` is outside `[0, 1)`.
        ValueError: If required rate columns are missing or no posterior samples remain.
    """
    if burnin < 0 or burnin >= 1:
        raise typer.BadParameter("--burnin must be in [0, 1).")

    dfs = []
    for path in posterior_samples:
        raw = _read_table(path)
        dfs.append(_drop_burnin(raw, burnin=burnin))

    combined = pd.concat(dfs, ignore_index=True)

    if {"partition", "background_rate", "clade_rate"}.issubset(combined.columns):
        plot_df = _extract_long(combined)
    else:
        plot_df = _extract_wide(combined)

    plot_df = plot_df.dropna(subset=["partition", "clock_state", "rate"])
    if plot_df.empty:
        raise ValueError("No posterior rate samples available after filtering and burn-in.")

    plot_df["partition"] = plot_df["partition"].astype(str)
    plot_df["clock_state"] = pd.Categorical(
        plot_df["clock_state"], categories=["background", "local"], ordered=True
    )

    partitions = sorted(plot_df["partition"].unique())
    palette = _build_palette(partitions)

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(11, 6.5))
    annotations = []

    for partition in partitions:
        for state in ("background", "local"):
            subset = plot_df[
                (plot_df["partition"] == partition) & (plot_df["clock_state"] == state)
            ]
            if subset.empty:
                continue
            sns.kdeplot(
                data=subset,
                x="rate",
                fill=True,
                common_norm=False,
                alpha=0.3,
                linewidth=1.1,
                color=palette[(partition, state)],
                label=f"{partition} ({state})",
                ax=ax,
            )

            samples = subset["rate"].to_numpy(dtype=float)
            peak_x, peak_density = _estimate_density_peak(samples)
            annotations.append((partition, state, peak_x, peak_density))

    ymin, ymax = ax.get_ylim()
    y_span = max(ymax - ymin, 1e-9)
    for idx, (partition, state, peak_x, peak_density) in enumerate(annotations):
        y_text = peak_density + y_span * (0.02 + (idx % 3) * 0.015)
        label = f"{partition} ({state})"
        ax.text(
            peak_x,
            y_text,
            label,
            color=palette[(partition, state)],
            fontsize=12,
            ha="left",
            va="bottom",
        )

    ax.set_xlabel("substitutions per site per year")
    ax.set_ylabel("density")
    ax.set_title("Partitioned fixed local clock posterior rates")
    ax.legend(frameon=True, fontsize=9)

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    output_table.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(output_plot)
    plt.close(fig)

    plot_df.to_csv(output_table, index=False)


if __name__ == "__main__":
    app()
