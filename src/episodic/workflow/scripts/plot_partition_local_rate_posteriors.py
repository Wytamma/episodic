"""Plot partition-level background vs foreground posterior substitution-rate histograms.

This script consumes one or more posterior sample tables (typically BEAST `.log` files)
and writes:

- an overlaid frequency histogram of per-partition background and foreground rates
- a long-format CSV table used to generate the plot

Supported input formats:

1. Wide BEAST-style columns, e.g.:
    - `<partition>.clock.rate` (background)
    - `<partition>.clade.rate`, `<partition>.stem.rate`, `<partition>.stem_and_clade.rate`, or equivalent

2. Long columns:
    - `partition`, `background_rate`, `clade_rate`

For current FLC logs used by Episodic, foreground-rate columns are treated as absolute rates.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, List, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
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


def _safe_label(label: Optional[str]) -> Optional[str]:
    if label is None:
        return None
    safe_label = re.sub(r"[^0-9A-Za-z_.-]+", "_", label.strip()).strip("._-")
    if not safe_label:
        return None
    if not safe_label[0].isalpha():
        safe_label = f"label_{safe_label}"
    return safe_label


def _strip_label_token(prefix: str, label: Optional[str]) -> str:
    safe_label = _safe_label(label)
    if safe_label is None:
        return prefix
    if prefix == safe_label:
        return ""
    if prefix.startswith(f"{safe_label}."):
        return prefix[len(safe_label) + 1 :]
    if prefix.endswith(f".{safe_label}"):
        return prefix[: -(len(safe_label) + 1)]
    return prefix


def _rate_key(column: str, suffix: str, label: Optional[str] = None) -> str:
    prefix = column[: -len(suffix)]
    return _strip_label_token(prefix, label)


def _foreground_rate_state(column: str, partition: str, suffix: str, foreground_label: Optional[str]) -> str:
    prefix = column[: -len(suffix)]
    safe_foreground = _safe_label(foreground_label)

    if safe_foreground and prefix.startswith(f"{safe_foreground}."):
        prefix = prefix[len(safe_foreground) + 1 :]

    if partition and prefix.startswith(f"{partition}."):
        state = prefix[len(partition) + 1 :]
    elif prefix == partition:
        state = ""
    else:
        state = prefix

    if safe_foreground and state.startswith(f"{safe_foreground}."):
        state = state[len(safe_foreground) + 1 :]

    return state or safe_foreground or "foreground"


def _display_state(state: str, foreground_label: Optional[str], background_label: Optional[str]) -> str:
    if state == "background" and background_label:
        return _safe_label(background_label) or "background"
    if state == "foreground" and foreground_label:
        return _safe_label(foreground_label) or "foreground"
    return state


def _read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".log", ".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, comment="#")


def _drop_burnin(df: pd.DataFrame, burnin: float) -> pd.DataFrame:
    start = int(len(df) * burnin)
    return df.iloc[start:].reset_index(drop=True)


def _extract_wide(
    df: pd.DataFrame,
    foreground_label: Optional[str] = None,
    background_label: Optional[str] = None,
) -> pd.DataFrame:
    background_cols: dict[str, str] = {}

    for column in df.columns:
        if column.endswith(".clock.rate"):
            background_cols[_rate_key(column, ".clock.rate", background_label or "background")] = column
        elif column.endswith(".background_rate"):
            background_cols[_rate_key(column, ".background_rate", background_label or "background")] = column

    foreground_suffixes = (".clade.rate", ".stem.rate", ".stem_and_clade.rate", ".clade_rate")
    foreground_by_partition: dict[str, list[str]] = {}
    for column in df.columns:
        for suffix in foreground_suffixes:
            if column.endswith(suffix) and not column.endswith(".clock.rate"):
                key = _rate_key(column, suffix, foreground_label)
                parts = key.split(".")
                while parts:
                    partition = ".".join(parts)
                    if partition in background_cols:
                        foreground_by_partition.setdefault(partition, []).append(column)
                        break
                    parts = parts[:-1]
                if not parts and len(background_cols) == 1:
                    partition = next(iter(background_cols))
                    foreground_by_partition.setdefault(partition, []).append(column)
                break

    partitions = sorted(partition for partition in background_cols if foreground_by_partition.get(partition))
    if not partitions:
        msg = (
            "Could not find matching partition background/foreground columns in wide format. "
            "Expected columns like '<partition>.clock.rate' with one of "
            "'<partition>.clade.rate', '<partition>.stem.rate', or "
            "'<partition>.stem_and_clade.rate', including group-specific variants."
        )
        raise ValueError(msg)

    rows: list[pd.DataFrame] = []
    for partition in partitions:
        background = df[background_cols[partition]].astype(float)

        rows.append(
            pd.DataFrame(
                {
                    "partition": partition,
                    "clock_state": "background",
                    "rate": background,
                }
            )
        )

        candidate_columns = sorted(set(foreground_by_partition[partition]))
        for column in candidate_columns:
            suffix = next(suffix for suffix in foreground_suffixes if column.endswith(suffix))
            foreground = df[column].astype(float)
            rows.append(
                pd.DataFrame(
                    {
                        "partition": partition,
                        "clock_state": _foreground_rate_state(column, partition, suffix, foreground_label),
                        "rate": foreground,
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
    long_df["foreground_rate"] = long_df["clade_rate"].astype(float)

    background = long_df[["partition", "background_rate"]].rename(columns={"background_rate": "rate"})
    background["clock_state"] = "background"

    foreground = long_df[["partition", "foreground_rate"]].rename(columns={"foreground_rate": "rate"})
    foreground["clock_state"] = "foreground"

    return pd.concat([background, foreground], ignore_index=True)


def _build_palette(
    partitions: Iterable[str],
    states: Iterable[str],
) -> dict[tuple[str, str], tuple[float, float, float]]:
    partition_list = list(partitions)
    state_list = list(states)
    state_palettes = ["Reds", "Blues", "Greens", "Purples", "Oranges", "Greys"]
    colors: dict[tuple[str, str], tuple[float, float, float]] = {}

    for state_idx, state in enumerate(state_list):
        palette_name = state_palettes[state_idx % len(state_palettes)]
        if len(partition_list) == 1:
            state_colors = sns.color_palette(palette_name, 3)[1:2]
        else:
            state_colors = sns.color_palette(palette_name, len(partition_list) + 3)[1:-1]
        for partition, color in zip(partition_list, state_colors):
            colors[(partition, state)] = color

    return colors


def _estimate_frequency_peak(samples: np.ndarray, bins: np.ndarray) -> tuple[float, float]:
    if len(samples) < 2:
        sample_val = float(samples[0]) if len(samples) == 1 else 0.0
        return sample_val, 0.0

    weights = np.ones(len(samples), dtype=float) / len(samples)
    frequencies, edges = np.histogram(samples, bins=bins, weights=weights)
    if len(frequencies) == 0:
        return float(np.median(samples)), 0.0

    return float(np.median(samples)), float(np.max(frequencies))


def _scientific_tick(value: float, _position: int) -> str:
    if abs(value) < 1e-15:
        return "0"
    return f"{value:.0e}"


def _histogram_bins(values: np.ndarray, bins: int) -> np.ndarray:
    positive_values = values[values > 0]
    if positive_values.size == values.size and positive_values.size > 0:
        xmin = float(positive_values.min()) * 0.75
        xmax = float(positive_values.max()) * 1.25
        if xmin > 0 and xmax > xmin:
            return np.logspace(np.log10(xmin), np.log10(xmax), bins + 1)
    return np.histogram_bin_edges(values, bins=bins)


def _frequency_ticks(ymax: float) -> tuple[float, list[float], list[str]]:
    tick_step = 0.05
    y_limit = max(tick_step, np.ceil(ymax / tick_step) * tick_step)
    ticks = np.arange(0, y_limit + tick_step / 2, tick_step).tolist()
    labels = ["0.0" if tick == 0 else f"{tick:.2f}" for tick in ticks]
    return float(y_limit), ticks, labels


def _rate_ticks(xmin: float, xmax: float) -> list[float]:
    preferred = [3e-4, 1e-3, 3e-3, 1e-2]
    ticks = [tick for tick in preferred if xmin < tick < xmax]
    if ticks:
        return ticks
    return []


def _legend_entries(
    plot_df: pd.DataFrame,
    partitions: list[str],
    states: list[str],
    palette: dict[tuple[str, str], tuple[float, float, float]],
    background_state: str,
) -> tuple[list[Patch], list[str]]:
    handles: list[Patch] = []
    labels: list[str] = []

    for state in states:
        handles.append(Patch(facecolor="none", edgecolor="none", alpha=0.0))
        labels.append(state)
        for partition in partitions:
            subset = plot_df[
                (plot_df["partition"] == partition) & (plot_df["clock_state"] == state)
            ]
            if subset.empty:
                continue
            color = palette[(partition, state)]
            handles.append(
                Patch(
                    facecolor=color,
                    edgecolor=color,
                    alpha=0.35 if state == background_state else 0.42,
                )
            )
            labels.append(partition)

    return handles, labels


@app.command()
def plot_partition_foreground_rates(
    posterior_samples: List[Path] = typer.Argument(..., help="Posterior samples files (BEAST .log or CSV/TSV)."),
    output_plot: Path = typer.Argument(..., help="Output SVG/PDF/PNG path for the frequency histogram."),
    output_table: Path = typer.Argument(..., help="Output CSV for long-format rates used in plotting."),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of each chain to discard as burn-in."),
    foreground_label: Optional[str] = typer.Option(None, help="Optional foreground/foreground-rate label for plot legends."),
    background_label: Optional[str] = typer.Option(None, help="Optional background-rate label for plot legends."),
    bins: int = typer.Option(220, "--bins", help="Number of histogram bins."),
):
    """Generate partitioned posterior frequency histograms for background vs foreground rates.

    Args:
        posterior_samples: One or more posterior sample files (`.log`, `.tsv`, `.csv`).
        output_plot: Output path for overlaid frequency histogram (`.svg`, `.png`, `.pdf`).
        output_table: Output path for long-format CSV with columns
            `partition`, `clock_state`, and `rate`.
        burnin: Fraction of each input chain to discard from the start.

    Raises:
        typer.BadParameter: If `burnin` is outside `[0, 1)`.
        ValueError: If required rate columns are missing or no posterior samples remain.
    """
    if burnin < 0 or burnin >= 1:
        raise typer.BadParameter("--burnin must be in [0, 1).")
    if bins < 1:
        raise typer.BadParameter("--bins must be >= 1.")

    dfs = []
    for path in posterior_samples:
        raw = _read_table(path)
        dfs.append(_drop_burnin(raw, burnin=burnin))

    combined = pd.concat(dfs, ignore_index=True)

    if {"partition", "background_rate", "clade_rate"}.issubset(combined.columns):
        plot_df = _extract_long(combined)
    else:
        plot_df = _extract_wide(
            combined,
            foreground_label=foreground_label,
            background_label=background_label,
        )

    plot_df = plot_df.dropna(subset=["partition", "clock_state", "rate"])
    if plot_df.empty:
        raise ValueError("No posterior rate samples available after filtering and burn-in.")

    plot_df["partition"] = plot_df["partition"].astype(str)
    plot_df["clock_state"] = plot_df["clock_state"].map(
        lambda state: _display_state(str(state), foreground_label, background_label)
    )
    background_state = _display_state("background", foreground_label, background_label)
    foreground_states = sorted(state for state in plot_df["clock_state"].unique() if state != background_state)
    state_order = [background_state, *foreground_states]
    plot_df["clock_state"] = pd.Categorical(plot_df["clock_state"], categories=state_order, ordered=True)

    partitions = sorted(plot_df["partition"].unique())
    palette = _build_palette(partitions, states=state_order)

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(11, 7.5))
    annotations = []
    finite_rates = plot_df["rate"].to_numpy(dtype=float)
    finite_rates = finite_rates[np.isfinite(finite_rates)]
    if finite_rates.size == 0:
        raise ValueError("No finite posterior rate samples available for plotting.")
    histogram_bins = _histogram_bins(finite_rates, bins=bins)

    for state in state_order:
        for partition in partitions:
            subset = plot_df[
                (plot_df["partition"] == partition) & (plot_df["clock_state"] == state)
            ]
            if subset.empty:
                continue
            samples = subset["rate"].to_numpy(dtype=float)
            samples = samples[np.isfinite(samples)]
            if samples.size == 0:
                continue
            weights = np.ones(samples.size, dtype=float) / samples.size
            ax.hist(
                samples,
                bins=histogram_bins,
                weights=weights,
                alpha=0.35 if state == background_state else 0.42,
                linewidth=0.18,
                color=palette[(partition, state)],
                edgecolor=palette[(partition, state)],
                label=f"{partition} ({state})",
            )

            peak_x, peak_frequency = _estimate_frequency_peak(samples, histogram_bins)
            annotations.append((partition, state, peak_x, peak_frequency))

    ymin, ymax = ax.get_ylim()
    y_span = max(ymax - ymin, 1e-9)
    for idx, (partition, state, peak_x, peak_frequency) in enumerate(annotations):
        y_text = peak_frequency + y_span * (0.02 + (idx % 3) * 0.015)
        ax.text(
            peak_x,
            y_text,
            partition,
            color=palette[(partition, state)],
            fontsize=12,
            fontweight="bold",
            ha="center",
            va="bottom",
        )

    ax.set_xlabel("substitutions/site/year", fontsize=19, fontweight="bold")
    ax.set_ylabel("Frequency", fontsize=16)
    ax.set_title("Partitioned fixed local clock posterior rates", fontsize=16)
    if np.all(finite_rates > 0):
        ax.set_xscale("log")
        ax.set_xlim(float(histogram_bins[0]), float(histogram_bins[-1]))
        ticks = _rate_ticks(float(histogram_bins[0]), float(histogram_bins[-1]))
        if ticks:
            ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(_scientific_tick))
    y_limit, y_ticks, y_tick_labels = _frequency_ticks(ax.get_ylim()[1])
    ax.set_ylim(0, y_limit)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels)
    ax.tick_params(axis="both", labelsize=13)
    ax.grid(True, which="major", alpha=0.22, linewidth=1.2)
    ax.grid(True, which="minor", axis="x", alpha=0.05)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    handles, labels = _legend_entries(plot_df, partitions, state_order, palette, background_state)
    legend = ax.legend(handles, labels, frameon=True, fontsize=9)
    for text in legend.get_texts():
        if text.get_text() in state_order:
            text.set_fontweight("bold")

    output_plot.parent.mkdir(parents=True, exist_ok=True)
    output_table.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(output_plot)
    plt.close(fig)

    plot_df.to_csv(output_table, index=False)


if __name__ == "__main__":
    app()
