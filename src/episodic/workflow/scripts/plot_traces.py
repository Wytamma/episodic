import re
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from matplotlib.ticker import MaxNLocator, ScalarFormatter


DEFAULT_MAX_POINTS = 10000
DEFAULT_DPI = 300
VALID_FORMATS = {"png", "pdf", "svg"}

FIGURE_WIDTH = 13.5
FIGURE_HEIGHT = 5.8
POSTERIOR_COLOR = "#1f77b4"
POSTERIOR_FILL = "#dceaf7"
BURNIN_COLOR = "#c44e52"
MEAN_COLOR = "#2f2f2f"
MEDIAN_COLOR = "#117a65"
GRID_COLOR = "#e6e9ef"
TEXT_COLOR = "#1f2937"

TRACE_LOG_ARGUMENT = typer.Argument(..., help="The path to the trace log from BEAST.")
OUTPUT_ARGUMENT = typer.Argument(..., help="A path to the output directory.")
MAX_POINTS_OPTION = typer.Option(
    DEFAULT_MAX_POINTS,
    min=0,
    help="Maximum number of points to plot per burn-in/posterior panel. Use 0 to disable downsampling.",
)
DPI_OPTION = typer.Option(DEFAULT_DPI, min=72, help="Image DPI for raster outputs.")
DEBUG_OPTION = typer.Option(False, "--debug", help="Print timing information for each plotting step.")
FORMAT_OPTION = typer.Option(
    "png",
    "--format",
    help="Output format for plots. Choose from png, pdf, or svg.",
)


def debug_log(enabled: bool, message: str) -> None:
    """Print a debug message when enabled."""
    if enabled:
        typer.echo(f"[debug] {message}")


def camel_to_title_case(column: str) -> str:
    """Transform camelCase or dotted names into readable title case."""
    column = column.replace(".", " ")
    return re.sub(r"((?<=[a-z])[A-Z]|(?<!\A)[A-Z](?=[a-z]))", r" \1", column).title()


def split_trace(df: pd.DataFrame, burnin: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a trace dataframe into burn-in and posterior segments."""
    burnin_rows = int(len(df) * burnin)
    burnin_rows = min(max(burnin_rows, 0), len(df))
    return df.iloc[:burnin_rows], df.iloc[burnin_rows:]


def downsample_xy(x: pd.Series, y: pd.Series, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    """Evenly downsample paired x/y series for fast plotting."""
    if max_points <= 0 or len(x) <= max_points:
        return x.to_numpy(), y.to_numpy()

    indices = np.linspace(0, len(x) - 1, num=max_points, dtype=int)
    indices = np.unique(indices)
    return x.iloc[indices].to_numpy(), y.iloc[indices].to_numpy()


def get_axis_limits(values: pd.Series) -> tuple[float, float]:
    """Compute padded y-axis limits for a series."""
    value_min = float(values.min())
    value_max = float(values.max())
    padding = 0.1 * abs(value_max - value_min)

    if padding == 0:
        padding = 0.1 * abs(value_max) if value_max != 0 else 0.1

    return value_min - padding, value_max + padding


def style_axis(ax: plt.Axes) -> None:
    """Apply a publication-style theme to an axis."""
    ax.set_facecolor("white")
    ax.grid(True, color=GRID_COLOR, linewidth=0.8)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.spines["left"].set_color(TEXT_COLOR)
    ax.spines["bottom"].set_color(TEXT_COLOR)

    ax.tick_params(direction="out", colors=TEXT_COLOR, labelsize=10)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 4))
    ax.yaxis.set_major_formatter(formatter)


def add_posterior_summary(ax: plt.Axes, mean: float, median: float, hdi_low: float, hdi_high: float) -> None:
    """Add a compact summary annotation to the histogram panel."""
    summary = (
        f"mean   {mean:.3g}\n"
        f"median {median:.3g}\n"
        f"95% CI [{hdi_low:.3g}, {hdi_high:.3g}]"
    )
    ax.text(
        0.98,
        0.98,
        summary,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color=TEXT_COLOR,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#d1d5db", "alpha": 0.95},
    )


def plot_single_variable(
    output_path: Path,
    variable: str,
    burnin_df: pd.DataFrame,
    posterior_df: pd.DataFrame,
    max_points: int,
    output_format: str,
    dpi: int,
    debug: bool,
) -> None:
    """Render and save one trace plot."""
    variable_start = time.perf_counter()
    variable_title = camel_to_title_case(variable)
    typer.echo(f"Plotting {variable_title}")
    debug_log(debug, f"[{variable}] start")

    stats_start = time.perf_counter()
    posterior_series = posterior_df[variable]
    burnin_series = burnin_df[variable]
    burnin_y_min, burnin_y_max = get_axis_limits(burnin_series) if not burnin_series.empty else get_axis_limits(posterior_series)
    posterior_y_min, posterior_y_max = get_axis_limits(posterior_series)
    mean = float(posterior_series.mean())
    median = float(posterior_series.median())
    hdi_low, hdi_high = np.quantile(posterior_series.to_numpy(), [0.025, 0.975])
    debug_log(debug, f"[{variable}] computed stats in {time.perf_counter() - stats_start:.2f}s")

    sample_start = time.perf_counter()
    burnin_x, burnin_y = downsample_xy(burnin_df["state"], burnin_series, max_points)
    posterior_x, posterior_y = downsample_xy(posterior_df["state"], posterior_series, max_points)
    debug_log(debug, f"[{variable}] downsampled data in {time.perf_counter() - sample_start:.2f}s")

    figure_start = time.perf_counter()
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(FIGURE_WIDTH, FIGURE_HEIGHT),
        gridspec_kw={"width_ratios": [1, 3, 1]},
        constrained_layout=True,
    )
    burnin_ax, posterior_ax, hist_ax = axes

    fig.patch.set_facecolor("white")
    fig.suptitle(variable_title, fontsize=18, color=TEXT_COLOR, fontweight="semibold")

    for ax in axes:
        style_axis(ax)

    burnin_ax.set_ylim(burnin_y_min, burnin_y_max)
    posterior_ax.set_ylim(posterior_y_min, posterior_y_max)
    hist_ax.set_ylim(posterior_y_min, posterior_y_max)

    burnin_ax.axhspan(burnin_y_min, burnin_y_max, color=POSTERIOR_FILL, zorder=0)
    posterior_ax.axhspan(posterior_y_min, posterior_y_max, color=POSTERIOR_FILL, zorder=0)
    posterior_ax.axhspan(hdi_low, hdi_high, color="#c6dbef", alpha=0.65, zorder=0)

    if len(burnin_x) > 0:
        burnin_ax.axvspan(float(burnin_x.min()), float(burnin_x.max()), color=BURNIN_COLOR, alpha=0.12, zorder=0)
        burnin_ax.plot(burnin_x, burnin_y, color=BURNIN_COLOR, linewidth=1.1, alpha=0.95)
        if len(burnin_x) <= 1000:
            burnin_ax.scatter(burnin_x, burnin_y, color=BURNIN_COLOR, s=7, alpha=0.45, linewidths=0)

    posterior_ax.plot(posterior_x, posterior_y, color=POSTERIOR_COLOR, linewidth=1.1, alpha=0.95)
    if len(posterior_x) <= 1000:
        posterior_ax.scatter(posterior_x, posterior_y, color=POSTERIOR_COLOR, s=7, alpha=0.4, linewidths=0)

    posterior_ax.axhline(mean, color=MEAN_COLOR, linewidth=1.2, linestyle="--", alpha=0.9)
    posterior_ax.axhline(median, color=MEDIAN_COLOR, linewidth=1.2, linestyle=":", alpha=0.9)

    hist_ax.hist(
        posterior_series.to_numpy(),
        bins=60,
        orientation="horizontal",
        density=True,
        color=POSTERIOR_COLOR,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.6,
    )
    hist_ax.axhline(mean, color=MEAN_COLOR, linewidth=1.2, linestyle="--")
    hist_ax.axhline(median, color=MEDIAN_COLOR, linewidth=1.2, linestyle=":")
    hist_ax.axhspan(hdi_low, hdi_high, color="#c6dbef", alpha=0.35, zorder=0)
    add_posterior_summary(hist_ax, mean, median, hdi_low, hdi_high)
    hist_ax.set_xlim(left=0)

    burnin_ax.set_title("Burn-in", fontsize=12, color=TEXT_COLOR, fontweight="semibold")
    posterior_ax.set_title("Posterior", fontsize=12, color=TEXT_COLOR, fontweight="semibold")
    hist_ax.set_title("Posterior density", fontsize=12, color=TEXT_COLOR, fontweight="semibold")

    burnin_ax.set_xlabel("MCMC state", color=TEXT_COLOR)
    posterior_ax.set_xlabel("MCMC state", color=TEXT_COLOR)
    hist_ax.set_xlabel("Density", color=TEXT_COLOR)
    burnin_ax.set_ylabel(variable_title)
    posterior_ax.set_ylabel("")
    hist_ax.set_ylabel("")

    burnin_ax.text(
        0.02,
        0.98,
        f"n={len(burnin_df):,}",
        transform=burnin_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=TEXT_COLOR,
    )
    posterior_ax.text(
        0.02,
        0.98,
        f"n={len(posterior_df):,}",
        transform=posterior_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=TEXT_COLOR,
    )

    debug_log(debug, f"[{variable}] built matplotlib figure in {time.perf_counter() - figure_start:.2f}s")

    write_start = time.perf_counter()
    save_kwargs = {"bbox_inches": "tight"}
    if output_format == "png":
        save_kwargs["dpi"] = dpi
    fig.savefig(output_path, **save_kwargs)
    plt.close(fig)
    debug_log(debug, f"[{variable}] wrote figure in {time.perf_counter() - write_start:.2f}s")
    debug_log(debug, f"[{variable}] total {time.perf_counter() - variable_start:.2f}s")


def plot_traces(
    trace_log: Path = TRACE_LOG_ARGUMENT,
    output: Path = OUTPUT_ARGUMENT,
    burnin: float = 0.1,
    max_points: int = MAX_POINTS_OPTION,
    dpi: int = DPI_OPTION,
    debug: bool = DEBUG_OPTION,
    output_format: str = FORMAT_OPTION,
) -> None:
    """Produce fast publication-ready trace plots using Matplotlib."""
    start_time = time.perf_counter()
    output.mkdir(exist_ok=True, parents=True)

    output_format = output_format.lower()
    if output_format not in VALID_FORMATS:
        msg = f"Unsupported output format '{output_format}'. Choose from: {', '.join(sorted(VALID_FORMATS))}."
        raise typer.BadParameter(msg)

    debug_log(debug, f"Output directory: {output}")
    debug_log(debug, f"Output format: {output_format}")

    read_start = time.perf_counter()
    df = pd.read_csv(trace_log, sep="\t", comment="#")
    debug_log(debug, f"Loaded trace log with shape {df.shape} in {time.perf_counter() - read_start:.2f}s")

    if "state" not in df.columns:
        msg = "Trace log must contain a 'state' column."
        raise typer.BadParameter(msg)

    split_start = time.perf_counter()
    burnin_df, posterior_df = split_trace(df, burnin)
    debug_log(
        debug,
        (
            f"Split trace in {time.perf_counter() - split_start:.2f}s: "
            f"burn-in rows={len(burnin_df)}, posterior rows={len(posterior_df)}"
        ),
    )

    if posterior_df.empty:
        msg = "Posterior trace is empty after burn-in removal. Reduce burn-in or provide a longer chain."
        raise typer.BadParameter(msg)

    variables = [column for column in df.columns if column != "state"]
    debug_log(debug, f"Preparing plots for {len(variables)} variables")

    for variable in variables:
        output_path = output / f"{variable}.{output_format}"
        plot_single_variable(
            output_path=output_path,
            variable=variable,
            burnin_df=burnin_df,
            posterior_df=posterior_df,
            max_points=max_points,
            output_format=output_format,
            dpi=dpi,
            debug=debug,
        )

    debug_log(debug, f"Finished all plots in {time.perf_counter() - start_time:.2f}s")


if __name__ == "__main__":
    typer.run(plot_traces)
