from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import List

import arviz as az
import bokeh.io.showing
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
import xarray as xr
from bokeh.embed import file_html
from bokeh.io.showing import _show_file_with_state
from bokeh.plotting import figure, output_file, save
from bokeh.resources import CDN

app = typer.Typer()


@contextmanager
def no_browser():
    """
    A context manager that temporarily replaces the Bokeh show function with a dummy one.

    Yields:
      None: Allows the code inside the 'with' block to run.

    Notes:
      This context manager is used to prevent the Bokeh show function from opening a browser window when saving a plot.
    """
    # Save the original show function
    original_show = bokeh.io.showing._show_file_with_state

    # Define a custom show function that does nothing
    def dummy_show(obj, state, *args, **kwargs):
        """
    A custom show function that does nothing.

    Args:
      obj (object): The object to show.
      state (object): The state of the object.
      *args: Additional positional arguments.
      **kwargs: Additional keyword arguments.

    Returns:
      None: Does not return anything.
    """
        filename = save(obj, state=state)

    # Replace the Bokeh show function with the dummy one
    bokeh.io.showing._show_file_with_state = dummy_show

    try:
        yield  # This allows the code inside the 'with' block to run
    finally:
        # Restore the original show function
        bokeh.io.showing._show_file_with_state = original_show

def load_log_files(logs: List[Path], burnin: float = 0.1) -> pd.DataFrame:
    """
    Loads BEAST log files into a single pandas DataFrame.

    Args:
      logs (List[Path]): A list of paths to the log files.
      burnin (float): The fraction of the chain to discard as burnin.

    Returns:
      pd.DataFrame: A DataFrame containing the log data.
    """
    model_groups = defaultdict(list)
    for path in logs:
        model = path.parent.parent.name # assumes that the model name is the parent of the parent of the log file
        model_groups[model].append(path)

    dfs = []
    for model, paths in model_groups.items():
        model_df = pd.DataFrame()
        chain_count = 0
        for trace_log in paths:
            print(trace_log)
            duplicate_df = pd.read_csv(trace_log, sep="\t", comment="#").rename(columns={"state": "draw"})
            posterior_df = duplicate_df.truncate(before=burnin * len(duplicate_df))
            posterior_df["chain"] = chain_count
            model_df = pd.concat([model_df, posterior_df])
            chain_count += 1
        if len(model_groups) > 1:
            # rename the columns to include the model name when there are multiple models
            var_names = [c for c in model_df.columns if c not in ("chain", "draw")]
            rename_dict = {col: f"{model}.{col}" for col in var_names}
            model_df.rename(columns=rename_dict, inplace=True)
        if len(dfs):
            # drop the chain and draw columns from all but the first model
            model_df.drop(columns=["chain", "draw"], inplace=True)
        dfs.append(model_df)
    df = pd.concat(dfs, axis=1)
    df = df.set_index(["chain", "draw"])
    return df


@app.command()
def rates(
    logs: List[Path] = typer.Argument(..., help="BEAST log files"),
    output_prefix: Path = typer.Option(..., help="Prefix for output files"),
    gamma_shape: float = typer.Option(..., help="Shape parameter for the gamma prior"),
    gamma_scale: float = typer.Option(..., help="Scale parameter for the gamma prior"),
    burnin: float = typer.Option(0.1, help="Fraction of the chain to discard as burnin"),
):
    """
    Plots the rates from BEAST log files.

    Args:
      logs (List[Path]): A list of paths to the log files.
      output_prefix (Path): The prefix for the output files.
      gamma_shape (float): The shape parameter for the gamma prior.
      gamma_scale (float): The scale parameter for the gamma prior.
      burnin (float): The fraction of the chain to discard as burnin.

    Returns:
      None: Does not return anything.
    """
    df = load_log_files(logs, burnin=burnin)

    # extract the rate columns
    var_names = [c for c in df.columns if c.endswith(".rate") or c.endswith(".ucgd.mean")]
    df = df[var_names]


    # add a prior rate column
    prior_rate = np.random.gamma(gamma_shape, gamma_scale, len(df))
    df["prior.rate"] = prior_rate

    # convert to xarray
    xdata = xr.Dataset.from_dataframe(df)
    dataset = az.InferenceData(posterior=xdata)

    for rug in (True, False):
        rug_str = "violin-rug" if rug else "violin"
        # plot the rates
        axs = az.plot_violin(
            dataset,
            figsize=(len(df.columns) * 5, 12),
            textsize=16,
            sharey=True,
            sharex=False,
            rug=rug,
            grid=(1, len(df.columns)),
        )
        plt.savefig(f"{output_prefix}-{rug_str}.svg")

        selected_columns = df.columns[~df.columns.isin(["prior.rate", "draw", "chain"])]

        for ax in axs.flatten():
            ymax = df[selected_columns].max().max() + df[selected_columns].min().std()
            ymin = df[selected_columns].min().min() - df[selected_columns].min().std()
            ax.set_ylim(ymin, ymax)

        plt.savefig(f"{output_prefix}-{rug_str}-trimmed.svg")

    # plot the trace
    az.plot_trace(
        dataset,
        figsize=(12, len(df.columns) * 4),
        )
    plt.savefig(f"{output_prefix}-trace.svg")

    # plot interactive trace
    output_file(filename=f"{output_prefix}-trace.html")
    with no_browser():
        # hack to save the html file without opening it in the browser
        # must set show=True to save the file
        az.plot_trace(dataset, backend="bokeh", show=True)

@app.command()
def trace(
        logs: List[Path] = typer.Argument(..., help="BEAST log files"),
        directory: Path = typer.Argument(..., help="Output directory"),
        burnin: float = typer.Option(0.1, help="Fraction of the chain to discard as burnin"),
):
    """
    Plots the trace from BEAST log files.

    Args:
      logs (List[Path]): A list of paths to the log files.
      directory (Path): The output directory.
      burnin (float): The fraction of the chain to discard as burnin.

    Returns:
      None: Does not return anything.
    """
    df = load_log_files(logs, burnin=burnin)

    # convert to xarray
    xdata = xr.Dataset.from_dataframe(df)
    dataset = az.InferenceData(posterior=xdata)

    output_file(filename=directory / f"{directory.name}-trace.html", title="Static HTML file")

    with no_browser():
        # hack to save the html file without opening it in the browser
        # must set show=True to save the file
        with az.rc_context(rc={'plot.max_subplots': None}):
            az.plot_trace(dataset, backend="bokeh", show=True)

@app.command()
def summary(
        logs: List[Path] = typer.Argument(..., help="BEAST log files"),
        output: Path = typer.Argument(..., help="Output csv file"),
        burnin: float = 0.1,
    ):
    """
    Generates a summary of the BEAST log files and saves it to a csv file.

    Args:
      logs (List[Path]): A list of paths to the log files.
      output (Path): The output csv file.
      burnin (float): The fraction of the chain to discard as burnin.

    Returns:
      None: Does not return anything.
    """

    df = load_log_files(logs, burnin=burnin)

    xdata = xr.Dataset.from_dataframe(df)
    dataset = az.InferenceData(posterior=xdata)
    summary = az.summary(dataset, round_to=6)

    # save the summary to csv
    summary.to_csv(output)


if __name__ == "__main__":
    app()
