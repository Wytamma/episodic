import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import typer
from plotly.subplots import make_subplots


def format_fig(fig):
    """
    Formats the given figure with custom layout and style.

    Args:
      fig (plotly.graph_objects.Figure): The figure to format.

    Returns:
      None.

    Examples:
      >>> fig = go.Figure()
      >>> format_fig(fig)
      None
    """
    fig.update_layout(
        width=1200,
        height=550,
        plot_bgcolor="white",
        title_font_color="black",
        font=dict(
            family="Linux Libertine Display O",
            size=18,
            color="black",
        ),
    )
    gridcolor = "#dddddd"
    fig.update_xaxes(gridcolor=gridcolor)
    fig.update_yaxes(gridcolor=gridcolor)

    fig.update_xaxes(showline=True, linewidth=1, linecolor="black", mirror=True, ticks="outside")
    fig.update_yaxes(showline=True, linewidth=1, linecolor="black", mirror=True, ticks="outside")
    fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor=gridcolor)
    fig.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor=gridcolor)


def camel_to_title_case(column: str):
    """
    Transforms camel case to title case with space delimited.

    Regex from https://stackoverflow.com/a/9283563
    """
    return re.sub(r"((?<=[a-z])[A-Z]|(?<!\A)[A-Z](?=[a-z]))", r" \1", column).title()


def plot_traces(
    trace_log: Path = typer.Argument(..., help="The path to the trace log from beast."),
    output: Path = typer.Argument(..., help="A path to the output directory."),
    burnin: float = 0.1,
):
    """
    Plots the traces from a trace log file.

    Args:
      trace_log (pathlib.Path): The path to the trace log from beast.
      output (pathlib.Path): A path to the output directory.
      burnin (float): The percentage of the trace to discard as burn-in.

    Returns:
      None.

    Examples:
      >>> trace_log = Path('trace.log')
      >>> output = Path('output')
      >>> burnin = 0.1
      >>> plot_traces(trace_log, output, burnin)
      None
    """
    output.mkdir(exist_ok=True, parents=True)

    df = pd.read_csv(trace_log, sep="\t", comment="#")
    burnin_df = df.truncate(after=burnin * len(df))
    posterior_df = df.truncate(before=burnin * len(df))
    posterior_color = "#1A32E6"
    posterior_background = "#DDE1FA"
    for variable in df.columns[1:]:
        print(variable)
        variable_title = camel_to_title_case(variable)
        print(f"Plotting {variable_title}")
        fig = make_subplots(
            rows=1, cols=3, subplot_titles=("Burn-in", "Posterior", "Histogram"), column_widths=[0.2, 0.6, 0.2]
        )
        posterior_max = posterior_df[variable].max()
        posterior_min = posterior_df[variable].min()
        y_range_min = posterior_min - 0.1 * np.abs(posterior_max - posterior_min)
        y_range_max = posterior_max + 0.1 * np.abs(posterior_max - posterior_min)
        fig.add_shape(
            type="rect",
            xref="x domain",
            yref="y domain",
            x0=0,
            y0=0,
            x1=1,
            y1=1,
            fillcolor="red",
            col=1,
            row=1,
            layer="below",
            opacity=0.25,
        )
        fig.add_shape(
            type="rect",
            xref="x domain",
            x0=0,
            y0=y_range_min,
            x1=1,
            y1=y_range_max,
            fillcolor=posterior_background,
            col=1,
            row=1,
            line={"width": 0},
            layer="below",
        )
        # fig.add_shape(
        #     type="rect",
        #     xref="x domain", yref="y domain",
        #     x0=0, y0=0,
        #     x1=1, y1=1.0,
        #     fillcolor=posterior_background,
        #     col=2, row=1,
        #     line={'width':0},
        #     layer="below",
        # )
        # fig.add_shape(
        #     type="rect",
        #     xref="x domain", yref="y domain",
        #     x0=0, y0=0,
        #     x1=1, y1=1.0,
        #     fillcolor=posterior_background,
        #     col=3, row=1,
        #     line={'width':0},
        #     layer="below",
        # )
        fig.add_trace(
            go.Scatter(
                x=burnin_df["state"],
                y=burnin_df[variable],
                marker_color="red",
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=posterior_df["state"], y=posterior_df[variable], marker_color=posterior_color),
            row=1,
            col=2,
        )
        fig.add_trace(
            go.Histogram(
                y=posterior_df[variable], marker_color=posterior_color, histnorm="probability density", opacity=0.5
            ),
            row=1,
            col=3,
        )
        mean = posterior_df[variable].mean()
        fig.add_annotation(
            x=1,
            y=mean,
            xref="x domain",
            yref="y",
            text=f"mean: {mean:.2g}",
            showarrow=False,
            align="right",
            xanchor="right",
            yanchor="bottom",
            row=1,
            col=3,
        )
        fig.add_shape(
            type="line",
            x0=0,
            y0=mean,
            x1=1,
            y1=mean,
            xref="x domain",
            yref="y",
            line=dict(
                color="#777777",
                width=1.5,
                dash="dot",
            ),
            row=1,
            col=3,
        )

        format_fig(fig)
        fig.data[0].update(mode="markers+lines")
        fig.update_layout(
            yaxis_title=variable_title,
            showlegend=False,
        )
        fig.update_layout(
            xaxis1_range=[burnin_df["state"].min(), burnin_df["state"].max()],
            yaxis2_range=[y_range_min, y_range_max],
            yaxis3_range=[y_range_min, y_range_max],
            xaxis1_title="states",
            xaxis2_title="states",
            xaxis3_title="Density",
        )
        print([y_range_min, y_range_max])
        fig.write_image(output / f"{variable}.svg")


if __name__ == "__main__":
    typer.run(plot_traces)
