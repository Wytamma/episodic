from pathlib import Path
import re
from typing import List, Optional

import numpy as np
import pandas as pd
import typer


def safe_log(x):
    """Returns the logarithm of x if x is positive, otherwise returns None."""
    if x > 0:
        return np.log(x)
    else:
        # Handle the case where x is not positive
        return None


def calculate_log_diff(pos_odds, p_odds):
    """
    Calculates the log difference of two odds ratios.

    Args:
      pos_odds (float): The posterior odds ratio.
      p_odds (float): The prior odds ratio.

    Returns:
      float: The log difference of the two odds ratios.

    Examples:
      >>> calculate_log_diff(2, 1)
      0.6931471805599453
    """
    log_pos_odds = safe_log(pos_odds)
    log_p_odds = safe_log(p_odds)

    if log_pos_odds is not None and log_p_odds is not None:
        log_diff = log_pos_odds - log_p_odds
    else:
        # Handle the case where one or both of the odds are not valid for log operation
        log_diff = float("inf")  # Or choose another appropriate response

    return log_diff


app = typer.Typer()


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


@app.command()
def calculate_odds(
    logs: List[Path] = typer.Argument(..., help="The path to the log CSV file"),
    output_file: str = typer.Argument(..., help="The path to the output CSV file where the results will be saved"),
    gamma_shape: float = typer.Option(..., help="The shape parameter for the gamma distribution"),
    gamma_scale: float = typer.Option(..., help="The scale parameter for the gamma distribution"),
    foreground_label: Optional[str] = typer.Option(None, help="Optional foreground/local-rate label prefix."),
    background_label: Optional[str] = typer.Option(None, help="Optional background-rate label prefix."),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of trees to discard as burn-in"),
):
    """
    Calculates the odds and log differences for a given set of data.

    Args:
      logs (List[Path]): The paths to the log CSV files.
      output_file (str): The path to the output CSV file where the results will be saved.
      gamma_shape (float): The shape parameter for the gamma distribution.
      gamma_scale (float): The scale parameter for the gamma distribution.
      burnin (float): Fraction of trees to discard as burn-in.

    Returns:
      None

    Examples:
      >>> calculate_odds(['log1.csv', 'log2.csv'], 'results.csv', 2, 1, 0.1)
      Calculated odds and log differences saved to results.csv
    """
    # Read the CSV file into a DataFrame df
    dfs = []
    for log_path in logs:
        duplicate = pd.read_csv(log_path, sep="\t", comment="#")
        # discard burn-in
        dfs.append(duplicate[int(burnin * len(duplicate)) :])

    df = pd.concat(dfs)

    # Identify local clock rate columns and compare them to the matching
    # partition background rate. For backward compatibility, unpartitioned
    # analyses still use `clock.rate`.
    partition_backgrounds = {
        _rate_key(column, ".clock.rate", background_label or "background"): column
        for column in df.columns
        if column.endswith(".clock.rate")
    }
    default_background = "clock.rate" if "clock.rate" in df.columns else None
    if default_background is None and len(partition_backgrounds) == 1:
        default_background = next(iter(partition_backgrounds.values()))

    rate_columns = []
    for col in df.columns:
        if not col.endswith(".rate"):
            continue
        if col == default_background or col.endswith(".clock.rate"):
            continue
        rate_columns.append(col)

    # Generate gamma-distributed random variables for p_fg and p_bg
    n_samples = len(df)
    p_fg = np.random.gamma(gamma_shape, gamma_scale, n_samples)
    p_bg = np.random.gamma(gamma_shape, gamma_scale, n_samples)

    # Calculate the prior odds (p_odds)
    p_p = np.mean(p_fg > p_bg)
    p_odds = p_p / (1 - p_p) if p_p < 1 else float("inf")

    # Results dictionary
    results = []

    # Calculate the odds for each .rate column compared to df['clock.rate']
    for rate_column in rate_columns:
        pos_fg = df[rate_column]

        rate_key = _rate_key(rate_column, ".rate", foreground_label)
        local_parts = rate_key.split(".")

        background_column = partition_backgrounds.get(rate_key, default_background)
        while background_column is None and local_parts:
            local_parts = local_parts[:-1]
            background_column = partition_backgrounds.get(".".join(local_parts))

        if background_column is None:
            raise ValueError(
                f"Could not find a background clock rate column for '{rate_column}'."
            )
        pos_bg = df[background_column]

        # Calculate the posterior odds (pos_odds)
        pos_p = np.mean(pos_fg > pos_bg)
        pos_odds = pos_p / (1 - pos_p) if pos_p < 1 else float("inf")

        # Calculate the log difference of the odds ratios, if possible
        log_diff = calculate_log_diff(pos_odds, p_odds)

        # Store the results
        results.append(
            {
                "Rate Column": rate_column,
                "Background Column": background_column,
                "p_p": p_p,
                "p_odds": p_odds,
                "pos_p": pos_p,
                "pos_odds": pos_odds,
                "bf": log_diff,
            }
        )

    # Create a DataFrame from the results dictionary
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)

    typer.echo(f"Calculated odds and log differences saved to {output_file}")


if __name__ == "__main__":
    app()
