from math import log
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import typer

app = typer.Typer()

@app.command()
def calculate_odds(
    logs: List[Path] = typer.Argument(..., help="The path to the log CSV file"),
    output_file: str = typer.Argument(..., help="The path to the output CSV file where the results will be saved"),
    gamma_shape: float = typer.Option(..., help="The shape parameter for the gamma distribution"),
    gamma_scale: float = typer.Option(..., help="The scale parameter for the gamma distribution"),
    burnin: float = typer.Option(0.1, "--burnin", "-b", help="Fraction of trees to discard as burn-in")
):
    # Read the CSV file into a DataFrame df
    dfs = []
    for log_path in logs:
        duplicate = pd.read_csv(log_path, sep="\t", comment="#")
        # discard burn-in
        dfs.append(duplicate[int(burnin * len(duplicate)):])

    df = pd.concat(dfs)

    # Identify columns that end with .rate and are not 'clock.rate'
    rate_columns = [col for col in df.columns if col.endswith(".rate") and col != "clock.rate"]

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
        pos_bg = df["clock.rate"]

        # Calculate the posterior odds (pos_odds)
        pos_p = np.mean(pos_fg > pos_bg)
        pos_odds = pos_p / (1 - pos_p) if pos_p < 1 else float("inf")

        # Calculate the log difference of the odds ratios, if possible
        log_diff = log(pos_odds) - log(p_odds) if np.isfinite(pos_odds) and np.isfinite(p_odds) else float("inf")

        # Store the results
        results.append({
            "Rate Column": rate_column,
            "p_p": p_p,
            "p_odds": p_odds,
            "pos_p": pos_p,
            "pos_odds": pos_odds,
            "bf": log_diff
        })

    # Create a DataFrame from the results dictionary
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)

    typer.echo(f"Calculated odds and log differences saved to {output_file}")

if __name__ == "__main__":
    app()
