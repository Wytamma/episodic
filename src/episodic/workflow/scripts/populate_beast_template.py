#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from jinja2 import StrictUndefined, Template


@dataclass
class Taxon:
    """
    Dataclass representing a taxon.

    Attributes:
      id (str): The id of the taxon.
      sequence (str): The sequence of the taxon.
      date (float): The date of the taxon.
      uncertainty (float): The uncertainty of the taxon's date.
    """
    id: str
    sequence: str
    date: float
    uncertainty: float = 0.0

@dataclass
class Log:
    """
    Dataclass representing a log file.

    Attributes:
      log_every (int): The frequency at which to log.
      file_name (str): The name of the log file.
    """
    log_every: int
    file_name: str

@dataclass
class MLE(Log):
    """
    Dataclass representing a marginal likelihood estimator.

    Attributes:
      log_every (int): The frequency at which to log.
      file_name (str): The name of the log file.
      results_file_name (str): The name of the results file.
      chain_length (int): The length of the MCMC chain.
      path_steps (int): The number of path steps for the MLE.
    """
    results_file_name: str
    chain_length: int
    path_steps: int


def taxa_from_fasta(fasta_path, date_delimiter="|", date_index=-1) -> List[Taxon]:
    """
    Parses a fasta file into a list of Taxon objects.

    Args:
      fasta_path (Path): The path to the fasta file.
      date_delimiter (str): The delimiter for the date in the fasta header.
      date_index (int): The index of the date in the fasta header.

    Returns:
      List[Taxon]: A list of Taxon objects representing the taxa in the fasta file.

    Raises:
      ValueError: If the fasta file is invalid.
    """
    # Read the fasta file
    with open(fasta_path) as fasta_file:
        fasta_lines = fasta_file.readlines()

    # check if valid fasta
    if not fasta_lines[0].startswith(">"):
        raise ValueError("Invalid fasta file.")
    # Parse the fasta file into Taxon objects. Support multi-line sequences.
    taxa = []
    for line in fasta_lines:
        if line.startswith(">"):
            header = line[1:].strip()
            # 1992/1 = 1992 to 1993
            date_with_uncertainty = header.split(date_delimiter)[date_index]
            date, *uncertainty = date_with_uncertainty.split("/")
            if uncertainty:
                uncertainty = float(uncertainty[0])
            else:
                uncertainty = 0.0
            taxa.append(Taxon(id=header, sequence="", date=float(date), uncertainty=uncertainty))
        else:
            taxa[-1].sequence += line.strip()

    return taxa

def populate_beast_template(
        work_dir: Path,
        name: str,
        template_path: Path,
        alignment_path: Path,
        groups: list,
        clock: str,
        rate_gamma_prior_shape: float = 0.5,
        rate_gamma_prior_scale: float = 0.1,
        chain_length: int = 100000000,
        samples: int = 10000,
        mle_chain_length: int = 1000000,
        mle_path_steps: int = 100,
        mle_log_every: int = 10000,
        date_delimiter="|",
        date_index=-1,
        fixed_tree: Optional[Path] = None,
        *, # Force the user to specify the following arguments by name
        trace: bool = True,
        trees: bool = True,
        mle: bool = True,
    ):
    """
    Populates a Beast XML template with an alignment file.

    Args:
      work_dir (Path): The path to the working directory.
      name (str): The name of the output file.
      template_path (Path): The path to the input Beast template file.
      alignment_path (Path): The path to the input alignment file.
      groups (list): A list of groups to include in the analysis.
      clock (str): The clock model to use in the analysis.
      rate_gamma_prior_shape (float): The shape parameter of the gamma prior on the rate.
      rate_gamma_prior_scale (float): The scale parameter of the gamma prior on the rate.
      chain_length (int): The length of the MCMC chain.
      samples (int): The number of samples to draw from the MCMC chain.
      mle_chain_length (int): The length of the MCMC chain for the marginal likelihood estimator.
      mle_path_steps (int): The number of path steps for the marginal likelihood estimator.
      mle_log_every (int): The log every for the marginal likelihood estimator.
      date_delimiter (str): The delimiter for the date in the fasta header.
      date_index (int): The index of the date in the fasta header.
      fixed_tree (Path): The path to the fixed tree file.

    Keyword Args:
      trace (bool): Whether to enable the trace log.
      trees (bool): Whether to enable the trees log.
      mle (bool): Whether to run the marginal likelihood estimator.

    Returns:
      str: The rendered Beast XML template.

    Examples:
      >>> populate_beast_template(
      ...     work_dir=Path("output"),
      ...     name="my_analysis",
      ...     template_path=Path("template.xml"),
      ...     alignment_path=Path("alignment.fasta"),
      ...     groups=["group1", "group2"],
      ...     clock="strict",
      ...     rate_gamma_prior_shape=0.5,
      ...     rate_gamma_prior_scale=0.1,
      ...     chain_length=100000000,
      ...     samples=10000,
      ...     mle_chain_length=1000000,
      ...     mle_path_steps=100,
      ...     mle_log_every=10000,
      ...     date_delimiter="|",
      ...     date_index=-1,
      ...     fixed_tree=Path("fixed_tree.nwk"),
      ...     trace=True,
      ...     trees=True,
      ...     mle=True,
      ... )
      <Rendered Beast XML template>
    """
    # Load the template
    template = Template(template_path.read_text(), undefined=StrictUndefined)

    # Parse the alignment file into Taxon objects
    taxa = taxa_from_fasta(alignment_path, date_delimiter=date_delimiter, date_index=date_index)

    if fixed_tree is not None:
        fixed_tree = fixed_tree.read_text()

    log_every = max(1, chain_length // samples)

    trace_log = None
    if trace:
        trace_log = Log(
            log_every=log_every,
            file_name=f"{name}.log",
        )

    tree_log = None
    if trees:
        tree_log = Log(
            log_every=log_every,
            file_name=f"{name}.trees",
        )

    mle_log = None
    if mle:
        mle_log = MLE(
            log_every = mle_log_every,
            file_name=f"{name}.mle.log",
            results_file_name=work_dir / f"{name}.mle.results.log",
            chain_length=mle_chain_length,
            path_steps=mle_path_steps,
        )

    # Render the template
    rendered_template = template.render(
        taxa=taxa,
        groups=groups,
        clock=clock,
        fixedTree=fixed_tree,
        rateGammaPriorShape=rate_gamma_prior_shape,
        rateGammaPriorScale=rate_gamma_prior_scale,
        chainLength=chain_length,
        screenLogEvery=log_every,
        traceLog=trace_log,
        treeLog=tree_log,
        marginalLikelihoodEstimator=mle_log,
    )

    # Write the rendered template to a file
    return(rendered_template)


if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Populate a Beast XML template with an alignment file.")
    parser.add_argument("template", type=Path, help="Path to the input Beast template file.")
    parser.add_argument("--alignment", type=Path, required=True, help="Path to the input alignment file.")
    parser.add_argument("--date-delimiter", type=str, default="|", help="Delimiter for the date in the fasta header.")
    parser.add_argument("--date-index", type=int, default=-1, help="Index of the date in the fasta header.")
    parser.add_argument("--groups", nargs="+", required=True, help="List of groups to include in the analysis. Space-separated list.")
    parser.add_argument("--clock", type=str, help="Clock model to use in the analysis.", choices=["strict", "relaxed", "flc-stem", "flc-shared-stem", "flc-clade", "flc-shared-clade"])
    parser.add_argument("--rate-gamma-prior-shape", type=float, default=0.5, help="Shape parameter of the gamma prior on the rate.")
    parser.add_argument("--rate-gamma-prior-scale", type=float, default=0.1, help="Scale parameter of the gamma prior on the rate.")
    parser.add_argument("--chain-length", type=int, default=100000000, help="Length of the MCMC chain.")
    parser.add_argument("--samples", type=int, default=10000, help="Number of samples to draw from the MCMC chain.")
    parser.add_argument("--no-trace", action="store_false", help="Whether to disable the trace log.")
    parser.add_argument("--no-trees", action="store_false", help="Whether to disable the trees log.")
    parser.add_argument("--mle", action="store_true", help="Whether to run the marginal likelihood estimator.")
    parser.add_argument("--mle-chain-length", type=int, default=1000000, help="Length of the MCMC chain for the marginal likelihood estimator.")
    parser.add_argument("--mle-path-steps", type=int, default=100, help="Number of path steps for the marginal likelihood estimator.")
    parser.add_argument("--mle-log-every", type=int, default=10000, help="Log every for the marginal likelihood estimator.")
    parser.add_argument("--fixed-tree", type=Path, help="Path to the fixed tree file.")
    parser.add_argument("--output", type=Path, required=True, help="Path to the output Beast XML file.")

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the function to populate the Beast template
    beast_xml = populate_beast_template(
        work_dir=args.output.parent,
        name=args.output.stem,
        template_path=args.template,
        alignment_path=args.alignment,
        date_delimiter=args.date_delimiter,
        date_index=args.date_index,
        groups=args.groups,
        clock=args.clock,
        chain_length=1 if args.mle else args.chain_length,
        samples=args.samples,
        fixed_tree=args.fixed_tree,
        rate_gamma_prior_shape=args.rate_gamma_prior_shape,
        rate_gamma_prior_scale=args.rate_gamma_prior_scale,
        trace=args.no_trace,
        trees=args.no_trees,
        mle=args.mle,
        mle_chain_length=args.mle_chain_length,
        mle_path_steps=args.mle_path_steps,
        mle_log_every=args.mle_log_every,
    )
    # save to file
    args.output.write_text(beast_xml)
