#!/usr/bin/env Rscript
# CLI R script using ggtree to create a ggdensitree from a BEAST trees log

# Load the necessary libraries
library(ggtree)
library(optparse)
library(treeio)  # Explicitly load the treeio library
library(ape)

# Set up the command-line options
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the BEAST tree log file", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = "tree_plot.pdf",
              help = "Output filename for the tree plot", metavar = "file")
)

# Parse the command line arguments
parser = OptionParser(option_list = option_list)
args = parse_args(parser)

# Check if input file is provided
if (is.null(args$input)) {
  stop("No input file provided. Use --input to specify the BEAST tree log file.", call. = FALSE)
}

read_trees_with_ape <- function(file_path) {
  # Read in the BEAST tree log using the ape package
  trees <- read.nexus(file_path)
  
#   # If the BEAST output is a list of trees, convert them to multiPhylo
#   if (is.list(trees)) {
#     trees <- do.call(combine, args = list(trees))
#   }
  
  # Convert trees to ggtree format if needed
  return(trees)
}

# Read the BEAST tree log
btrees <- read_trees_with_ape(args$input)

# Create the ggdensitree
p <- ggdensitree(btrees, alpha = .3, colour = 'steelblue') +
  geom_tiplab(size = 3)

# Save the plot to file
ggsave(args$output, p, width = 8, height = 6)

cat("Plot saved to", args$output, "\n")
