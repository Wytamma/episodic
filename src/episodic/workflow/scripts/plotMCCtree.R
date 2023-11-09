#!/usr/bin/env Rscript

# Load the necessary libraries
library(ggplot2)
library(ggtree)
library(treeio)
library(optparse)

# Set up the command line argument parser
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the BEAST MCC tree file", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="tree_plot.png",
              help="Output file for the tree plot (e.g., tree_plot.png)", metavar="file"),
  make_option(c("-m", "--mrsd"), type="character", default=NULL,
              help="Most recent sampling date (MRSD) for the tree", metavar="date")
)

# Parse command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if the file argument is provided
if (is.null(opt$input)) {
  stop("No input file provided. Use --input to specify the BEAST MCC tree file.", call.=FALSE)
}

# Read in the BEAST tree output
beast <- read.beast(opt$input)



# Create the plot
p <- ggtree(beast, aes(color=rate), size=0.7, mrsd = opt$mrsd) +
    geom_range(range='height_0.95_HPD', color='#3498db', alpha=.6, size=2) +
    geom_nodelab(aes(label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_color_continuous(type = "viridis") +
    scale_size_continuous(range=c(1, 3)) +
    theme_tree2() +
    theme(legend.position=c(.1, .8)) 

if(is.null(opt$mrsd)) {
    p <- revts(p) + scale_x_continuous(labels=abs)
}

# Save the plot to the specified output file
ggsave(opt$output, plot=p)

