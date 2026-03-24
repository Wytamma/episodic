#!/usr/bin/env Rscript
#SNK --env ggtree

# Load the necessary libraries
library(ggplot2)
library(ggtree)
library(treeio)
library(optparse)
library(ape)

options(ignore.negative.edge = TRUE)

get_descendant_tips <- function(tree, node) {
  if (node <= Ntip(tree)) return(node)
  subtree <- extract.clade(tree, node)
  match(subtree$tip.label, tree$tip.label)
}

decimal_year_to_date <- function(decimal_year) {
  year <- floor(decimal_year)
  year_start <- as.Date(sprintf("%d-01-01", year))
  next_year_start <- as.Date(sprintf("%d-01-01", year + 1))
  days_in_year <- as.numeric(next_year_start - year_start)

  year_start + round((decimal_year - year) * days_in_year)
}

format_decimal_year_dates <- function(x) {
  vapply(
    x,
    function(value) {
      if (!is.finite(value)) {
        return(NA_character_)
      }
      format(decimal_year_to_date(value), "%Y-%m-%d")
    },
    character(1)
  )
}

group_to_clade_rows <- function(phylo_tree, groups_df) {
  split_groups <- split(groups_df$label, groups_df$group)

  do.call(
    rbind,
    lapply(names(split_groups), function(group) {
      tip_labels <- split_groups[[group]]
      tip_labels <- intersect(tip_labels, phylo_tree$tip.label)
      if (length(tip_labels) == 0) {
        return(NULL)
      }

      node <- if (length(tip_labels) == 1) {
        match(tip_labels[[1]], phylo_tree$tip.label)
      } else {
        getMRCA(phylo_tree, tip_labels)
      }

      if (is.null(node) || is.na(node)) {
        return(NULL)
      }

      data.frame(
        group = group,
        node = node,
        stringsAsFactors = FALSE
      )
    })
  )
}

# Set up the command line argument parser
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to the BEAST MCC tree file",
    metavar = "file"
  ),
  make_option(c("-g", "--groups-file"),
    type = "character", default = NULL,
    help = "Optional TSV file with taxon and group columns",
    metavar = "file"
  ),
  make_option(c("-o", "--output-prefix"),
    type = "character", default = "mcc.tree",
    help = "Output file for the tree plot (e.g., tree_plot.png)",
    metavar = "file"
  ),
  make_option(c("-m", "--mrsd"),
    type = "character",
    default = NULL,
    help = "Most recent sampling date (MRSD) for the tree", metavar = "date"
  ),
  make_option(c("-e", "--ext"),
    type = "character",
    default = ".svg",
    help = "Output file extension (e.g., .png, .pdf, .svg)", metavar = "ext"
  )
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the file argument is provided
if (is.null(opt$input)) {
  stop(
    "No input file provided. Use --input to specify the BEAST MCC tree file.",
    call. = FALSE
  )
}

# Read in the BEAST tree output
beast <- read.beast(opt$input)

# Extract the file name without the extension
output_without_ext <- opt$`output-prefix`
plot_title <- tools::file_path_sans_ext(basename(opt$input))

# Extract the file extension
output_extension <- opt$ext

groups_df <- NULL
if (!is.null(opt$`groups-file`)) {
  groups_df <- read.delim(opt$`groups-file`, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (!all(c("taxon", "group") %in% colnames(groups_df))) {
    stop("Groups file must contain 'taxon' and 'group' columns.", call. = FALSE)
  }
  groups_df <- groups_df[!is.na(groups_df$group) & groups_df$group != "", c("taxon", "group")]
  colnames(groups_df)[colnames(groups_df) == "taxon"] <- "label"
}

# Create the plot
p <- ggtree(beast, aes(color = rate), linewidth = 0.7, mrsd = opt$mrsd) +
  scale_color_continuous(type = "viridis") +
  scale_size_continuous(range = c(1, 3)) +
  labs(title = plot_title, x = "Date") +
  theme_tree2(legend.position = c(0.1, 0.8))
  

# Apply time axis formatting
if (is.null(opt$mrsd)) {
  p <- revts(p) + scale_x_continuous(labels = abs) + labs(x = "Time before present")
} else {
  p <- p + scale_x_continuous(labels = format_decimal_year_dates)
}

# Compute tree bounds (after potential revts)
tree_max <- max(p$data$x, na.rm = TRUE)
tree_min <- min(p$data$x, na.rm = TRUE)
tree_span <- tree_max - tree_min
if (!is.finite(tree_span) || tree_span <= 0) tree_span <- 1

# Default right margin (no clade labels)
right_margin_pt <- 5.5

if (!is.null(groups_df) && nrow(groups_df) > 0) {
  group_ranges <- group_to_clade_rows(as.phylo(beast), groups_df)

  if (!is.null(group_ranges) && nrow(group_ranges) > 0) {
    phylo_tree <- as.phylo(beast)
    tree_data <- p$data

    annot_rows <- lapply(seq_len(nrow(group_ranges)), function(i) {
      node_id <- group_ranges$node[i]
      if (node_id <= length(phylo_tree$tip.label)) {
        tip_ys <- tree_data$y[tree_data$node == node_id]
      } else {
        desc_tips <- get_descendant_tips(phylo_tree, node_id)
        tip_ys <- tree_data$y[tree_data$node %in% desc_tips]
      }
      data.frame(
        group = group_ranges$group[i],
        ymin = min(tip_ys),
        ymax = max(tip_ys),
        ymid = mean(range(tip_ys)),
        stringsAsFactors = FALSE
      )
    })
    annot_df <- do.call(rbind, annot_rows)

    bar_x <- tree_max + tree_span * 0.015
    label_x <- tree_max + tree_span * 0.03

    p <- p +
      geom_segment(
        data = annot_df,
        aes(x = bar_x, xend = bar_x, y = ymin, yend = ymax),
        inherit.aes = FALSE, color = "#6B7280", linewidth = 0.7
      ) +
      geom_text(
        data = annot_df,
        aes(x = label_x, y = ymid, label = group),
        inherit.aes = FALSE, hjust = 0, size = 3, color = "#374151"
      )

    max_label_nchar <- max(nchar(annot_df$group))
    right_margin_pt <- max_label_nchar * 7 + 30
  }
}

# Constrain x-axis to tree range; clip = "off" lets labels render in the margin
p <- p +
  coord_cartesian(
    xlim = c(tree_min - tree_span * 0.02, tree_max),
    clip = "off"
  ) +
  theme(plot.margin = margin(5.5, right_margin_pt, 5.5, 5.5, "pt"))

plot_width <- 7 + right_margin_pt / 72

ggsave(paste0(output_without_ext, output_extension), plot = p, width = plot_width)

p <- p +
  geom_range(range = "height_0.95_HPD", color = "#3498db", alpha = 0.6, linewidth = 2)

ggsave(
  paste0(output_without_ext, ".height_0.95_HPD", output_extension),
  plot = p, width = plot_width
)

p <- p +
  geom_nodelab(aes(label = round(posterior, 2)), vjust = -.5, size = 3)

ggsave(
  paste0(output_without_ext, ".posterior", output_extension),
  plot = p, width = plot_width
)
