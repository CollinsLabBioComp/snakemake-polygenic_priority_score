#!/usr/bin/env Rscript

# disable stringsAsFactors and re-enable upon exit
old <- options(stringsAsFactors = FALSE)
on.exit(options(old), add = TRUE)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option("--observed",
                        type = "character",
                        default = "",
                        help = "Observed input"
  ),
  
  optparse::make_option("--expected",
                        type = "character",
                        default = "",
                        help = "Expected input"
  ),
  
  optparse::make_option("--score_col",
                        type = "character",
                        default = "",
                        help = "Score column"
  ),
  
  optparse::make_option("--output",
                        type = "character",
                        default = "",
                        help = "Out file"
  ),
  
  optparse::make_option(c("-v", "--verbose"),
                        action = "store_true",
                        default = FALSE,
                        help = ""
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Plots results from differential gene expression."
  )
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
  optionStrings <- character()
  for (item in parserObj@options) {
    optionStrings <- append(optionStrings,
                            c(item@short_flag, item@long_flag))
  }
  optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
################################################################################

######################## Required Packages #####################################
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
################################################################################

######################## Read Data & Manipulate ################################
out_base <- arguments$options$output

# Read in data
observed <- read.csv(arguments$options$observed, header = T, sep = '\t')
expected <- read.csv(arguments$options$expected, header = T, sep = '\t')

# label and score
observed$set <- 'observed'
observed$score <- observed[[arguments$options$score_col]]
expected$set <- 'permuted'
expected$score <- expected[[arguments$options$score_col]]

cols <- intersect(
  colnames(observed),
  colnames(expected)
)
df <- rbind(observed[cols], expected[cols])

# Make a qqplot of the -log10(p-values)
p <- ggplot2::ggplot(df, ggplot2::aes(sample = score, color = set)) +
  ggplot2::geom_qq(size = 0.5) +
  ggplot2::scale_color_brewer(palette = "Dark2")

ggsave(paste0(out_base, "_qq.png"), p, width = 6, height = 4)


max_score <- max(df$score)
p95 <- quantile(expected$score, 0.95)
p99 <- quantile(expected$score, 0.99)
p999 <- quantile(expected$score, 0.999)

# Same data, but run a density plot
p <- ggplot2::ggplot(df, ggplot2::aes(x = score, fill = set)) +
  ggplot2::geom_density(alpha = 0.25) +
  ggplot2::scale_fill_manual(values=c('observed'="yellow", 'permuted'="blue")) +
  ggplot2::lims(x = c(0, max_score)) +
  ggplot2::theme_bw() +
  ggplot2::geom_vline(xintercept = p95, linetype="dashed") +
  ggplot2::geom_vline(xintercept = p99, linetype="dashed") +
  ggplot2::geom_vline(xintercept = p999, linetype="dashed") +
  ggplot2::geom_vline(xintercept = max(expected$score), linetype="dashed")

ggsave(paste0(out_base, "_density.png"), p, width = 6, height = 4)


# Make a second plot, but zoomed in 
p <- ggplot2::ggplot(df, ggplot2::aes(x = score, fill = set)) +
  ggplot2::geom_density(alpha = 0.25) +
  ggplot2::theme_bw() +
  ggplot2::lims(x = c(0, max_score)) +
  ggplot2::scale_fill_manual(values=c('observed'="yellow", 'permuted'="blue")) +
  ggplot2::geom_vline(xintercept = p95, linetype="dashed") +
  ggplot2::geom_vline(xintercept = p99, linetype="dashed") +
  ggplot2::geom_vline(xintercept = p999, linetype="dashed") +
  ggplot2::geom_vline(xintercept = max(expected$score), linetype="dashed") +
  ggplot2::scale_x_continuous(expand = expansion(mult = c(0, 0)))

stats <- ggplot2::ggplot_build(p)
p95_max <- stats$data[[1]] %>% filter(x >= p95)
p95_max <- max(p95_max$y)

p_zoomed_95th <- p + ggplot2::coord_cartesian(xlim = c(p95, max(df$score)), ylim = c(0, p95_max))
ggsave(paste0(out_base, "_density_zoom_95th.png"), p_zoomed_95th, width = 6, height = 4)


p99_max <- stats$data[[1]] %>% filter(x >= p99)
p99_max <- max(p99_max$y)
p_zoomed_99th <- p + ggplot2::coord_cartesian(xlim = c(p99, max(df$score)), ylim = c(0, p99_max))
ggsave(paste0(out_base, "_density_zoom_99th.png"), p_zoomed_99th, width = 6, height = 4)


p999_max <- stats$data[[1]] %>% filter(x >= p999)
p999_max <- max(p999_max$y)
p_zoomed_999th <- p + ggplot2::coord_cartesian(xlim = c(p999, max(df$score)), ylim = c(0, p999_max))
ggsave(paste0(out_base, "_density_zoom_999th.png"), p_zoomed_999th, width = 6, height = 4)


# Make a third plot for the cumulative density function 
p <- ggplot2::ggplot(df, aes(x = score, color = set)) +
  ggplot2::stat_ecdf() +
  ggplot2::scale_color_manual(values=c('observed'="yellow", 'permuted'="blue")) +
  ggplot2::theme_bw()

ggsave(paste0(out_base, "_ecdf.png"), p, width = 6, height = 4)
