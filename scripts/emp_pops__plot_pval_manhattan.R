#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option("--in_file",
                        type = "character",
                        default = "",
                        help = "Base dir."
  ),
  
  optparse::make_option("--gene_loc",
                        type = "character",
                        default = "",
                        help = "Base dir."
  ),
  
  optparse::make_option("--output_file",
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
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(qvalue))
################################################################################

################################ Functions #####################################
################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
theme_size <- 60

MIN_PVALUE <- 1 / 10001
df <- read.csv(arguments$options$in_file, header = T, sep='\t')
df$qval <- p.adjust(df$pvalue, method="BH")

print("Number of genes with qval < 0.05:")
print(sum(df$qval < 0.05))

print(paste0("Number of genes with pval:", MIN_PVALUE))
print(sum(df$pvalue <= MIN_PVALUE))

# Add chromosome location to genes
gene_info <- read.csv(
  arguments$options$gene_loc,
  header = T,
  sep = '\t',
  row.names = 'ENSGID'
)
df$chr <- gene_info[df$ENSGID, 'CHR']
df$pos <- (gene_info[df$ENSGID, 'START'] + gene_info[df$ENSGID, 'END']) / 2
df$index <- df$symbol

don <- df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarise(chr_len=max(pos)) %>% # get size
  dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>% 
  dplyr::select(-chr_len) %>%
  dplyr::left_join(df, ., by=c("chr"="chr")) %>%
  dplyr::arrange(chr, pos) %>%
  dplyr::mutate(pos_cum=pos+tot) %>% # cumulative pos
  dplyr::mutate(
    is_highlight = ifelse(qval <= 0.05, "yes", "no"),
    is_annotate = ifelse(qval <= 0.05, "yes", "no")
  )


# Prepare X axis
axisdf <- don %>%
  dplyr::group_by(chr) %>%
  dplyr::summarize(center=(max(pos_cum) + min(pos_cum))/2)

options(ggrepel.max.overlaps = Inf)

# Plot now
p <- ggplot2::ggplot(don, ggplot2::aes(x=pos_cum, y=-log10(pvalue))) +
  ggplot2::geom_point(
    data = subset(don, is_highlight=="no"),
    mapping = ggplot2::aes(color=as.factor(chr)),
    alpha=0.8,
    size=1.3
  ) +
  ggplot2::geom_point(
    data=subset(don, is_highlight=="yes"),
    colour = "black",
    size = 2,
    stroke = 0
  ) +
  ggplot2::scale_color_manual(
    values = rep(c("grey", "darkgrey"), 22 ),
    guide='none'
  ) +
  ggplot2::scale_fill_brewer(palette='Dark2', direction=-1) +
  ggrepel::geom_label_repel(
    data=subset(don, is_annotate=="yes"),
    mapping = ggplot2::aes(label=index),
    size=2
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = 'Gene Position', y='-log10(pvalue)', color='FDR < 5%?') +
  ggplot2::scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  ggplot2::scale_y_continuous(expand = c(0, 0) ) +
  ggplot2::theme( 
    legend.position = 'bottom',
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ggplot2::scale_shape_manual(
    values=c("above:yes" = 24, "inside:yes" = 21, "below:no" = 25,
             "above:no" = 2, "inside:no" = 19, "below:no" = 6),
    guide = "none"
  ) +
  ggplot2::guides(fill = guide_legend(override.aes=list(shape = 21), nrow=2)) + 
  ggplot2::expand_limits(y = c(0, -log10(min(don$pvalue)) * 1.05)) + 
  ggplot2::ggtitle(paste0("Manhattan Plot for trait: ", basename(arguments$options$in_file)))
    
ggplot2::ggsave(file.path(arguments$options$output_file), p, height=11.75, width=9.5)
