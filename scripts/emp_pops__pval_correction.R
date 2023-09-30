#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(0)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
  optparse::make_option("--in_file",
                        type = "character",
                        default = "qtl_results",
                        help = "TSV containing pvalue to correct"
  ),
  
  optparse::make_option("--pvalue_col",
                        type = "character",
                        default = "pvalue",
                        help = "Column to correct over."
  ),
  
  optparse::make_option("--correction_method",
                        type = "character",
                        default = "bh",
                        help = "Method to correct over"
  ),
  
  optparse::make_option("--out_file",
                        type = "character",
                        default = "corrected_pvalue.tsv",
                        help = "Output file."
  )
)

parser <- optparse::OptionParser(
  usage = "%prog",
  option_list = optionList,
  description = paste0(
    "Corrects diffential results using IHW."
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
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(dplyr))
################################################################################

################################ Functions #####################################

################################################################################

######################## Read Data & Manipulate ################################
output_file <- arguments$options$out_file

df <- read.csv(
  arguments$options$in_file,
  sep = "\t",
  header = T
)
pval_col <- arguments$options$pvalue_col
correction_method <- arguments$options$correction_method

# For ease, set pval column
orig_columns <- colnames(df)
df$pval__UNIQUE <- df[[pval_col]]

out_col <- sprintf("%s_%s", correction_method, pval_col)
if (correction_method == 'storey_fdr') {
  df <- df %>%
    dplyr::mutate(
      !!out_col := qvalue(pval__UNIQUE)$qvalues
    )
} else if (correction_method == 'bh') {
  df <- df %>%
    dplyr::mutate(
      !!out_col := p.adjust(pval__UNIQUE, method='BH')
    )
} else if (correction_method == 'bonf') {
  df <- df %>%
    dplyr::mutate(
      !!out_col := p.adjust(pval__UNIQUE, method='bonferroni')
    )
}

orig_columns <- c(orig_columns, out_col)
write.table(
  x = df[orig_columns],
  file = output_file,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

################################################################################