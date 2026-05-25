### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "02-absolute_table"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
kin <- "16S"
path_prefix <- paste0("../00-rawdata/", kin, "/")

# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
feature_table <- read.delim(paste0(path_prefix, '/table/fin_feature_table.tsv'), check.names = F)
copies_df <- read.csv('01-absolute_abundance/all_sample_copies.csv', row.names = 1)

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
absolute_table <- feature_table[c('FeatureID', row.names(copies_df))]

for (i in names(absolute_table)[-1]) {
    absolute_table[i] <- round(absolute_table[i] / copies_df[i, 'Bacteria'] * copies_df[i, 'Copies'])
}

name <- paste0(dir_name, "/absolute_table")
write.table(absolute_table, paste0(name, ".tsv"), sep = '\t', row.names = F, quote = F)
# ------------------------------------------------------------------------------

