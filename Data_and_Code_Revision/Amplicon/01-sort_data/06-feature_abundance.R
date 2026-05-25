### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "06-feature_abundance"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Define variable -----------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (amp in amplicon) {
    raw_df <- read.delim(paste0("../00-rawdata/feature_table/", amp, "/fin_feature_table_absolute.tsv"), row.names = 1)
    sum_df <- data.frame(rowSums(raw_df))
    names(sum_df) <- "Abundance"
    
    name <- paste0(dir_name, "/", amp, "_feature_absolute_abundance")
    write.csv(sum_df, paste0(name, ".csv"), quote = F, row.names = T)
    
    sample_df <- data.frame(colSums(raw_df))
    names(sample_df) <- "Abundance"
    
    name <- paste0(dir_name, "/", amp, "_sample_absolute_abundance")
    write.csv(sample_df, paste0(name, ".csv"), quote = F, row.names = T)
    
}
# ------------------------------------------------------------------------------