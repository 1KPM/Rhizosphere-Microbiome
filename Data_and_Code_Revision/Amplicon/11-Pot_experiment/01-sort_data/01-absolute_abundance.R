### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "01-absolute_abundance"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
copies_16s_number <- 42332891

kin <- "16S"
path_prefix <- paste0("../00-rawdata/", kin, "/")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
clean_tax_table <- read.csv(paste0(path_prefix, "taxonomy/clean_bar_plots/data/level-4.csv"), row.names = 1)
fin_tax_table <- read.csv(paste0(path_prefix, "taxonomy/fin_bar_plots/data/level-1.csv"), row.names = 1)
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
spike_df <- clean_tax_table %>%
    select(contains("Spike")) %>%
    rename_with(~ "SynSpike")

bacteria_df <- fin_tax_table %>%
    rename(Bacteria = d__Bacteria)

copies_df <- merge(spike_df, bacteria_df, by = 'row.names')
names(copies_df)[1] <- "SampleID"
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
copies_df$Copies <- round((copies_df$Bacteria / copies_df$SynSpike) * copies_16s_number * (1000 / copies_df$Weight) 
                           * (100 / 30), digits = 0)
name <- paste0(dir_name, "/all_sample_copies")
write.csv(copies_df, paste0(name, ".csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------