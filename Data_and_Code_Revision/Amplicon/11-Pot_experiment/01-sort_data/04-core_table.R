### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "04-core_table"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
kin <- "16S"

threshold <- 0.2
removed_plant <- c("Lespedeza bicolor", "Medicago sativa", "Tephrosia candida", 
                   "Sesbania cannabina", "Indigofera amblyantha", "Robinia pseudoacacia")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
metadata <- read.delim(paste0("../00-rawdata/", kin, "/metadata/metadata.tsv"))
absolute_table <- read.delim("02-absolute_table/absolute_table.tsv", check.names = F, row.names = 1)

asv_tax_16s <- read.csv("03-taxonomy/16S_ASV_taxonomy.csv", header = T)
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
sub_metadata <- metadata %>%
    filter(Order != "" & !Species %in% removed_plant)

name <- paste0(dir_name, "/sub_metadata")
write.csv(sub_metadata, paste0(name, ".csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
sub_table <- absolute_table[sub_metadata$SampleID]

occurrence_table <- 1 * ((sub_table > 0) == 1)
occurrence_df <- data.frame(Occurrence = rowSums(occurrence_table) / ncol(occurrence_table))

core_table <- sub_table[row.names(occurrence_df)[occurrence_df$Occurrence >= threshold],]


core_tax_16s <- asv_tax_16s %>%
    filter(FeatureID %in% row.names(core_table))
name <- paste0(dir_name, "/16S_core_taxonomy")
write.csv(core_tax_16s, paste0(name, ".csv"), row.names = F, quote = F, na = "")


row.names(core_table) <- core_tax_16s[match(row.names(core_table), core_tax_16s$FeatureID), 'ASVID']

name <- paste0(dir_name, "/core_feature_table_absolute")
write.csv(core_table, paste0(name, ".csv"), quote = F)

# ------------------------------------------------------------------------------