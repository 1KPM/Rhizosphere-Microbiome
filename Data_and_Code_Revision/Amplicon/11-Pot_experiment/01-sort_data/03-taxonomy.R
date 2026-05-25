### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "03-taxonomy"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

tax_16s <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
feature_16s <- read.csv('../00-rawdata/16S/table/fin_table/data/feature-frequency-detail.csv', header = F)
taxonomy_16s <- read.delim("../00-rawdata/16S/taxonomy/fin_taxonomy/data/taxonomy.tsv", row.names = 1)

# ------------------------------------------------------------------------------



### Get results ----------------------------------------------------------------
names(feature_16s) <- c("FeatureID", "Number")
feature_16s$ASVID <- paste0("bASV", row.names(feature_16s))
asv_tax_16s <- merge(feature_16s, taxonomy_16s, by.x = "FeatureID", by.y = "row.names")
asv_tax_16s$Taxonomy <- asv_tax_16s$Taxon
asv_tax_16s <- separate(asv_tax_16s, Taxonomy, into = tax_16s, sep = ";")
asv_tax_16s[tax_16s] <- lapply(asv_tax_16s[tax_16s], function(x) {gsub(".*__", "", x)})
asv_tax_16s <- asv_tax_16s[order(as.numeric(sub("bASV", "", asv_tax_16s$ASVID))),]

name <- paste0(dir_name, "/16S_ASV_taxonomy")
write.csv(asv_tax_16s, paste0(name, ".csv"), row.names = F, quote = F, na = "")

# ------------------------------------------------------------------------------