### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "02-taxonomy"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)


### Define variable -----------------------------------------------------------
tax_16s <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_its <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_pro <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
header <- c("FeatureID", "ASVID", "Number", "Taxon", "Kingdom", "Supergroup", "Phylum", "Class", 
            "Order", "Family", "Genus", "Species", "Clade", "Label", "ASVLabel")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
# 1. 加载Feature数量数据
feature_16s <- read.csv('../00-rawdata/taxonomy/16S/feature-frequency-detail.csv', header = F)
feature_its <- read.csv('../00-rawdata/taxonomy/ITS/feature-frequency-detail.csv', header = F)
feature_pro <- read.csv('../00-rawdata/taxonomy/Protist/feature-frequency-detail.csv', header = F)

# 2. 加载Feature注释数据
taxonomy_16s <- read.csv("../00-rawdata/taxonomy/16S/taxonomy.tsv", row.names = 1, sep = "\t")
taxonomy_its <- read.csv("../00-rawdata/taxonomy/ITS/taxonomy.tsv", row.names = 1, sep = "\t")
taxonomy_pro <- read.csv("../00-rawdata/taxonomy/Protist/taxonomy.tsv", row.names = 1, sep = "\t")

# 3. 加载Core ASV数据
core_table_16s <- read.csv("../00-rawdata/feature_table/16S/core_feature_table.tsv", sep = "\t", row.names = 1)
core_table_its <- read.csv("../00-rawdata/feature_table/ITS/core_feature_table.tsv", sep = "\t", row.names = 1)
core_table_pro <- read.csv("../00-rawdata/feature_table/Protist/core_feature_table.tsv", sep = "\t", row.names = 1)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
names(feature_16s) <- c("FeatureID", "Number")
feature_16s$ASVID <- paste0("bASV", row.names(feature_16s))
asv_tax_16s <- merge(feature_16s, taxonomy_16s, by.x = "FeatureID", by.y = "row.names")
asv_tax_16s$Taxonomy <- asv_tax_16s$Taxon
asv_tax_16s <- separate(asv_tax_16s, Taxonomy, into = tax_16s, sep = ";")
asv_tax_16s[tax_16s] <- lapply(asv_tax_16s[tax_16s], function(x) {gsub(".*__", "", x)})
asv_tax_16s <- asv_tax_16s[order(as.numeric(sub("bASV", "", asv_tax_16s$ASVID))),]
name <- paste0(dir_name, "/16S_ASV_taxonomy.csv")
write.csv(asv_tax_16s, name, row.names = F, quote = F, na = "")

core_tax_16s <- asv_tax_16s[asv_tax_16s$FeatureID %in% row.names(core_table_16s),]
core_tax_16s$Clade <- "Bacteria"
core_tax_16s$Supergroup <- "Bacteria"
core_tax_16s[core_tax_16s == "unidentified" | core_tax_16s == "uncultured"] <- NA
core_tax_16s$Label <- 
    ifelse(!is.na(core_tax_16s$Species), core_tax_16s$Species,
           ifelse(!is.na(core_tax_16s$Genus), core_tax_16s$Genus, 
                  ifelse(!is.na(core_tax_16s$Family), core_tax_16s$Family, 
                         ifelse(!is.na(core_tax_16s$Order), core_tax_16s$Order,
                                ifelse(!is.na(core_tax_16s$Class), core_tax_16s$Class, 
                                       ifelse(!is.na(core_tax_16s$Phylum), core_tax_16s$Phylum, "Bacteria"))))))
core_tax_16s$ASVLabel <- paste0(core_tax_16s$ASVID, "|", core_tax_16s$Label)
core_tax_16s <- core_tax_16s[header]


# ITS
names(feature_its) <- c("FeatureID", "Number")
feature_its$ASVID <- paste0("fASV", row.names(feature_its))
asv_tax_its <- merge(feature_its, taxonomy_its, by.x = "FeatureID", by.y = "row.names")
asv_tax_its$Taxonomy <- asv_tax_its$Taxon
asv_tax_its <- separate(asv_tax_its, Taxonomy, into = tax_its, sep = ";")
asv_tax_its[tax_its] <- lapply(asv_tax_its[tax_its], function(x) {gsub(".*__", "", x)})
asv_tax_its <- asv_tax_its[order(as.numeric(sub("fASV", "", asv_tax_its$ASVID))),]
name <- paste0(dir_name, "/ITS_ASV_taxonomy.csv")
write.csv(asv_tax_its, name, row.names = F, quote = F, na = "")

core_tax_its <- asv_tax_its[asv_tax_its$FeatureID %in% row.names(core_table_its),]
core_tax_its$Clade <- "Fungi"
core_tax_its$Supergroup <- "Fungi"
core_tax_its[core_tax_its == "unidentified" | core_tax_its == "uncultured"] <- NA
core_tax_its$Label <- 
    ifelse(!is.na(core_tax_its$Species), core_tax_its$Species,
           ifelse(!is.na(core_tax_its$Genus), core_tax_its$Genus, 
                  ifelse(!is.na(core_tax_its$Family), core_tax_its$Family, 
                         ifelse(!is.na(core_tax_its$Order), core_tax_its$Order,
                                ifelse(!is.na(core_tax_its$Class), core_tax_its$Class, 
                                       ifelse(!is.na(core_tax_its$Phylum), core_tax_its$Phylum, "Fungi"))))))
core_tax_its$ASVLabel <- paste0(core_tax_its$ASVID, "|", core_tax_its$Label)
core_tax_its <- core_tax_its[header]

# Protist
names(feature_pro) <- c("FeatureID", "Number")
feature_pro$ASVID <- paste0("pASV", row.names(feature_pro))
asv_tax_pro <- merge(feature_pro, taxonomy_pro, by.x = "FeatureID", by.y = "row.names")
asv_tax_pro$Taxonomy <- asv_tax_pro$Taxon
asv_tax_pro <- separate(asv_tax_pro, Taxonomy, into = tax_pro, sep = ";")
asv_tax_pro <- asv_tax_pro[order(as.numeric(sub("pASV", "", asv_tax_pro$ASVID))),]

# modify taxonomy name
asv_tax_pro$Order[asv_tax_pro$Order == "ATCC50593-Flamella-WIM80-lineage"] <- "ATCC50593-Flamella-WIM80"
name <- paste0(dir_name, "/Protist_ASV_taxonomy.csv")
write.csv(asv_tax_pro, name, row.names = F, quote = F, na = "")

core_tax_pro <- asv_tax_pro[asv_tax_pro$FeatureID %in% row.names(core_table_pro),]
core_tax_pro$Clade <- "Protist"
core_tax_pro[core_tax_pro == "unidentified" | core_tax_pro == "uncultured"] <- NA
core_tax_pro$Label <- 
    ifelse(!is.na(core_tax_pro$Species), core_tax_pro$Species,
           ifelse(!is.na(core_tax_pro$Genus), core_tax_pro$Genus, 
                  ifelse(!is.na(core_tax_pro$Family), core_tax_pro$Family, 
                         ifelse(!is.na(core_tax_pro$Order), core_tax_pro$Order,
                                ifelse(!is.na(core_tax_pro$Class), core_tax_pro$Class, 
                                       ifelse(!is.na(core_tax_pro$Division), core_tax_pro$Division, 
                                              ifelse(!is.na(core_tax_pro$Supergroup), core_tax_pro$Supergroup, "Protist")))))))
core_tax_pro$ASVLabel <- paste0(core_tax_pro$ASVID, "|", core_tax_pro$Label)
names(core_tax_pro)[8] <- "Phylum"
core_tax_pro <- core_tax_pro[header]

### All core taxonomy
core_taxonomy <- rbind(core_tax_16s, core_tax_its, core_tax_pro)
row.names(core_taxonomy) <- core_taxonomy$ASVID
name <- paste0(dir_name, "/All_core_ASV_taxonomy.csv")
write.csv(core_taxonomy, name, quote = F, na = "")
# ------------------------------------------------------------------------------
