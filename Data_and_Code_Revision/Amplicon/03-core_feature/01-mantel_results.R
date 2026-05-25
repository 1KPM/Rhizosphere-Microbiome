### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "01-mantel_results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(vegan)

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
occ <- 0.2
kingdom <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


all_df <- data.frame()
### Get results ----------------------------------------------------------------
for (kin in kingdom) {
    otu_path <- paste0("../00-rawdata/feature_table/", kin, "/rarefied_feature_table.tsv")
    otu <- read.delim(otu_path, row.names = 1)
    
    bc_matrix_path <- paste0("../00-rawdata/diversity/", kin, "/distance-matrix.tsv")
    bc_matrix_df <- read.delim(bc_matrix_path, row.names = 1)
    

    core_feature_path <- paste0("../00-rawdata/core_feature/", kin, "/core-features-", occ, "00.tsv")
    core_feature <- read.delim(core_feature_path)
    feature_list <- core_feature$Feature.ID
    
    valid_features <- intersect(rownames(otu), feature_list)
    
    common_samples <- intersect(rownames(bc_matrix_df), colnames(otu))
    
    bc_matrix_sorted <- bc_matrix_df[common_samples, common_samples]
    bc_full <- as.dist(bc_matrix_sorted)
    
    comm_subset <- t(otu[valid_features, common_samples, drop=FALSE])
    
    bc_sub <- vegdist(comm_subset, method = "bray")
    
    mantel_res <- mantel(bc_full, bc_sub, permutations = 999)
    
    lm_res <- summary(lm(as.vector(bc_full) ~ as.vector(bc_sub)))
    
    res_df <- data.frame(
        Kingdom = kin,
        occurrence = occ,
        mantel_r = mantel_res$statistic,
        mantel_p = mantel_res$signif,
        linear_r2 = lm_res$r.squared,
        feature_count = length(valid_features)
    )
    
    all_df <- rbind(all_df, res_df)
    
    name <- paste0(dir_name, "/bray_curtis_", kin, "_", occ)
    write.csv(as.data.frame(as.matrix(bc_full)), paste0(name, "_full.csv"), quote = F)
    write.csv(as.data.frame(as.matrix(bc_sub)), paste0(name, "_sub.csv"), quote = F)
}
name <- paste0(dir_name, "/core_feature_mantel_results")
write.csv(all_df, paste0(name, ".csv"), quote = F, row.names = F)

# ------------------------------------------------------------------------------
