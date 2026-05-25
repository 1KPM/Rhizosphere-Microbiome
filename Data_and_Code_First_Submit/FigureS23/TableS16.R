### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(1994)

# Import Packages
library(ggplot2)
library(vegan)

# Create Directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
meta <- c("Bacteria", "Fungi", "Protist")
feature <- c("Class", "Order", "Family", "Genus", "Weather", "Habitat", "Moisture", "Environment", "Morphology", 
             "RootDepth", "GrowthStage", "SoilType", "SoilColor", "Location")
# ------------------------------------------------------------------------------


### Import Data ----------------------------------------------------------------
metadata_rs <- read.csv('data/rhizosphere_metadata_merge_info.csv')
metadata_tree <- read.csv('data/tree_metadata_merge_info.csv')
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
permanova_results <- data.frame()
for (amp in amplicon) {
    matrix_df <- read.table(paste0("data/diversity/", amp, "/distance-matrix.tsv"), row.names = 1)

    group_df <- metadata_rs[metadata_rs$FileID %in% rownames(matrix_df),]
    rownames(group_df) <- group_df$FileID

    for (fea in feature) {
        if (fea == "Location") {
            sub_group <- group_df[group_df$Location != 'Wanding Town', fea, drop = F]
        } else {
            sub_group <- group_df[fea]
        }
        names(sub_group) <- "Group"
        tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
        tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
        tmp_dist <- as.dist(tmp_matrix)
        set.seed(1994)
        pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
        variance <- pcoa_sig['Model', 'R2']
        p.val <- pcoa_sig['Model', 'Pr(>F)']
        
        group_name <- ifelse(amp == "Protist", "Amplicon (18S)", paste0("Amplicon (", amp, ")"))
        tmp_reults <- data.frame(Group = group_name, Feature = fea, Variance = variance, Pvalue = p.val)
        permanova_results <- rbind(permanova_results, tmp_reults)
    }
}


for (met in meta) {
    matrix_df <- read.csv(paste0("data/diversity/Meta/", met, "_KO_bray.csv"), row.names = 1)
    
    group_df <- metadata_tree[metadata_tree$TreeID %in% rownames(matrix_df),]
    rownames(group_df) <- group_df$TreeID 
    
    for (fea in feature) {
        if (fea == "Location") {
            sub_group <- group_df[group_df$Location != 'Wanding Town', fea, drop = F]
        } else {
            sub_group <- group_df[fea]
        }
        names(sub_group) <- "Group"
        tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
        tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
        tmp_dist <- as.dist(tmp_matrix)
        set.seed(1994)
        pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
        variance <- pcoa_sig['Model', 'R2']
        p.val <- pcoa_sig['Model', 'Pr(>F)']
        
        group_name <- paste0("Metagenome (", met, ")")
        tmp_reults <- data.frame(Group = group_name, Feature = fea, Variance = variance, Pvalue = p.val)
        permanova_results <- rbind(permanova_results, tmp_reults)
    }
}
name <- paste0(dir_name, "/Supplemental table 16")
write.csv(permanova_results, paste0(name, ".csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------