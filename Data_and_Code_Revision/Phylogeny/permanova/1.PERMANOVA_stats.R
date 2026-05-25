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
metadata_rs <- read.csv('../metadata/rhizosphere_metadata_merge_info.csv')
metadata_tree <- read.csv('../metadata/tree_metadata_merge_info.csv')
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
permanova_results <- data.frame()
for (amp in amplicon) {
    matrix_df <- read.csv(paste0("../asv_qa_alpha_beta/02-get_absolute_beta_diversity/", 
                                   amp, 
                                   "_beta_diversity_absolute__bray_curtis.csv"), 
                            fileEncoding = "GBK", row.names = 1)
    

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
        set.seed(12306)
        pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
        variance <- pcoa_sig$'R2'
        p.val <- pcoa_sig$`Pr(>F)`
        
        group_name <- ifelse(amp == "Protist", "Amplicon (18S)", paste0("Amplicon (", amp, ")"))
        tmp_reults <- data.frame(Group = group_name, Feature = fea, Variance = variance, Pvalue = p.val)
        permanova_results <- rbind(permanova_results, tmp_reults)
    }
}


for (met in meta) {
    matrix_df <- read.csv(paste0("../ko_qa_alpha_beta/diversity/",
                                 met,
                                 ".quantitative.bray.csv"), 
                          check.names = F, fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
    names(matrix_df) <- gsub("RS-", "", names(matrix_df))
    row.names(matrix_df) <- gsub("RS-", "", row.names(matrix_df))
    
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
        set.seed(12306)
        pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
        variance <- pcoa_sig$'R2'
        p.val <- pcoa_sig$`Pr(>F)`
        
        group_name <- paste0("Metagenome (", met, ")")
        tmp_reults <- data.frame(Group = group_name, Feature = fea, Variance = variance, Pvalue = p.val)
        permanova_results <- rbind(permanova_results, tmp_reults)
    }
}

permanova_results$Distance <- "Bray Curtis"

name <- paste0(dir_name, "/PERMANOVA_stats_table")
write.csv(permanova_results, paste0(name, ".csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------



#### weighted unifrac

permanova_results <- data.frame()
for (amp in amplicon) {
  matrix_df <- read.csv(paste0("../asv_qa_alpha_beta/02-get_absolute_beta_diversity/", 
                               amp, 
                               "_beta_diversity_absolute__weighted.csv"),
                        fileEncoding = "GBK", row.names = 1)
  
  
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
    set.seed(12306)
    pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
    variance <- pcoa_sig$'R2'
    p.val <- pcoa_sig$`Pr(>F)`
    
    group_name <- ifelse(amp == "Protist", "Amplicon (18S)", paste0("Amplicon (", amp, ")"))
    tmp_reults <- data.frame(Group = group_name, Feature = fea, Variance = variance, Pvalue = p.val)
    permanova_results <- rbind(permanova_results, tmp_reults)
  }
}

permanova_results$Distance <- "Weighted Unifrac"

name <- paste0(dir_name, "/PERMANOVA_stats_table_Weighted_Unifrac")
write.csv(permanova_results, paste0(name, ".csv"), quote = F, row.names = F)


