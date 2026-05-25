### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(12306)

# Import Packages
library(ggplot2)
library(ape)
library(phytools)

# Create Directory
dir_name <- "core_asv"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Import Data ----------------------------------------------------------------
core_taxonomy <- read.csv('../../feature_table/All_core_ASV_taxonomy.csv', row.names = 1)
metadata_rs <- read.csv('../../metadata/rhizosphere_metadata_merge_info.csv')
tree_file <- read.tree('../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
### Step 1: Relative

for (a in amplicon) {
    k <- ifelse(a == '16S', 'Bacteria', ifelse(a == 'ITS', 'Fungi', 'Protists'))
    core_table <- read.csv(
        paste0('../../feature_table/', a, '/core_feature_table_relative.tsv'), sep = '\t', row.names = 1)
    
    tmp_df <- data.frame(t(core_table), check.names = F)
    raw_df <- merge(tmp_df, metadata_rs[c('FileID', 'TreeID')], by.x = 'row.names', by.y = 'FileID', all.x = T)
    res_df <- raw_df[-1]
    fin_df <- aggregate(res_df[-ncol(res_df)], by = list(TreeID = res_df$TreeID), FUN = mean)
    row.names(fin_df) <- fin_df$TreeID
    fin_df <- fin_df[-1]
    
    treeid_list <- intersect(row.names(fin_df), tree_file$tip.label)
    drop_list <- tree_file$tip.label[-match(treeid_list, tree_file$tip.label)]
    fin_tree <- drop.tip(tree_file, drop_list)
    fin_table <- fin_df[fin_tree$tip.label,]
    fin_table <- fin_table[, colSums(fin_table) > 0]
    
    write.csv(fin_table, file = paste0('../../feature_table/', a, '/core_feature_table_relative_mean.csv'))
    
    # 为当前物种创建临时存储
    species_res <- data.frame()
    
    for (f in names(fin_table)) {
      print(f)
      res <- try(phylosig(fin_tree, fin_table[, f], method = 'lambda', test = T), silent = TRUE)
      
      if (!inherits(res, "try-error")) {
        tmp_res <- data.frame(
          'KOID' = f,
          'Kingdom' = k,
          'lambda' = res$lambda,
          'logL' = res$logL,
          'logL0' = res$logL0,
          'P' = res$P
        )
        species_res <- rbind(species_res, tmp_res)
      }
    }
    
    # 对当前物种的结果进行FDR校正
    if (nrow(species_res) > 0) {
      species_res$padj <- p.adjust(species_res$P, method = 'fdr')
      
      # 保存当前物种的结果
      file_name <- paste0(dir_name, '/', k, '_relative_core_ASV_phylogenetic_signal.csv')
      write.csv(species_res, file_name, row.names = F)
    }
}


### Step 2: log(Absolute + 1)

for (a in amplicon) {
  k <- ifelse(a == '16S', 'Bacteria', ifelse(a == 'ITS', 'Fungi', 'Protists'))
    core_table <- read.csv(
        paste0('../../feature_table/', a, '/core_feature_table_absolute.tsv'), sep = '\t', row.names = 1)
    
    tmp_df <- data.frame(t(core_table), check.names = F)
    raw_df <- merge(tmp_df, metadata_rs[c('FileID', 'TreeID')], by.x = 'row.names', by.y = 'FileID', all.x = T)
    res_df <- raw_df[-1]
    fin_df <- aggregate(res_df[-ncol(res_df)], by = list(TreeID = res_df$TreeID), FUN = mean)
    row.names(fin_df) <- fin_df$TreeID
    fin_df <- fin_df[-1]
    
    tree_list <- intersect(row.names(fin_df), tree_file$tip.label)
    drop_list <- tree_file$tip.label[-match(tree_list, tree_file$tip.label)]
    fin_tree <- drop.tip(tree_file, drop_list)
    fin_table <- fin_df[fin_tree$tip.label,]
    
    fin_table <- fin_table[, colSums(fin_table) > 0]
    # write.csv(fin_table, file = paste0('../../feature_table/', a, '/core_feature_table_absolute_mean.csv'))
    
    fin_table <- log(fin_table + 1)   ###不能先log再取均值，之前的计算有点问题
    write.csv(fin_table, file = paste0('../../feature_table/', a, '/core_feature_table_absolute_mean_log.csv'))
    
    # 为当前物种创建临时存储
    species_res <- data.frame()
    
    for (f in names(fin_table)) {
      print(f)
      res <- try(phylosig(fin_tree, fin_table[, f], method = 'lambda', test = T), silent = TRUE)
      
      if (!inherits(res, "try-error")) {
        tmp_res <- data.frame(
          'KOID' = f,
          'Kingdom' = k,
          'lambda' = res$lambda,
          'logL' = res$logL,
          'logL0' = res$logL0,
          'P' = res$P
        )
        species_res <- rbind(species_res, tmp_res)
      }
    }
    
    # 对当前物种的结果进行FDR校正
    if (nrow(species_res) > 0) {
      species_res$padj <- p.adjust(species_res$P, method = 'fdr')
      
      # 保存当前物种的结果
      file_name <- paste0(dir_name, '/', k, '_absolute_format_core_ASV_phylogenetic_signal.csv')
      write.csv(species_res, file_name, row.names = F)
    }
}
