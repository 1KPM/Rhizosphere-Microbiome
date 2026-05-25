
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
dir_name <- "core_KO"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
meta <- c("Bacteria", "Fungi", "Protists")
# ------------------------------------------------------------------------------


### Import Data ----------------------------------------------------------------
ko_annotation <- read.csv('../../kegg_abundance/kegg_pathway_for_1KPM.csv', row.names = 1)
metadata_tree <- read.csv('../../kegg_abundance/tree_metadata_merge_info.csv')
tree_file <- read.tree('../../kegg_abundance/tree_metadata_merge_info_final_align_tree.nwk')
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
### Step 1: Relative (TPM)
all_res <- data.frame()

for (a in meta) {
  
  core_table <- read.csv(
        paste0('../../kegg_abundance/relative/', a, '_core0.2_tpm.csv'), row.names = 1, check.names = F)
  
  names(core_table) <- gsub("RS-", "", names(core_table))
  
  fin_df <- data.frame(t(core_table), check.names = F)
  
  treeid_list <- intersect(row.names(fin_df), tree_file$tip.label)
  
  drop_list <- tree_file$tip.label[-match(treeid_list, tree_file$tip.label)]
  
  fin_tree <- drop.tip(tree_file, drop_list)
  
  fin_table <- fin_df[fin_tree$tip.label, ]
  fin_table <- fin_table[, colSums(fin_table) > 0]
  
  # 为当前物种创建临时存储
  species_res <- data.frame()
  
  for (f in names(fin_table)) {
    print(f)
    res <- try(phylosig(fin_tree, fin_table[, f], method = 'lambda', test = T), silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      tmp_res <- data.frame(
        'KOID' = f,
        'Kingdom' = a,
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
    file_name <- paste0(dir_name, '/', a, '_core0.2_tpm_phylogenetic_signal.csv')
    write.csv(species_res, file_name, row.names = F)
  }
}

###############################################################################################
### Step 2: log(Absolute + 1)
all_res_logqa <- data.frame()

for (a in meta) {
  
  core_table <- read.csv(
    paste0('../../kegg_abundance/absolute/', a, '_core0.2_quantitative.csv'), row.names = 1, check.names = F)
  
  names(core_table) <- gsub("RS-", "", names(core_table))
  
  fin_df <- data.frame(t(core_table), check.names = F)
  fin_df <- log(fin_df + 1)
  
  treeid_list <- intersect(row.names(fin_df), tree_file$tip.label)
  
  drop_list <- tree_file$tip.label[-match(treeid_list, tree_file$tip.label)]
  
  fin_tree <- drop.tip(tree_file, drop_list)
  
  fin_table <- fin_df[fin_tree$tip.label,]
  fin_table <- fin_table[, colSums(fin_table) > 0]
  
  # 为当前物种创建临时存储
  species_res <- data.frame()
  
  for (f in names(fin_table)) {
    print(f)
    res <- try(phylosig(fin_tree, fin_table[, f], method = 'lambda', test = T), silent = TRUE)
    
    if (!inherits(res, "try-error")) {
      tmp_res <- data.frame(
        'KOID' = f,
        'Kingdom' = a,
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
    file_name <- paste0(dir_name, '/', a, '_core0.2_qa_log_phylogenetic_signal.csv')
    write.csv(species_res, file_name, row.names = F)
  }
}







