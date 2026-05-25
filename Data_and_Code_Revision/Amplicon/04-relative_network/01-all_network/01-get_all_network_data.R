### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "01-get_all_network_data"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(microeco)
library(rgexf)
library(servr)
library(igraph)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

r_threshold <- 0.48
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
core_table <- read.csv("../../01-sort_data/05-core_table/All_core_feature_table_relative.csv", row.names = 1)
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
rs_metadata <- read.csv("../../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv")

# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
core_metadata <- rs_metadata %>% 
    filter(FileID %in% names(core_table)) %>%
    column_to_rownames("FileID")

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 构建microtable网络对象
dataset <- microtable$new(sample_table = core_metadata, otu_table = core_table, tax_table = core_taxonomy)
dataset$tidy_dataset()

# 计算Spearman相关性
# The parameter cor_method in trans_network is used to select correlation calculation method.
# Default Pearson or Spearman correlation invoke R base cor.test, a little slow.
# Spearman correlation based on WGCNA package is applied in all the following operations.
network <- trans_network$new(
    dataset = dataset, cor_method = 'spearman', use_WGCNA_pearson_spearman = T, nThreads = 4)

r_table <- network$res_cor_p$cor


network$cal_network(COR_p_thres = 0.01, COR_cut = r_threshold)


# 计算网络模块
network$cal_module(method = 'cluster_fast_greedy')

# 计算node属性
network$get_node_table()
node_property <- network$res_node_table
write.csv(node_property, paste0(dir_name, '/all_node_property_raw.csv'), row.names = F)

# 计算edge属性
network$get_edge_table()
edge_property <- network$res_edge_table
edge_property$weight <- ifelse(edge_property$label == '+', edge_property$weight, 0-edge_property$weight)
write.csv(edge_property, paste0(dir_name, '/all_edge_property_raw.csv'), row.names = F)       

# 获取最终相关性table
fin_r_table <- r_table * 0
for (i in 1:nrow(edge_property)) {
    fin_r_table[edge_property$node1[i], edge_property$node2[i]] <- edge_property$weight[i]
    fin_r_table[edge_property$node2[i], edge_property$node1[i]] <- edge_property$weight[i]
}
for (i in 1:nrow(fin_r_table)) {fin_r_table[i, i] <- 1}
write.csv(fin_r_table, paste0(dir_name, '/WGCNA_spearman_r_table.csv'))
# ------------------------------------------------------------------------------
