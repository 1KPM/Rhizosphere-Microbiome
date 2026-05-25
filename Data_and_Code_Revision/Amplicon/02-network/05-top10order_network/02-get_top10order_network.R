### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "02-get_top10order_network"
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
metadata_rs <- read.csv('../../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv')
core_table <- read.csv("../../01-sort_data/05-core_table/All_core_feature_table_absolute.csv", row.names = 1)
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)

top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
row.names(metadata_rs) <- metadata_rs$FileID
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (ord in names(top10order_list)) {
    selected_sample <- na.omit(top10order_list[, ord])

    sub_core_table <- core_table[selected_sample]
    sub_core_table <- sub_core_table[rowSums(sub_core_table) > 0,]
    
    sub_core_taxonomy <- core_taxonomy[row.names(sub_core_table),]
    
    sub_metadata <- metadata_rs[selected_sample,]
    
    write.csv(sub_core_table, paste0(dir_name, '/', ord, '_core_table.csv'), quote = F)

    ### 构建microtable网络对象
    dataset <- microtable$new(sample_table = sub_metadata, otu_table = sub_core_table, tax_table = sub_core_taxonomy)
    dataset$tidy_dataset()

    network <- trans_network$new(
        dataset = dataset, cor_method = 'spearman', use_WGCNA_pearson_spearman = T, nThreads = 4)

    r_table <- network$res_cor_p$cor

    network$cal_network(COR_p_thres = 0.01, COR_cut = r_threshold)


    # 计算网络模块
    network$cal_module(method = 'cluster_fast_greedy')

    # 计算node属性
    network$get_node_table()
    node_property <- network$res_node_table
    write.csv(node_property, paste0(dir_name, '/', ord, '_all_node_property_raw.csv'), row.names = F)

    # 计算edge属性
    network$get_edge_table()
    edge_property <- network$res_edge_table
    edge_property$weight <- ifelse(edge_property$label == '+', edge_property$weight, 0-edge_property$weight)
    write.csv(edge_property, paste0(dir_name, '/', ord, '_all_edge_property_raw.csv'), row.names = F)
    # 获取最终相关性table
    fin_r_table <- r_table * 0
    for (i in 1:nrow(edge_property)) {
        fin_r_table[edge_property$node1[i], edge_property$node2[i]] <- edge_property$weight[i]
        fin_r_table[edge_property$node2[i], edge_property$node1[i]] <- edge_property$weight[i]
    }
    for (i in 1:nrow(fin_r_table)) {fin_r_table[i, i] <- 1}
    write.csv(fin_r_table, paste0(dir_name, '/', ord, '_WGCNA_spearman_r_table.csv'))
}
