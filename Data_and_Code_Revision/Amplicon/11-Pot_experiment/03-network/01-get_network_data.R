### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "01-get_network_data"
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
treatment <- c("B", "BF", "BFP", "All")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
core_table <- read.csv("../01-sort_data/04-core_table/core_feature_table_absolute.csv", row.names = 1, check.names = F)
core_taxonomy <- read.csv("../01-sort_data/04-core_table/16S_core_taxonomy.csv", row.names = 1)
core_metadata <- read.csv("../01-sort_data/04-core_table/sub_metadata.csv", row.names = 1)

# ------------------------------------------------------------------------------



### Get results ----------------------------------------------------------------
for (tre in treatment) {
    tmp_metadata <- core_metadata %>%
        filter(Inoculation == tre)
    
    tmp_table <- core_table %>%
        select(all_of(row.names(tmp_metadata))) %>%
        filter(rowSums(across(where(is.numeric))) > 0)
    
    tmp_taxonomy <- core_taxonomy %>%
        rownames_to_column("FeatuereID") %>%
        column_to_rownames("ASVID")
        
    tmp_taxonomy <- tmp_taxonomy[row.names(tmp_table),]
    
    name <- paste0(dir_name, "/", tre)
    write.csv(tmp_metadata, paste0(name, "_metadata.csv"), quote = F)
    write.csv(tmp_table, paste0(name, "_table.csv"), quote = F)
    write.csv(tmp_taxonomy, paste0(name, "_taxonomy.csv"), quote = F)
    
    # 构建microtable网络对象
    dataset <- microtable$new(sample_table = tmp_metadata, otu_table = tmp_table, tax_table = tmp_taxonomy)
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
    write.csv(node_property, paste0(dir_name, '/', tre, '_node_property_raw.csv'), row.names = F)
    
    # 计算edge属性
    network$get_edge_table()
    edge_property <- network$res_edge_table
    edge_property$weight <- ifelse(edge_property$label == '+', edge_property$weight, 0-edge_property$weight)
    write.csv(edge_property, paste0(dir_name, '/', tre, '_edge_property_raw.csv'), row.names = F)
    # 获取最终相关性table
    fin_r_table <- r_table * 0
    for (i in 1:nrow(edge_property)) {
        fin_r_table[edge_property$node1[i], edge_property$node2[i]] <- edge_property$weight[i]
        fin_r_table[edge_property$node2[i], edge_property$node1[i]] <- edge_property$weight[i]
    }
    for (i in 1:nrow(fin_r_table)) {fin_r_table[i, i] <- 1}
    write.csv(fin_r_table, paste0(dir_name, '/', tre, '_WGCNA_spearman_r_table.csv'))
    
}


# ------------------------------------------------------------------------------
