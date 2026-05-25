### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "03-get_top10order_network_property"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(igraph)
library(ggClusterNet)
library(dplyr)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
network_color <- read.csv("../../01-sort_data/08-network_color/top_5percent_taxa.csv", header = T)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (ord in names(top10order_list)) {
    edge <- read.csv(paste0("02-get_top10order_network/", ord, "_all_edge_property_raw.csv"), header = T)
    node <- read.csv(paste0("02-get_top10order_network/", ord, "_all_node_property_raw.csv"), header = T)
    r_table <- read.csv(paste0("02-get_top10order_network/", ord, "_WGCNA_spearman_r_table.csv"), header = T, row.names = 1)
    core_table <- read.csv(paste0("02-get_top10order_network/", ord, "_core_table.csv"), row.names = 1)
    
    
    ### 1. All network
    tmp_r_table <- r_table[node$name, node$name]
    otu_table <- core_table[names(tmp_r_table),]
    tax_table <- core_taxonomy[names(tmp_r_table),]
    
    write.csv(tmp_r_table, paste0(dir_name, '/', ord, '_all_r_table.csv'))
    write.csv(otu_table, paste0(dir_name, '/', ord, '_all_otu_table.csv'))
    write.csv(tax_table, paste0(dir_name, '/', ord, '_all_tax_table.csv'))
    
    igraph <- graph_from_data_frame(edge, directed = F, vertices = node)
    
    # 计算网络特征
    net_property <- net_properties(igraph)
    write.csv(net_property, paste0(dir_name, '/', ord, '_all_netwrok_property.csv'))

    node_property <- node_properties(igraph)
    zipi <- ZiPiPlot(igraph = igraph, method = 'cluster_fast_greedy')
    ggsave(paste0(dir_name, '/', ord, '_all_network_zipi.pdf'), zipi[[1]], width = 8, height = 6)

    hub_info <- data.frame(node_property,zipi[[2]][row.names(node_property),])
    write.csv(hub_info, paste0(dir_name, '/', ord, '_all_network_hub_info.csv'))
    
    ### 2. InterKingdom network
    # 只关注跨界，去除界内相关数据
    tmp_edge <- edge[gsub('ASV.*', '', edge$node1) != gsub('ASV.*', '', edge$node2),]
    write.csv(tmp_edge, paste0(dir_name, '/', ord, '_inter_edge_property_raw.csv'), row.names = F)
    
    tmp_node <- node[node$name %in% c(tmp_edge$node1, tmp_edge$node2),]
    write.csv(tmp_node, paste0(dir_name, '/', ord, '_inter_node_property_raw.csv'), row.names = F)
    
    # 整理相关性表格
    tmp_r_table <- r_table[tmp_node$name, tmp_node$name]
    tmp_r_table <- tmp_r_table * 0
    for (i in 1:nrow(tmp_edge)) {
        tmp_r_table[tmp_edge$node1[i], tmp_edge$node2[i]] <- tmp_edge$weight[i]
        tmp_r_table[tmp_edge$node2[i], tmp_edge$node1[i]] <- tmp_edge$weight[i]
    }
    for (i in names(tmp_r_table)) {tmp_r_table[i, i] <- 1}
    write.csv(tmp_r_table, paste0(dir_name, '/', ord, '_inter_r_table.csv'))
    
    otu_table <- core_table[names(tmp_r_table),]
    tax_table <- core_taxonomy[names(tmp_r_table),]
    
    write.csv(otu_table, paste0(dir_name, '/', ord, '_inter_otu_table.csv'))
    write.csv(tax_table, paste0(dir_name, '/', ord, '_inter_tax_table.csv'))
    
    igraph <- graph_from_data_frame(tmp_edge, directed = F, vertices = tmp_node)
    
    # 计算网络特征
    net_property <- net_properties(igraph)
    write.csv(net_property, paste0(dir_name, '/', ord, '_inter_netwrok_property.csv'))

    node_property <- node_properties(igraph)
    zipi <- ZiPiPlot(igraph = igraph, method = 'cluster_fast_greedy')
    ggsave(paste0(dir_name, '/', ord, '_inter_network_zipi.pdf'), zipi[[1]], width = 8, height = 6)

    hub_info <- data.frame(node_property,zipi[[2]][row.names(node_property),])
    write.csv(hub_info, paste0(dir_name, '/', ord, '_inter_network_hub_info.csv'))
    
    ### 3. IntraKingdom network
    # 只关注界内，去除跨界相关数据
    tmp_edge <- edge[gsub('ASV.*', '', edge$node1) == gsub('ASV.*', '', edge$node2),]
    write.csv(tmp_edge, paste0(dir_name, '/', ord, '_intra_edge_property_raw.csv'), row.names = F)
    
    tmp_node <- node[node$name %in% c(tmp_edge$node1, tmp_edge$node2),]
    write.csv(tmp_node, paste0(dir_name, '/', ord, '_intra_node_property_raw.csv'), row.names = F)
    
    # 整理相关性表格
    tmp_r_table <- r_table[tmp_node$name, tmp_node$name]
    tmp_r_table <- tmp_r_table * 0
    for (i in 1:nrow(tmp_edge)) {
        tmp_r_table[tmp_edge$node1[i], tmp_edge$node2[i]] <- tmp_edge$weight[i]
        tmp_r_table[tmp_edge$node2[i], tmp_edge$node1[i]] <- tmp_edge$weight[i]
    }
    for (i in names(tmp_r_table)) {tmp_r_table[i, i] <- 1}
    write.csv(tmp_r_table, paste0(dir_name, '/', ord, '_intra_r_table.csv'))
    
    otu_table <- core_table[names(tmp_r_table),]
    tax_table <- core_taxonomy[names(tmp_r_table),]
    
    write.csv(otu_table, paste0(dir_name, '/', ord, '_intra_otu_table.csv'))
    write.csv(tax_table, paste0(dir_name, '/', ord, '_intra_tax_table.csv'))
    
    igraph <- graph_from_data_frame(tmp_edge, directed = F, vertices = tmp_node)
    
    # 计算网络特征
    net_property <- net_properties(igraph)
    write.csv(net_property, paste0(dir_name, '/', ord, '_intra_netwrok_property.csv'))

    node_property <- node_properties(igraph)
    zipi <- ZiPiPlot(igraph = igraph, method = 'cluster_fast_greedy')
    ggsave(paste0(dir_name, '/', ord, '_intra_network_zipi.pdf'), zipi[[1]], width = 8, height = 6)

    hub_info <- data.frame(node_property,zipi[[2]][row.names(node_property),])
    write.csv(hub_info, paste0(dir_name, '/', ord, '_intra_network_hub_info.csv'))
}

# ------------------------------------------------------------------------------
