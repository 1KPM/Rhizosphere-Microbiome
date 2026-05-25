### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "02-get_network_property"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(igraph)
library(ggClusterNet)
library(dplyr)

# Define function
source('../00-rawdata/scripts/get_network_robustness.R')
# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

treatment <- c("B", "BF", "BFP", "All")

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (tre in treatment) {
    edge <- read.csv(paste0("01-get_network_data/", tre, "_edge_property_raw.csv"), header = T)
    node <- read.csv(paste0("01-get_network_data/", tre, "_node_property_raw.csv"), header = T)
    r_table <- read.csv(paste0("01-get_network_data/", tre, "_WGCNA_spearman_r_table.csv"), header = T, row.names = 1)
    core_table <- read.csv(paste0("01-get_network_data/", tre, "_table.csv"), row.names = 1, check.names = F)
    core_taxonomy <- read.csv(paste0("01-get_network_data/", tre, "_taxonomy.csv"), row.names = 1)
    
    
    tmp_r_table <- r_table[node$name, node$name]
    otu_table <- core_table[names(tmp_r_table),]
    tax_table <- core_taxonomy[names(tmp_r_table),]
    
    write.csv(tmp_r_table, paste0(dir_name, '/', tre, '_r_table.csv'))
    write.csv(otu_table, paste0(dir_name, '/', tre, '_otu_table.csv'))
    write.csv(tax_table, paste0(dir_name, '/', tre, '_tax_table.csv'))

    igraph <- graph_from_data_frame(edge, directed = F, vertices = node)
    
    # 计算网络特征
    net_property <- net_properties(igraph)
    robustness <- get_random_remove_robustness(tmp_r_table, otu_table, table_type = "absolute")
    
    net_property <- as.data.frame(net_property) %>%
        rownames_to_column("Property") %>%
        rename("Value" = "value") %>%
        bind_rows(data.frame(Property = "robustness", Value = robustness$remove_simulation[10, 2])) %>%
        remove_rownames()
    write.csv(net_property, paste0(dir_name, '/', tre, '_netwrok_property.csv'))

    
    node_property <- node_properties(igraph)
    zipi <- ZiPiPlot(igraph = igraph, method = 'cluster_fast_greedy')
    ggsave(paste0(dir_name, '/', tre, '_all_network_zipi.pdf'), zipi[[1]], width = 8, height = 6)
    
    hub_info <- data.frame(node_property, zipi[[2]][row.names(node_property),])
    write.csv(hub_info, paste0(dir_name, '/', tre, '_network_hub_info.csv'))
}

# ------------------------------------------------------------------------------


