# ******************************************************************************
# @File: 05-sub_network.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-16 16:16:48
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "05-sub_network"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(igraph)
library(ggClusterNet)

# Define function
source('../../00-rawdata/scripts/get_network_robustness.R')

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (typ in type) {
    if (typ == "all") {
        dir_main <- "../01-all_network"
        dir_data <- "01-get_all_network_data"  # edge/node 所在目录
        dir_prop <- "02-get_all_network_property"  # table 所在目录
    } else if (typ == "inter") {
        dir_main <- "../02-inter_network"
        dir_data <- "01-get_inter_network_property"
        dir_prop <- "01-get_inter_network_property"
    } else if (typ == "intra") {
        dir_main <- "../03-intra_network"
        dir_data <- "01-get_intra_network_property"
        dir_prop <- "01-get_intra_network_property"
    }
    
    otu_table_path <- file.path(dir_main, dir_prop, paste0(typ, "_otu_table.csv"))
    r_table_path   <- file.path(dir_main, dir_prop, paste0(typ, "_r_table.csv"))
    tax_table_path <- file.path(dir_main, dir_prop, paste0(typ, "_tax_table.csv"))
    
    edge_path <- file.path(dir_main, dir_data, paste0(typ, "_edge_property_raw.csv"))
    node_path <- file.path(dir_main, dir_data, paste0(typ, "_node_property_raw.csv"))
    
    otu_table <- read.csv(otu_table_path, row.names = 1)
    r_table <- read.csv(r_table_path, row.names = 1)
    tax_table <- read.csv(tax_table_path, row.names = 1)
    edge <- read.csv(edge_path, header = T)
    node <- read.csv(node_path, header = T)
    
    igraph <- graph_from_data_frame(edge, directed = F, vertices = node)

    tmp_data_df <- NULL
    for (sam in names(otu_table)) {
        sub_asv <- row.names(otu_table)[otu_table[sam] != 0]
        sub_asv <- intersect(sub_asv, V(igraph)$name)    
        sub_igraph <- induced_subgraph(igraph, vids = sub_asv)
        net_property <- net_properties(sub_igraph)
        
        tmp_r_table <- r_table[sub_asv, sub_asv]
        tmp_otu_table <- otu_table[sub_asv, ]
        
        robustness <- get_random_remove_robustness(tmp_r_table, tmp_otu_table, table_type = "absolute")
        
        data_df <- as.data.frame(net_property) %>%
            rownames_to_column("Property") %>%
            rename("Value" = "value") %>%
            bind_rows(data.frame(Property = "robustness", Value = robustness$remove_simulation[10, 2])) %>%
            mutate(Sample = sam) %>%
            remove_rownames()
        
        tmp_data_df <- tmp_data_df %>%
            bind_rows(data_df)
    }
    name <- paste0(dir_name, "/", typ, "_network_property")
    write.csv(tmp_data_df, paste0(name, ".csv"), quote = F, row.names = F)
}

# ------------------------------------------------------------------------------