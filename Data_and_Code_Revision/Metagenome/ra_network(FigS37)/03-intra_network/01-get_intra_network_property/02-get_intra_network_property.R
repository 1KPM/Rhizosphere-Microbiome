### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "02-get_intra_network_property"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(igraph)
library(ggClusterNet)
library(dplyr)
library(MetaNet)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
type <- "intra"
# ------------------------------------------------------------------------------

### Import data ----------------------------------------------------------------
edge <- read.csv("intra_network_edge.csv", header = T)
node <- read.csv("intra_network_node.csv", header = T)
r_table <- read.csv("intra_r_table.csv", header = T, row.names = 1)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
igraph <- graph_from_data_frame(edge, directed = F, vertices = node)
node_property <- node_properties(igraph)
r_matrix <- as.matrix(r_table)
maptree <- model_maptree2(cor = r_matrix, method = 'cluster_fast_greedy')
maptree_node <- maptree[[1]]
maptree_edge <- edgeBuild(cor = r_matrix, node = maptree_node)
names(maptree_edge)[names(maptree_edge) == 'cor'] <- 'Correlation'
write.csv(maptree_edge,"Core0.2_intra_edge_mt2.csv",row.names = F)

maptree_module <- maptree[[2]]
write.csv(maptree_module, "Core0.2_intra_module_mt2.csv", quote = F, row.names = F)

node_property <- data.frame(node_property)
names(node_property) <- c('Degree', 'Closeness', 'Betweenness', 'CenDegree')
maptree_node <- merge(maptree_node, node_property, by.x = 'row.names', by.y = 'row.names')
write.csv(maptree_node,file = "Core0.2_intra_node_mt2.csv")

all_hub <- data.frame(HubScore = hub_score(igraph)$vector)

all_igraph_m <- module_detect(igraph, method = "cluster_fast_greedy")
all_igraph_m_zp <- zp_analyse(all_igraph_m)
zp_result_df <- get_v(all_igraph_m_zp)

#p_zipi <- zp_plot(all_igraph_m_zp)
#ggsave('intra_network_zipi.pdf', p_zipi, width = 8, height = 6)
all_hub_info <- cbind(zp_result_df, all_hub)
write.csv(all_hub_info, 'intra_network_hub_info.csv')