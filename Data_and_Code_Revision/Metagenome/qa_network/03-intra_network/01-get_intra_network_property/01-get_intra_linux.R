library(igraph)
library(dplyr)
#library(tidyfst)

##输入数据来自单个目总网络
KO_network_cor_all <- read.csv('../network/All_KO_WGCNA_spearman_cor_raw.csv', row.names = 1)
KO_network_node_all <- read.csv('../network/KO_network_node_property_All.csv')
KO_network_edge_all <- read.csv('../network/KO_network_edge_property_All.csv')
KO_subset <- read.csv('../data/All_core0.2.csv', row.names = 1,check.names = F)

# 只关注界内，去除跨界相关数据
within_network_edge <- KO_network_edge_all[gsub('K.*', '', KO_network_edge_all$node1) == gsub('K.*', '', KO_network_edge_all$node2),]
write.csv(within_network_edge, 'Core0.2_within_network_edge.csv', row.names = F)

within_network_node <- KO_network_node_all[KO_network_node_all$name %in% c(within_network_edge$node1, within_network_edge$node2),]
write.csv(within_network_node, 'Core0.2_within_network_node.csv', row.names = F)

# 过滤network相关性数据
within_network_cor <- KO_network_cor_all[within_network_node$name, within_network_node$name]
within_network_cor <- within_network_cor * 0
for (i in 1:nrow(within_network_edge)) {
  within_network_cor[within_network_edge$node1[i], within_network_edge$node2[i]] <- within_network_edge$weight[i]
  within_network_cor[within_network_edge$node2[i], within_network_edge$node1[i]] <- within_network_edge$weight[i]
}
for (i in names(within_network_cor)) {within_network_cor[i, i] <- 1}
write.csv(within_network_cor, 'Core0.2_within_kingdom_network_cor.csv')




within_network_cor <- read.csv('Core0.2_within_kingdom_network_cor.csv',row.names = 1) 
within_network_node <- read.csv('Core0.2_within_network_node.csv')
within_network_edge <- read.csv('Core0.2_within_network_edge.csv')


### 按ASV过滤table
within_igraph <- graph_from_data_frame(within_network_edge, directed = F, vertices = within_network_node)

# 计算网络特征
within_net_property <- net_properties(within_igraph)
write.csv(within_net_property, 'Core0.2_within_kingdom_network_property.csv')

