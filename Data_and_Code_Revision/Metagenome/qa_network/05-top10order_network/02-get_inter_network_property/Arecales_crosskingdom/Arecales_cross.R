pwd <- "/public/home/wangxiaolin_lab/group/core0.2ko_network/Arecales_crosskingdom"
setwd(pwd)

library(igraph)
library(ggClusterNet)
library(dplyr)
library(ggplot2)
library(tidyfst)

##输入数据来自单个目总网络
KO_network_cor_all <- read.csv('../Arecales_network/Arecales_KO_WGCNA_spearman_cor_raw.csv', row.names = 1)
KO_network_node_all <- read.csv('../Arecales_network/Arecales_KO_network_node_property.csv')
KO_network_edge_all <- read.csv('../Arecales_network/Arecales_KO_network_edge_property.csv')
KO_subset <- read.csv('../Arecales_network/Arecales_KO_abundance.csv', row.names = 1,check.names = F)


cross_network_edge <- KO_network_edge_all[gsub('K.*', '', KO_network_edge_all$node1) != gsub('K.*', '', KO_network_edge_all$node2),]
write.csv(cross_network_edge, 'cross_network_edge.csv', row.names = F)

cross_network_node <- KO_network_node_all[KO_network_node_all$name %in% c(cross_network_edge$node1, cross_network_edge$node2),]
write.csv(cross_network_node, 'cross_network_node.csv', row.names = F)

cross_network_cor <- KO_network_cor_all[cross_network_node$name, cross_network_node$name]
cross_network_cor[abs(cross_network_cor) < min(cross_network_edge$weight)] <- 0

for(i in names(cross_network_cor)) {
  cross_network_cor[gsub('K.*', '', row.names(cross_network_cor)) == gsub('K.*', '', i), i] <- 0
}
for (i in names(cross_network_cor)) {cross_network_cor[i, i] <- 1}
write.csv(cross_network_cor, 'Cross_kingdom_netwrok_cor.csv')

cross_KO <- KO_subset[row.names(cross_network_cor),]

cross_igraph <- graph_from_data_frame(cross_network_edge, directed = F, vertices = cross_network_node)

cross_net_property <- net_properties(cross_igraph)
write.csv(cross_net_property, 'Cross_kingdom_netwrok_property.csv')

cross_node_property <- node_properties(cross_igraph)
write.csv(cross_node_property, 'Cross_kingdom_node_property.csv')

cross_hub <- data.frame(HubScore = hub_score(cross_igraph)$vector)
cross_hub_top50 <- head(arrange(cross_hub, desc(HubScore)), 50)

cross_hub_top50_p <- ggplot(cross_hub_top50, aes(x = HubScore, y = reorder(row.names(cross_hub_top50), HubScore))) + 
  geom_bar(stat = 'identity',fill = '#4DAF4A') +labs(x = 'Hub scores', y = '')
ggsave('Cross_kingdom_network_top50_hubscore.pdf', cross_hub_top50_p, width = 6, height = 8)

############计算Hub，最费时间的一步##################
cross_zipi <- ZiPiPlot(igraph = cross_igraph, method = 'cluster_fast_greedy')
ggsave('Cross_kingdom_network_zipi.pdf', cross_zipi[[1]], width = 8, height = 6)
cross_zipi[[1]]

cross_hub_info <- data.frame(cross_node_property,HubScore = cross_hub[row.names(cross_node_property),],cross_zipi[[2]][row.names(cross_node_property),])
write.csv(cross_hub_info, 'Cross_kingdom_network_hub_info.csv')
######################################################

cross_network_cor <- as.matrix(cross_network_cor)

cross_mt2 <- model_maptree2(cor = cross_network_cor, method = 'cluster_fast_greedy')
node_mt2 <- cross_mt2[[1]]

node_mt2 <- merge(node_mt2, cross_node_property, by.x = 'row.names', by.y = 'row.names')
model_group <- cross_mt2[[2]][1:2]
names(model_group)[2] <- 'Group'
names(node_mt2)[1] <- 'KO' 
node_mt2 <- merge(node_mt2, model_group, by.x = 'KO', by.y = 'ID')
names(node_mt2)[5] <- 'Degree'
write.csv(node_mt2,file = "node_mt2.csv",row.names = F)

edge_mt2 <- edgeBuild(cor = cross_network_cor, node = node_mt2)
names(edge_mt2)[20] <- 'Correlation'
write.csv(edge_mt2,"edge_mt2.csv",row.names = F)
