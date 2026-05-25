library(igraph)
library(dplyr)
#library(tidyfst)

##输入数据来自单个目总网络
KO_network_cor_all <- read.csv('../network/All_KO_WGCNA_spearman_cor_raw.csv', row.names = 1)
KO_network_node_all <- read.csv('../network/KO_network_node_property_All.csv')
KO_network_edge_all <- read.csv('../network/KO_network_edge_property_All.csv')
KO_subset <- read.csv('../data/All_core0.2.csv', row.names = 1,check.names = F)


cross_network_edge <- KO_network_edge_all[gsub('K.*', '', KO_network_edge_all$node1) != gsub('K.*', '', KO_network_edge_all$node2),]
write.csv(cross_network_edge, 'cross_network_edge.csv', row.names = F)

cross_network_node <- KO_network_node_all[KO_network_node_all$name %in% c(cross_network_edge$node1, cross_network_edge$node2),]
write.csv(cross_network_node, 'cross_network_node.csv', row.names = F)

cross_network_node <- read.csv("cross_network_node.csv")
cross_network_edge <- read.csv("cross_network_edge.csv")

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

#cross_hub <- data.frame(HubScore = hub_score(cross_igraph)$vector)
#cross_hub_top50 <- head(arrange(cross_hub, desc(HubScore)), 50)

#cross_hub_top50_p <- ggplot(cross_hub_top50, aes(x = HubScore, y = reorder(row.names(cross_hub_top50), HubScore))) + 
#  geom_bar(stat = 'identity',fill = '#4DAF4A') +labs(x = 'Hub scores', y = '')
#ggsave('Cross_kingdom_network_top50_hubscore.pdf', cross_hub_top50_p, width = 6, height = 8)
