library(microeco)
#library(rgexf)
library(servr)
library(ggplot2)

#KO_Bacteria <- read.csv('../data/Bacteria_core0.2.csv', row.names = 1,check.names = F)
#KO_Fungi <- read.csv('../data/Fungi_core0.2.csv', row.names = 1,check.names = F)
#KO_Protist <- read.csv('../data/Protist_core0.2.csv', row.names = 1,check.names = F)

#KO_all <- rbind(KO_Bacteria, KO_Fungi, KO_Protist)

metadata <- read.csv("/data/work/data/metadata.csv", row.names = 1, header = T)
KO_all <- read.csv("/data/work/data/All_core0.2.csv", row.names = 1, check.names = F)
#cor_mat <- read.csv("../data/All_KO_WGCNA_spearman_cor_raw.csv", row.names = 1, check.names = F)
#p_mat <- read.csv("../data/All_KO_WGCNA_spearman_p_raw.csv", row.names = 1, check.names = F)
KO_dataset_all <- microtable$new(sample_table = metadata, otu_table = KO_all)
KO_dataset_all$tidy_dataset()

KO_network_All <- trans_network$new(dataset = KO_dataset_all, cor_method = 'spearman',use_WGCNA_pearson_spearman = T, filter_thres = 0, nThreads = 18)
write.csv(KO_network_All$res_cor_p$cor, '/data/work/network/All_KO_WGCNA_spearman_cor_raw.csv')
write.csv(KO_network_All$res_cor_p$p, '/data/work/network/All_KO_WGCNA_spearman_p_raw.csv')

KO_network_All$cal_network(
  COR_p_thres = 0.01,
#  COR_optimization = TRUE,
  COR_cut = 0.8
)

KO_network_All$cal_module(method = 'cluster_fast_greedy')

KO_network_All$cal_network_attr()
write.csv(KO_network_All$res_network_attr, '/data/work/network/KO_network_property_All.csv')

KO_network_All$get_node_table()
write.csv(KO_network_All$res_node_table, '/data/work/network/KO_network_node_property_All.csv', row.names = F)
KO_network_All$get_edge_table()
write.csv(KO_network_All$res_edge_table, '/data/work/network/KO_network_edge_property_All.csv', row.names = F)
KO_network_All$get_adjacency_matrix()
write.csv(KO_network_All$res_adjacency_matrix, '/data/work/network/KO_network_adjacency_matrix_All')

#KO_network_All$res_node_table <- na.omit(KO_network_All$res_node_table)
#KO_network_hub_All <- KO_network_All$plot_taxa_roles(use_type = 1, add_label = T)
#ggsave('/data/work/network/KO_network_hub_All.pdf', KO_network_hub_All, width = 8, height = 6)
