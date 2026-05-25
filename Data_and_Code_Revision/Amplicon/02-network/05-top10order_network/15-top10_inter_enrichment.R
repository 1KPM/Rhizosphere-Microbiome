### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "15-top10_inter_enrichment"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(clusterProfiler)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

order_list <-  c("Fabales","Arecales","Malpighiales","Sapindales","Lamiales",
                 "Rosales","Asparagales","Myrtales","Malvales","Gentianales")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)

order2asv <- read.csv('../../01-sort_data/09-enrichment_database/order2asv.csv')
order2kingdom <- read.csv('../../01-sort_data/09-enrichment_database/order2kingdom.csv')

# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
net <- "inter"

# 1. Hub ASV 富集分析
hub_list <- list()
for (ord in names(top10order_list)) {
    hub_info <- read.csv(paste0("03-get_top10order_network_property/", ord, "_", net, "_network_hub_info.csv"), row.names = 1)
    hub_asv <- row.names(hub_info)[hub_info$roles != "Peripherals"]
    hub_list[[ord]] <- data.frame(t(hub_asv))
}

hub_df <- data.frame(t(rbind.fill(hub_list)))
names(hub_df) <- names(top10order_list)
name <- paste0(dir_name, '/top10order_hub_asv_raw')
write.csv(hub_df, paste0(name, '_list.csv'), row.names = F)

res <- compareCluster(hub_df, fun = 'enricher', TERM2GENE = order2asv)
res@compareClusterResult <- merge(res@compareClusterResult, order2kingdom, by.x = 'ID', by.y = 'Order', all.x = T)
write.csv(res@compareClusterResult, paste0(name, '_enrichment_results.csv'), row.names = F)

remove_list <- c('Unidentified Bacteria', 'Unidentified Fungi', 'Unidentified Protist')
res@compareClusterResult <- res@compareClusterResult[!(res@compareClusterResult$ID %in% remove_list),]
data_df <- res@compareClusterResult

data_df[data_df$ID == "bacteriap25", "Description"] <- "Myxococcota bacteriap25"
data_df[data_df$ID == "Subgroup_2", "Description"] <- "Acidobacteriota subgroup_2"
data_df[data_df$ID == "MB-A2-108", "Description"] <- "Actinobacteriota MB-A2-108"

data_df <- data_df %>%
    mutate(
        GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )

data_df$Cluster <- factor(data_df$Cluster, levels = names(top10order_list))


p <- ggplot(data_df, aes(x = Description, y = Cluster)) +
    geom_point(aes(size = GeneRatio, fill = p.adjust), shape = 21) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) + 
    scale_fill_gradient(low = "#e17674", high = "#3b7eb8", name = "p.adjust") + # 颜色映射p值
    scale_y_discrete(drop = FALSE, limits = rev) +
    scale_size_continuous(range = c(1, 3), name = "ASV Ratio") +
    facet_grid(cols = vars(Kingdom), scales = "free", space = "free_x") +
    theme_bw() +
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 7, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          # legend.box.margin = margin(0, 0, 0, -40),
          legend.key.size = unit(0.25, 'cm'),
          panel.spacing = unit(0.1, "cm"),
          legend.box.spacing = unit(0.1,"cm"),
          legend.position = "top") +
    guides(
        fill = guide_colorbar(
            label.theme = element_text(angle = -90, size = 6, color = 'black', vjust = 0.5)
        ),
        size = guide_legend(nrow = 2) 
    )

width = 8.5
height = 8.6
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


# 2. Stable Hub ASV 富集分析 (level 1：与hub相连的ASV)
stable_hub_asv <- read.csv(paste0("05-get_stable_hub_asv/inter_stable_hub_asv.csv"), header = T)
stable_hub_asv <- stable_hub_asv$x

asv_list <- list()
for (asv in stable_hub_asv) {
    network_nodes <- list()
    
    for (ord in order_list) {
        file_path <- paste0("../05-top10order_network/03-get_top10order_network_property/", ord, "_inter_edge_property_raw.csv")
        
        if (file.exists(file_path)) {
            raw_edge <- read.csv(file_path)
            
            tmp_edge <- raw_edge[raw_edge$node1 %in% asv | raw_edge$node2 %in% asv, ]
            
            if (nrow(tmp_edge) > 0) {
                # 提取所有连接的节点（排除ASV本身）
                connected_nodes <- unique(c(tmp_edge$node1[tmp_edge$node1 != asv], 
                                            tmp_edge$node2[tmp_edge$node2 != asv]))
                
                if (length(connected_nodes) > 0) {
                    network_nodes[[ord]] <- connected_nodes
                }
            }
        }
    }
    
    # 取所有网络中连接节点的交集
    if (length(network_nodes) > 0) {
        # 使用Reduce函数取多个向量的交集
        common_nodes <- Reduce(intersect, network_nodes)
        asv_list[[asv]] <- common_nodes
    } else {
        asv_list[[asv]] <- character(0)  # 如果没有连接，返回空向量
    }
}

# 将列表转换为数据框（每个ASV一列）
# 首先找到最长的向量长度
max_len <- max(sapply(asv_list, length))

# 创建数据框，用NA填充缺失值
asv_df <- data.frame(lapply(asv_list, function(x) {
    if (length(x) == 0) {
        return(rep(NA, max_len))
    } else if (length(x) < max_len) {
        return(c(x, rep(NA, max_len - length(x))))
    } else {
        return(x)
    }
}))

# 设置列名
colnames(asv_df) <- stable_hub_asv

# 保存结果
name <- paste0(dir_name, "/", 'top10_inter_stable_hub_asv_level1')
write.csv(asv_df, paste0(name, '_list.csv'), row.names = F)

res <- compareCluster(asv_df, fun = 'enricher', TERM2GENE = order2asv)
res@compareClusterResult <- merge(res@compareClusterResult, order2kingdom, by.x = 'ID', by.y = 'Order', all.x = T)
write.csv(res@compareClusterResult, paste0(name, '_enrichment_results.csv'), row.names = F)

remove_list <- c('Unidentified Bacteria', 'Unidentified Fungi', 'Unidentified Protist')
res@compareClusterResult <- res@compareClusterResult[!(res@compareClusterResult$ID %in% remove_list),]
data_df <- res@compareClusterResult

data_df[data_df$ID == "bacteriap25", "Description"] <- "Myxococcota bacteriap25"
data_df[data_df$ID == "Subgroup_2", "Description"] <- "Acidobacteriota subgroup_2"
data_df[data_df$ID == "MB-A2-108", "Description"] <- "Actinobacteriota MB-A2-108"
data_df[data_df$ID == "Gammaproteobacteria_Incertae_Sedis", "Description"] <- "Gammaproteobacteria"

data_df <- data_df %>%
    mutate(
        GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
data_df <- merge(data_df, core_taxonomy[names(core_taxonomy) != "Kingdom"], by.x = "Cluster", by.y = "ASVID")
data_df$Cluster <- factor(data_df$Cluster, levels = stable_hub_asv)
data_df$ASVLabel <- gsub("_.*", "", data_df$ASVLabel)
ASVLabel_level <- gsub("_.*", "", core_taxonomy[core_taxonomy$ASVID %in% stable_hub_asv, "ASVLabel"])
data_df$ASVLabel <-factor(data_df$ASVLabel, ASVLabel_level)

p <- ggplot(data_df, aes(x = Description, y = ASVLabel)) +
    geom_point(aes(size = GeneRatio, fill = p.adjust), shape = 21) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) + 
    scale_fill_gradient(low = "#e17674", high = "#3b7eb8", name = "p.adjust") + # 颜色映射p值
    scale_y_discrete(drop = FALSE, limits = rev) +
    scale_size_continuous(range = c(1, 3), name = "ASV Ratio") +
    facet_grid(cols = vars(Kingdom), scales = "free", space = "free_x") +
    theme_bw() +
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 7, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          # legend.box.margin = margin(0, 0, 0, -40),
          legend.key.size = unit(0.25, 'cm'),
          panel.spacing = unit(0.1, "cm"),
          legend.box.spacing = unit(0.1,"cm"),
          legend.position = "top") +
    guides(
        fill = guide_colorbar(
            label.theme = element_text(angle = -90, size = 6, color = 'black', vjust = 0.5)
        ),
        size = guide_legend(nrow = 2) 
    )

width = 8.5
height = 8
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

# ------------------------------------------------------------------------------
