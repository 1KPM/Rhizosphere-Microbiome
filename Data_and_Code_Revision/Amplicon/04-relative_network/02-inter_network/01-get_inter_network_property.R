### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "01-get_inter_network_property"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(igraph)
library(ggClusterNet)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- "inter"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
edge <- read.csv("../01-all_network/01-get_all_network_data/all_edge_property_raw.csv", header = T)
node <- read.csv("../01-all_network/01-get_all_network_data/all_node_property_raw.csv", header = T)
r_table <- read.csv("../01-all_network/01-get_all_network_data/WGCNA_spearman_r_table.csv", header = T, row.names = 1)
core_table <- read.csv("../../01-sort_data/05-core_table/All_core_feature_table_relative.csv", row.names = 1)
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
network_color <- read.csv("../../01-sort_data/08-network_color/top_5percent_taxa.csv", header = T)
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
# 只关注跨界，去除界内相关数据
tmp_edge <- edge[gsub('ASV.*', '', edge$node1) != gsub('ASV.*', '', edge$node2),]
write.csv(tmp_edge, paste0(dir_name, '/', type, '_edge_property_raw.csv'), row.names = F)

tmp_node <- node[node$name %in% c(tmp_edge$node1, tmp_edge$node2),]
write.csv(tmp_node, paste0(dir_name, '/', type, '_node_property_raw.csv'), row.names = F)

# 整理相关性表格
tmp_r_table <- r_table[tmp_node$name, tmp_node$name]
tmp_r_table <- tmp_r_table * 0
for (i in 1:nrow(tmp_edge)) {
    tmp_r_table[tmp_edge$node1[i], tmp_edge$node2[i]] <- tmp_edge$weight[i]
    tmp_r_table[tmp_edge$node2[i], tmp_edge$node1[i]] <- tmp_edge$weight[i]
}
for (i in names(tmp_r_table)) {tmp_r_table[i, i] <- 1}
write.csv(tmp_r_table, paste0(dir_name, '/', type, '_r_table.csv'))

otu_table <- core_table[names(tmp_r_table),]
tax_table <- core_taxonomy[names(tmp_r_table),]

write.csv(otu_table, paste0(dir_name, '/', type, '_otu_table.csv'))
write.csv(tax_table, paste0(dir_name, '/', type, '_tax_table.csv'))

# 定义一个简易清理函数：去除向量中的 NA 和 空字符串 ("")
clean_taxa <- function(x) {
    x[!is.na(x) & x != ""]
}

# 批量清理并直接释放为全局变量！
# lapply(network_color, clean_taxa) 会返回一个干净的列表。
# list2env() 会把这个列表中的元素，直接按照列名（如 top_16s_phylum）变成你环境里的独立变量。
list2env(lapply(network_color, clean_taxa), envir = .GlobalEnv)

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
igraph <- graph_from_data_frame(tmp_edge, directed = F, vertices = tmp_node)

# 计算网络特征
net_property <- net_properties(igraph)
write.csv(net_property, paste0(dir_name, '/', type, '_netwrok_property.csv'))

node_property <- node_properties(igraph)
zipi <- ZiPiPlot(igraph = igraph, method = 'cluster_fast_greedy')
ggsave(paste0(dir_name, '/', type, '_network_zipi.pdf'), zipi[[1]], width = 8, height = 6)

hub_info <- data.frame(node_property,zipi[[2]][row.names(node_property),])
write.csv(hub_info, paste0(dir_name, '/', type, '_network_hub_info.csv'))


# 按照maptree内聚算法改进离散点排布
r_matrix <- as.matrix(tmp_r_table)
maptree <- model_maptree2(cor = r_matrix, method = 'cluster_fast_greedy')
maptree_node <- maptree[[1]]
maptree_edge <- edgeBuild(cor = r_matrix, node = maptree_node)
names(maptree_edge)[names(maptree_edge) == 'cor'] <- 'Correlation'

# node节点注释
maptree_node <- nodeadd(plotcord = maptree_node, otu_table = otu_table, tax_table = tax_table)
node_property <- data.frame(node_property)
names(node_property) <- c('Degree', 'Closeness', 'Betweenness', 'CenDegree')
maptree_node <- merge(maptree_node, data.frame(node_property), by.x = 'ASVID', by.y = 'row.names')

maptree_node <- maptree_node %>%
    mutate(
        # 1. 整理 Phylum 因子：如果不在 top 列表中，则归为 Other
        Phylum = case_when(
            Clade == 'Bacteria' & !(Phylum %in% top_16s_phylum) ~ 'Other bacteria',
            Clade == 'Fungi'    & !(Phylum %in% top_its_phylum) ~ 'Other fungi',
            Clade == 'Protist'  & !(Phylum %in% top_pro_phylum) ~ 'Other protists',
            TRUE ~ Phylum  # 其余情况（即在 top 列表中的）保持原名
        ),
        
        # 2. 整理 Order 因子：如果不在 top 列表中，则归为 Other
        Order = case_when(
            Clade == 'Bacteria' & !(Order %in% top_16s_order) ~ 'Other bacteria',
            Clade == 'Fungi'    & !(Order %in% top_its_order) ~ 'Other fungi',
            Clade == 'Protist'  & !(Order %in% top_pro_order) ~ 'Other protists',
            TRUE ~ Order   # 其余情况保持原名
        ),
        
        # 3. 转换为因子并严格指定 levels 顺序
        Phylum = factor(Phylum, levels = c(
            top_16s_phylum, 'Other bacteria', 
            top_its_phylum, 'Other fungi',
            top_pro_phylum, 'Other protists'
        )),
        
        Order = factor(Order, levels = c(
            top_16s_order, 'Other bacteria', 
            top_its_order, 'Other fungi',
            top_pro_order, 'Other protists'
        ))
    )

# 按Phylum着色
color_phylum <- colorRampPalette(brewer.pal(12, 'Paired'))(13)

p <- ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Correlation), 
                 alpha = 0.01, linewidth = 0.1, data = maptree_edge) + 
    geom_point(aes(X1, X2, fill = Phylum, size = Degree), pch = 21, data = maptree_node) +
    scale_fill_manual(values = color_phylum) +
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(breaks = NULL) +
    theme_classic() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'right')

width <- 22
height <- 18
name <- paste0(dir_name, '/', type, '_network_maptree_Phylum')

ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------
