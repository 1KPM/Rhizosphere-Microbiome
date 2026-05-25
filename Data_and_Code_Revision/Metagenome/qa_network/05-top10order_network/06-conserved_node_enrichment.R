### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Create directory
dir_name <- "06-conserved_node_enrichment"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# 1. Stable Hub ko 富集分析 (level 1：与hub相连的ko)
stable_hub_ko <- read.csv(paste0("05-get_stable_hub_ko/inter_stable_hub_ko.csv"), header = T)
stable_hub_ko <- stable_hub_ko$x

### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
kegg_name <- read.csv("00-data/K_gene_name.csv",row.names = 1)

pathway2ko <- read.csv('00-data/pathway2ko.csv')
pathway2superpathway <- read.csv('00-data/pathway2superpathway.csv')



# 加载必要的包
library(plyr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(data.table)

# 定义所有网络列表
All_list <- c("Fabales","Arecales","Malpighiales","Sapindales","Lamiales",
              "Rosales","Asparagales","Myrtales","Malvales","Gentianales")

# 假设 stable_hub_ko 已经定义
# stable_hub_ko <- c("KO1", "KO2", ...)  # 示例

# 创建一个列表来存储每个KO在所有网络中的连接节点
ko_list <- list()

# 遍历每个稳定的hub KO
for (ko in stable_hub_ko) {
  # 存储每个网络中与该KO连接的节点
  network_nodes <- list()
  
  # 遍历每个网络
  for (i in All_list) {
    # 读取网络边文件
    file_path <- paste0("../05-top10order_network/03-get_top10order_network_property/", 
                        i, "_inter_edge_property_raw.csv")
    
    if (file.exists(file_path)) {
      raw_edge <- read.csv(file_path)
      
      # 找出与该KO连接的节点
      tmp_edge <- raw_edge[raw_edge$node1 %in% ko | raw_edge$node2 %in% ko, ]
      
      if (nrow(tmp_edge) > 0) {
        # 提取所有连接的节点（排除KO本身）
        connected_nodes <- unique(c(tmp_edge$node1[tmp_edge$node1 != ko], 
                                    tmp_edge$node2[tmp_edge$node2 != ko]))
        
        if (length(connected_nodes) > 0) {
          network_nodes[[i]] <- connected_nodes
        }
      }
    }
  }
  
  # 取所有网络中连接节点的交集
  if (length(network_nodes) > 0) {
    # 使用Reduce函数取多个向量的交集
    common_nodes <- Reduce(intersect, network_nodes)
    ko_list[[ko]] <- common_nodes
  } else {
    ko_list[[ko]] <- character(0)  # 如果没有连接，返回空向量
  }
}

# 将列表转换为数据框（每个KO一列）
# 首先找到最长的向量长度
max_len <- max(sapply(ko_list, length))

# 创建数据框，用NA填充缺失值
ko_df <- data.frame(lapply(ko_list, function(x) {
  if (length(x) == 0) {
    return(rep(NA, max_len))
  } else if (length(x) < max_len) {
    return(c(x, rep(NA, max_len - length(x))))
  } else {
    return(x)
  }
}))

# 设置列名
colnames(ko_df) <- stable_hub_ko

# 保存结果
name <- paste0(dir_name, "/", 'top10_inter_stable_hub_ko_level1')

write.csv(ko_df, paste0(name, '_list.csv'), row.names = F)


##仅保留与其相连的细菌KO
ko_df <- ko_df %>%
  mutate(across(everything(), ~{
    cell <- as.character(.)
    ifelse(grepl("^b", cell), cell, "")
  })) %>%
  filter(!if_all(everything(), ~. == ""))

ko_df <- ko_df %>% mutate(across(everything(), ~substr(as.character(.), 2, nchar(.))))

ko_res <- compareCluster(ko_df, fun = 'enricher', TERM2GENE = pathway2ko)

ko_res@compareClusterResult <- merge(ko_res@compareClusterResult, pathway2superpathway, by.x = 'ID', by.y = 'Pathway', all.x = T)
write.csv(ko_res@compareClusterResult, paste0(name, '_enrichment_results.csv'), row.names = F)


############################################################################
###从这里开始
#################################
name <- paste0(dir_name, "/", 'top10_inter_stable_hub_ko_level1')

data_ko_df <- read.csv(paste0(name, '_enrichment_results.csv'))

data_ko_df <- data_ko_df %>%
  mutate(
    GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )

data_ko_df <- data_ko_df %>% filter(p.adjust < 0.01)

# ##本图片调整
# data_ko_df <- data_ko_df %>%
#   group_by(ID) %>%
#   filter(n() >= 3) %>%
#   ungroup()
# 


stable_hub_ko <- unique(data_ko_df$Cluster)
data_ko_df$Cluster <- factor(data_ko_df$Cluster, levels = stable_hub_ko)
data_ko_df$Cluster <- as.character(data_ko_df$Cluster)
data_ko_df$Cluster_temp <- substr(data_ko_df$Cluster, 2, nchar(data_ko_df$Cluster))
data_ko_df$gene_match <- kegg_name$gene[match(data_ko_df$Cluster_temp, kegg_name$KO)]
data_ko_df$KOLabel <- paste(data_ko_df$Cluster, data_ko_df$gene_match, sep = "|")
data_ko_df$Cluster_temp <- NULL
data_ko_df$gene_match <- NULL
data_ko_df$log.p.adjust <- -log10(data_ko_df$p.adjust)

#############行列聚类优化绘图顺序###############################################
# 将数据转换为矩阵格式
data_ko_df <- data_ko_df %>%
  mutate(KOLabel = gsub(" ", "", KOLabel) %>%
           gsub("fK02913|RP-L33,MRPL33,rpmG", "fK02913|RP-L33,rpmG", ., fixed = TRUE) %>%
           gsub("fK02950|RP-S12,MRPS12,rpsL", "fK02950|RP-S12,rpsL", ., fixed = TRUE) %>%
           gsub("fK00525|E1.17.4.1A,nrdA,nrdE", "fK00525|nrdA,nrdE", ., fixed = TRUE) %>%
           gsub("fK01676|E4.2.1.2A,fumA,fumB", "fK01676|fumA,fumB", ., fixed = TRUE))


data_ko_df <- setDT(data_ko_df)
cluster_mat <- dcast(data_ko_df, KOLabel ~ Description, value.var = "Count", fill = 0)
cluster_mat <- as.data.frame(cluster_mat)
rownames(cluster_mat) <- cluster_mat$KOLabel
cluster_mat <- cluster_mat[, -1]
cluster_mat <- (cluster_mat > 0) * 1 # 基于有无

# 行聚类：使用欧氏距离和ward.D2方法
row_dist <- dist(cluster_mat, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")
y_order <- rownames(cluster_mat)[row_hclust$order]

plot(row_hclust)

# 列聚类：使用欧氏距离和ward.D2方法
col_dist <- dist(t(cluster_mat), method = "euclidean")
col_hclust <- hclust(col_dist, method = "ward.D2")
x_order <- colnames(cluster_mat)[col_hclust$order]

# 重新排序数据框中的因子水平
data_ko_df <- data_ko_df %>%
  filter(KOLabel %in% y_order, Description %in% x_order) %>%
  mutate(
    KOLabel = factor(KOLabel, levels = rev(y_order)),      # rev() 可根据需要调整
    Description = factor(Description, levels = x_order)  # rev() 可根据需要调整
  )

################################################################################

p <- ggplot(data_ko_df, aes(x = Description, y = KOLabel)) +
  geom_point(
    aes(size = Count * 1.2),
    position = position_nudge(x = 0.12, y = -0.12),
    color = "gray20",
    alpha = 0.4,
    shape = 19,
    show.legend = FALSE
  ) +
  geom_point(
    aes(size = Count, fill = log.p.adjust),
    color = "black",
    shape = 21,
    stroke = 0.8,
    alpha = 1 
  )+
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) + 
  scale_fill_gradientn(
    name = expression(-log[10](p.adjust)),
    colors = c("#FFFFE0", "#E41A1C"),
    values = scales::rescale(c(0, 5, 10, 20, 50)),
    limits = c(0, 50)
  ) +  
  scale_size_continuous(
    range = c(1, 4),
    breaks = c(50, 100, 150,200),
    name = "KO Count"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 6, color = 'black', hjust = 0.5),
    plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 6, color = 'black'), 
    panel.grid.major = element_line(color = "grey90",linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 0.6),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
    legend.title = element_text(size = 6, color = 'black'), 
    legend.text = element_text(size = 6,  color = 'black'), 
    legend.box.margin = margin(0, 0, 0, -20),
    legend.key.size = unit(0.25, 'cm'),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.position = "none"
  )

width = 18
height = 21
name <- paste0(dir_name, '/stable_hub_ko_level1_ALL')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
