
### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

dir.create("./phylo_stats", showWarnings = FALSE)
dir.create("./phylo_plots", showWarnings = FALSE)

# Import package
library(dplyr)
library(tidyr)
library(rstatix)
library(agricolae)
library(stats)
library(ggplot2)
library(ggtree)
library(ape)
library(geosphere)

### Import data ----------------------------------------------------------------
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

### Keep sample with tree_ITS_reads
metadata_rs <- metadata_rs[metadata_rs$TreeID %in% tree_file$tip.label, ]
match(tree_file$tip.label, unique(metadata_rs$TreeID)) ###说明有的树木样本没有扩增子测序数据

tree_dist <- as.data.frame(cophenetic.phylo(tree_file)) 
###cophenetic.phylo computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.

######################################################################################################
# 读取top10植物目信息
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]

ggtheme <- 
  theme(
    text = element_text(color = "black", size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
  )

######################################################################################
###（1）不同植物目ITS序列的相似性（距离）比较
# 循环处理top10目
# 初始化一个空的数据框来存储所有组的结果

order_dist_data <- data.frame(Group = character(), Dist = numeric())

# 处理每个分类组
for (group in c(top10orders,"Others")) {
  cat("正在处理:", group, "\n")
  
  ###########################################
  # 1. 筛选该目的植物ID
  ###########################################
  if (group == "Others") {
    # 对于"Others"组，选择不属于top10orders的样本
    group_ids <- subset(metadata_tree, !(Order %in% top10orders))$TreeID
    group_ids <- intersect(group_ids, names(tree_dist))
    
  } else if (group %in% metadata_tree$Order) {
    # 对于具体的目，选择该目下的样本
    group_ids <- subset(metadata_tree, Order == group)$TreeID
    group_ids <- intersect(group_ids, names(tree_dist))
  }
  else {
    next
  }
  
  ###########################################
  # 2. 提取系统发育距离
  ###########################################
  # 直接提取距离矩阵
  group_tree_dist <- tree_dist[group_ids, group_ids]
  
  # 提取下三角部分（不包括对角线）
  extract_lower_tri <- function(mat) {
    mat[lower.tri(mat)]
  }
  
  # 提取距离值
  dist_values <- extract_lower_tri(as.matrix(group_tree_dist))
  
  # 创建当前组的数据框
  temp_df <- data.frame(
    Group = group,
    Dist = dist_values
  )
  
  # 合并到总数据框中
  order_dist_data <- rbind(order_dist_data, temp_df)
}

write.csv(order_dist_data, file = "./phylo_stats/plant_ITS_top10orders_dist_data.csv")

####可视化
order_dist_data$Group <- factor(order_dist_data$Group, levels = c(top10orders, "Others"))

### 计算组间距离统计和多重比较结果
calculate_group_stats <- function(data, value_col, group_col, alpha = 0.05, p_adjust_method = "fdr") {
  
  # 检查列名是否存在
  if (!value_col %in% names(data)) {
    stop(paste("列名", value_col, "不存在于数据框中"))
  }
  if (!group_col %in% names(data)) {
    stop(paste("列名", group_col, "不存在于数据框中"))
  }
  
  # 重命名列以便于处理
  names(data)[names(data) == value_col] <- "Dist"
  names(data)[names(data) == group_col] <- "Group"
  
  # 确保Group是因子
  if (!is.factor(data$Group)) {
    data$Group <- factor(data$Group)
  }
  
  ###########################################
  # 1. 计算每组的基本统计量
  ###########################################
  summarise_df <- data %>%
    group_by(Group) %>%
    summarise(
      n = length(Dist),
      mean = mean(Dist, na.rm = TRUE),
      min = min(Dist, na.rm = TRUE),
      max = max(Dist, na.rm = TRUE),
      sd = sd(Dist, na.rm = TRUE),
      se = sd / sqrt(n),
      ci = qt(1 - alpha/2, df = n - 1) * se
    ) %>%
    mutate(allmax = max(max))
  
  ###########################################
  # 2. 非参数检验
  ###########################################
  # Kruskal-Wallis检验
  kruskal_test <- kruskal.test(Dist ~ Group, data = data)
  sig_value <- kruskal_test$p.value
  
  # Dunn事后检验
  dunn_test <- rstatix::dunn_test(data, Dist ~ Group, p.adjust.method = p_adjust_method)
  
  # 添加显著性标记
  pvalue_df <- data.frame(
    KruskalWallis = sig_value, 
    dunn_test[c("group1", "group2", "p.adj")]
  )
  pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                            ifelse(pvalue_df$p.adj < 0.01, "**", 
                                   ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
  
  ###########################################
  # 3. 生成p值矩阵
  ###########################################
  n <- nrow(summarise_df)
  pvalue_matrix <- matrix(1, ncol = n, nrow = n)
  k <- 0
  
  for(i in 1:(n - 1)) { 
    for(j in (i + 1):n) { 
      k <- k + 1
      pvalue_matrix[i, j] <- pvalue_df$p.adj[k]
      pvalue_matrix[j, i] <- pvalue_df$p.adj[k]
    }
  }
  
  ###########################################
  # 4. 生成显著性字母标记
  ###########################################
  # 使用agricolae包
  letter_df <- agricolae::orderPvalue(
    summarise_df$Group, 
    summarise_df$mean, 
    alpha, 
    pvalue_matrix, 
    console = FALSE
  )
  
  # 确保字母标记的顺序与summarise_df一致
  letter_df <- letter_df[levels(data$Group), ]
  summarise_df$label <- letter_df$groups
  
  ###########################################
  # 5. 返回结果列表
  ###########################################
  result_list <- list(
    summary_df = summarise_df,
    kruskal_test = kruskal_test,
    dunn_test = dunn_test,
    letter_df = letter_df,
    pvalue_matrix = pvalue_matrix
  )
  
  return(result_list)
}

top10orders_dist_summarise <- calculate_group_stats(data = order_dist_data, value_col = "Dist", group_col = "Group")
###由于R包的更新问题，重新打开跑脚本

summary_df <- top10orders_dist_summarise$summary_df
write.csv(summary_df, file = "./phylo_stats/plant_ITS_top10orders_dist_summarise_df.csv")

tree_dist_g <- ggplot(data = order_dist_data, mapping = aes(x = Group, y = Dist)) +
  geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = summary_df,
    mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
    position = position_dodge(0.9),
    size = 8 / 2.835
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = paste0("Plant phylogenetic distance")
  ) + 
  scale_color_manual(values = top_tree$Color2) +
  theme_bw() + ggtheme

name <- "./phylo_plots/Plant_phylogenetic_distance_top10orders"
width <- 8
height <- 6
ggsave(paste0(name, ".png"), tree_dist_g, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), tree_dist_g, width = width, height = height, units = "cm")