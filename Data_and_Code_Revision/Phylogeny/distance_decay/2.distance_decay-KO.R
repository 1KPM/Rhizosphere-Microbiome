### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
library(ggplot2)
library(reshape2)
library(agricolae)
library(picante)
library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(gridExtra)
library(ggplotify)
library(geosphere)

### Import data ----------------------------------------------------------------
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

tree_dist <- as.data.frame(cophenetic.phylo(tree_file))

####################################################################################### 
###(1) bacteria KO
### Import and filter data
###
bray_bac_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Bacteria.quantitative.bray.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(bray_bac_KO) <- gsub("RS.", "", names(bray_bac_KO))
row.names(bray_bac_KO) <- gsub("RS-", "", row.names(bray_bac_KO))
bray_bac_KO <- bray_bac_KO[row.names(bray_bac_KO) %in% metadata_tree$TreeID, names(bray_bac_KO) %in% metadata_tree$TreeID]

###
jaccard_bac_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Bacteria.quantitative.jaccard.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(jaccard_bac_KO) <- gsub("RS.", "", names(jaccard_bac_KO))
row.names(jaccard_bac_KO) <- gsub("RS-", "", row.names(jaccard_bac_KO))
jaccard_bac_KO <- jaccard_bac_KO[row.names(jaccard_bac_KO) %in% metadata_tree$TreeID, names(jaccard_bac_KO) %in% metadata_tree$TreeID]

#######################################################################################
###(2) fungi KO
### Import and filter data
###
bray_fun_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Fungi.quantitative.bray.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(bray_fun_KO) <- gsub("RS.", "", names(bray_fun_KO))
row.names(bray_fun_KO) <- gsub("RS-", "", row.names(bray_fun_KO))
bray_fun_KO <- bray_fun_KO[row.names(bray_fun_KO) %in% metadata_tree$TreeID, names(bray_fun_KO) %in% metadata_tree$TreeID]

###
jaccard_fun_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Fungi.quantitative.jaccard.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(jaccard_fun_KO) <- gsub("RS.", "", names(jaccard_fun_KO))
row.names(jaccard_fun_KO) <- gsub("RS-", "", row.names(jaccard_fun_KO))
jaccard_fun_KO <- jaccard_fun_KO[row.names(jaccard_fun_KO) %in% metadata_tree$TreeID, names(jaccard_fun_KO) %in% metadata_tree$TreeID]

#######################################################################################
###(3) protist KO
### Import and filter data
###
bray_pro_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Protist.quantitative.bray.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(bray_pro_KO) <- gsub("RS.", "", names(bray_pro_KO))
row.names(bray_pro_KO) <- gsub("RS-", "", row.names(bray_pro_KO))
bray_pro_KO <- bray_pro_KO[row.names(bray_pro_KO) %in% metadata_tree$TreeID, names(bray_pro_KO) %in% metadata_tree$TreeID]

###
jaccard_pro_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Protist.quantitative.jaccard.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(jaccard_pro_KO) <- gsub("RS.", "", names(jaccard_pro_KO))
row.names(jaccard_pro_KO) <- gsub("RS-", "", row.names(jaccard_pro_KO))
jaccard_pro_KO <- jaccard_pro_KO[row.names(jaccard_pro_KO) %in% metadata_tree$TreeID, names(jaccard_pro_KO) %in% metadata_tree$TreeID]

########################################################################################
# (10) 绘制distance decay曲线
# 准备数据：提取距离矩阵的下三角部分
prepare_distance_pairs <- function(tree_dist_mat, microbe_mat_list) {
  # 确保所有矩阵有相同的顺序
  all_names <- Reduce(intersect, lapply(microbe_mat_list, function(x) row.names(x)))
  all_names <- intersect(all_names, row.names(tree_dist_mat))
  
  # 排序所有矩阵
  tree_dist_aligned <- tree_dist_mat[all_names, all_names]
  microbe_mats_aligned <- lapply(microbe_mat_list, function(mat) mat[all_names, all_names])
  
  # 提取下三角部分（不包括对角线）
  extract_lower_tri <- function(mat) {
    mat[lower.tri(mat)]
  }
  
  # 创建数据框
  tree_vec <- extract_lower_tri(as.matrix(tree_dist_aligned))
  
  # 创建结果数据框
  result_df <- data.frame(tree_dist = tree_vec)
  
  # 添加每个微生物距离
  for (i in seq_along(microbe_mats_aligned)) {
    mat_name <- names(microbe_mats_aligned)[i]
    microbe_vec <- extract_lower_tri(as.matrix(microbe_mats_aligned[[i]]))
    
    # 提取距离度量名称
    dist_metric <- gsub(".*_", "", mat_name)
    if (grepl("^bray", mat_name, ignore.case = TRUE)) {
      dist_metric <- "Bray Curtis"
    } else if (grepl("^jaccard", mat_name, ignore.case = TRUE)) {
      dist_metric <- "Jaccard"
    }
    
    result_df[[dist_metric]] <- microbe_vec
  }
  
  return(result_df)
}

# 计算Spearman相关系数和P值
calculate_spearman_correlation <- function(x, y) {
  # 计算Spearman相关系数
  cor_result <- cor.test(x, y, method = "spearman")
  
  return(list(
    spearman_rho = cor_result$estimate,
    spearman_p = cor_result$p.value
  ))
}

# 准备细菌KO数据
microbe_groups_bac_KO <- list(
  bray_bac_KO = bray_bac_KO,
  jaccard_bac_KO = jaccard_bac_KO
)

data_bac_KO <- prepare_distance_pairs(tree_dist, microbe_groups_bac_KO)

# 计算KO的Spearman相关系数
spearman_results_bac_KO <- list()
for (dist_name in c("Bray Curtis", "Jaccard")) {
  if (dist_name %in% names(data_bac_KO)) {
    x <- data_bac_KO$tree_dist
    y <- data_bac_KO[[dist_name]]
    spearman_results_bac_KO[[dist_name]] <- calculate_spearman_correlation(x, y)
  }
}

###########################################################
# 准备KO绘图数据
plot_data_bac_KO <- data_bac_KO %>%
  pivot_longer(cols = -tree_dist, 
               names_to = "distance_metric", 
               values_to = "microbe_dist")

# 设置颜色方案
distance_colors <- c("Bray Curtis" = "#609ac6", 
                     "Jaccard" = "#db6786")

distance_fill_colors <- c("Bray Curtis" = "#d6dde2", 
                          "Jaccard" = "#fbf3d5")

# 创建散点图
p_bac_KO <- plot_data_bac_KO %>%
  ggplot(aes(x = tree_dist, y = 1 - microbe_dist)) +
  geom_point(aes(color = distance_metric, fill = distance_metric),
             shape = 21, size = 1, alpha = 0.6) +
  scale_color_manual(values = distance_colors) +
  scale_fill_manual(values = distance_fill_colors) +
  geom_smooth(aes(color = distance_metric),
              method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 8)) +
  labs(x = "Phylogenetic distance",
       y = "Functional similarity (Bacteria)")

###
name <- "./meta_plots/Scatter_tree_vs_bac_KO_distances_spearman"
width <- 6
height <- 5
ggsave(paste0(name, ".png"), p_bac_KO, width = width, height = height, dpi = 300, units = "cm")
ggsave(paste0(name, ".pdf"), p_bac_KO, width = width, height = height, units = "cm")


# 准备fun_KO数据
microbe_groups_fun_KO <- list(
  bray_fun_KO = bray_fun_KO,
  jaccard_fun_KO = jaccard_fun_KO
)

data_fun_KO <- prepare_distance_pairs(tree_dist, microbe_groups_fun_KO)

# 计算fun_KO的Spearman相关系数
spearman_results_fun_KO <- list()
for (dist_name in c("Bray Curtis", "Jaccard")) {
  if (dist_name %in% names(data_fun_KO)) {
    x <- data_fun_KO$tree_dist
    y <- data_fun_KO[[dist_name]]
    spearman_results_fun_KO[[dist_name]] <- calculate_spearman_correlation(x, y)
  }
}

# 准备fun_KO绘图数据
plot_data_fun_KO <- data_fun_KO %>%
  pivot_longer(cols = -tree_dist, 
               names_to = "distance_metric", 
               values_to = "microbe_dist")

# 创建散点图
p_fun_KO <- plot_data_fun_KO %>%
  ggplot(aes(x = tree_dist, y = 1 - microbe_dist)) +
  geom_point(aes(color = distance_metric, fill = distance_metric),
             shape = 21, size = 1, alpha = 0.6) +
  scale_color_manual(values = distance_colors) +
  scale_fill_manual(values = distance_fill_colors) +
  geom_smooth(aes(color = distance_metric),
              method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 8)) +
  labs(x = "Phylogenetic distance",
       y = "Functional similarity (Fungi)")

###
name <- "./meta_plots/Scatter_tree_vs_fun_KO_distances_spearman"
width <- 6
height <- 5
ggsave(paste0(name, ".png"), p_fun_KO, width = width, height = height, dpi = 300, units = "cm")
ggsave(paste0(name, ".pdf"), p_fun_KO, width = width, height = height, units = "cm")


# 准备18S数据
microbe_groups_pro_KO <- list(
  bray_pro_KO = bray_pro_KO,
  jaccard_pro_KO = jaccard_pro_KO
)

data_pro_KO <- prepare_distance_pairs(tree_dist, microbe_groups_pro_KO)

# 计算18S的Spearman相关系数
spearman_results_pro_KO <- list()
for (dist_name in c("Bray Curtis", "Jaccard")) {
  if (dist_name %in% names(data_pro_KO)) {
    x <- data_pro_KO$tree_dist
    y <- data_pro_KO[[dist_name]]
    spearman_results_pro_KO[[dist_name]] <- calculate_spearman_correlation(x, y)
  }
}

# 准备18S绘图数据
plot_data_pro_KO <- data_pro_KO %>%
  pivot_longer(cols = -tree_dist, 
               names_to = "distance_metric", 
               values_to = "microbe_dist")

# 创建散点图
p_pro_KO <- plot_data_pro_KO %>%
  ggplot(aes(x = tree_dist, y = 1 - microbe_dist)) +
  geom_point(aes(color = distance_metric, fill = distance_metric),
             shape = 21, size = 1, alpha = 0.6) +
  scale_color_manual(values = distance_colors) +
  scale_fill_manual(values = distance_fill_colors) +
  geom_smooth(aes(color = distance_metric),
              method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 8)) +
  labs(x = "Phylogenetic distance",
       y = "Functional similarity (Protists)")

###
name <- "./meta_plots/Scatter_tree_vs_pro_KO_distances_spearman"
width <- 6
height <- 5
ggsave(paste0(name, ".png"), p_pro_KO, width = width, height = height, dpi = 300, units = "cm")
ggsave(paste0(name, ".pdf"), p_pro_KO, width = width, height = height, units = "cm")

######################################################################################
# 创建一个带图例的图形用于提取图例
p_legend <- plot_data_bac_KO %>%
  ggplot(aes(x = tree_dist, y = 1 - microbe_dist)) +
  geom_point(aes(color = distance_metric, fill = distance_metric),
             shape = 21, size = 1, alpha = 0.6) +
  scale_color_manual(values = distance_colors, name = "Distance Metric") +
  scale_fill_manual(values = distance_fill_colors, name = "Distance Metric") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.spacing.y = unit(0.1, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 21)))

# 提取图例
legend_grob <- get_legend(p_legend)

# 将图例转换为ggplot对象
legend_plot <- as.ggplot(legend_grob)

# 保存图例单独的图形
ggsave("./meta_plots/KO_Distance_metric_legend.png", legend_plot, 
       width = 4, height = 3, dpi = 300, units = "cm")
ggsave("./meta_plots/KO_Distance_metric_legend.pdf", legend_plot, 
       width = 4, height = 3, dpi = 300, units = "cm")

###
# 使用patchwork组合图形
combined_plot <- p_bac_KO + p_fun_KO + p_pro_KO + legend_plot +
  plot_layout(widths = c(3,3,3,3))

# 保存组合图形
name <- "./meta_plots/Combined_scatter_tree_vs_functional_distances"
width <- 20
height <- 6
ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 300, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")


# 将相关性检验结果保存为表格
# 整理bac_KO结果
results_bac_KO <- data.frame(
  Microbe_Group = "Bacteria (KO)",
  Distance_Metric = names(spearman_results_bac_KO),
  Spearman_rho = sapply(spearman_results_bac_KO, function(x) x$spearman_rho),
  P_value = sapply(spearman_results_bac_KO, function(x) x$spearman_p)
)

# 整理fun_KO结果
results_fun_KO <- data.frame(
  Microbe_Group = "Fungi (KO)",
  Distance_Metric = names(spearman_results_fun_KO),
  Spearman_rho = sapply(spearman_results_fun_KO, function(x) x$spearman_rho),
  P_value = sapply(spearman_results_fun_KO, function(x) x$spearman_p)
)

# 整理pro_KO结果
results_pro_KO <- data.frame(
  Microbe_Group = "Protist (KO)",
  Distance_Metric = names(spearman_results_pro_KO),
  Spearman_rho = sapply(spearman_results_pro_KO, function(x) x$spearman_rho),
  P_value = sapply(spearman_results_pro_KO, function(x) x$spearman_p)
)

# 合并所有结果
all_results <- rbind(results_bac_KO, results_fun_KO, results_pro_KO)

# 保存结果到CSV文件
write.csv(all_results, "./meta_stats/KO_Spearman_correlation_results.csv", row.names = FALSE, fileEncoding = "GBK")

