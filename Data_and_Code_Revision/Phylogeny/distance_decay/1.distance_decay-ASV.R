### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
library(ggplot2)
library(reshape2)
library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(ggpubr)

### Import data ----------------------------------------------------------------
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

### Keep sample with tree_ITS_reads
metadata_rs <- metadata_rs[metadata_rs$TreeID %in% tree_file$tip.label, ]

tree_dist <- as.data.frame(cophenetic.phylo(tree_file))

####################################################################################### 
### 距离矩阵处理函数
calculate_mean_distance <- function(dist_matrix, metadata) {
  dist_matrix$ID <- row.names(dist_matrix)
  dist_matrix_double <- melt(dist_matrix)
  dist_matrix_double <- merge(dist_matrix_double, unique(metadata[,1:2]), 
                              by.x = "ID", by.y = "FileID")
  dist_matrix_double <- merge(dist_matrix_double, unique(metadata[,1:2]), 
                              by.x = "variable", by.y = "FileID")
  dist_matrix_mean <- aggregate(value ~ TreeID.x + TreeID.y, 
                                data = dist_matrix_double, FUN = mean)
  dist_matrix_mean <- dcast(dist_matrix_mean, TreeID.x ~ TreeID.y)
  row.names(dist_matrix_mean) <- dist_matrix_mean$TreeID.x
  dist_matrix_mean <- dist_matrix_mean[,-1]
  return(dist_matrix_mean)
}

### (1) bacteria 16S距离矩阵
bray_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__bray_curtis.csv", fileEncoding = "GBK", row.names = 1)
bray_16S <- bray_16S[row.names(bray_16S) %in% metadata_rs$FileID, names(bray_16S) %in% metadata_rs$FileID]
bray_16S_mean <- calculate_mean_distance(bray_16S, metadata_rs)

jaccard_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__jaccard.csv", fileEncoding = "GBK", row.names = 1)
jaccard_16S <- jaccard_16S[row.names(jaccard_16S) %in% metadata_rs$FileID, names(jaccard_16S) %in% metadata_rs$FileID]
jaccard_16S_mean <- calculate_mean_distance(jaccard_16S, metadata_rs)

weighted_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__weighted.csv", fileEncoding = "GBK", row.names = 1)
weighted_16S <- weighted_16S[row.names(weighted_16S) %in% metadata_rs$FileID, names(weighted_16S) %in% metadata_rs$FileID]
weighted_16S_mean <- calculate_mean_distance(weighted_16S, metadata_rs)

unweighted_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__unweighted.csv", fileEncoding = "GBK", row.names = 1)
unweighted_16S <- unweighted_16S[row.names(unweighted_16S) %in% metadata_rs$FileID, names(unweighted_16S) %in% metadata_rs$FileID]
unweighted_16S_mean <- calculate_mean_distance(unweighted_16S, metadata_rs)

### (2) fungi ITS距离矩阵
bray_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__bray_curtis.csv", fileEncoding = "GBK", row.names = 1)
bray_ITS <- bray_ITS[row.names(bray_ITS) %in% metadata_rs$FileID, names(bray_ITS) %in% metadata_rs$FileID]
bray_ITS_mean <- calculate_mean_distance(bray_ITS, metadata_rs)

jaccard_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__jaccard.csv", fileEncoding = "GBK", row.names = 1)
jaccard_ITS <- jaccard_ITS[row.names(jaccard_ITS) %in% metadata_rs$FileID, names(jaccard_ITS) %in% metadata_rs$FileID]
jaccard_ITS_mean <- calculate_mean_distance(jaccard_ITS, metadata_rs)

weighted_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__weighted.csv", fileEncoding = "GBK", row.names = 1)
weighted_ITS <- weighted_ITS[row.names(weighted_ITS) %in% metadata_rs$FileID, names(weighted_ITS) %in% metadata_rs$FileID]
# 标准化ITS的weighted距离矩阵到0-1范围
normalize_distance_matrix <- function(dist_matrix) {
  max_val <- max(dist_matrix, na.rm = TRUE)
  if (max_val > 0) {
    return(dist_matrix / max_val)
  } else {
    return(dist_matrix)
  }
}
weighted_ITS_mean <- calculate_mean_distance(weighted_ITS, metadata_rs)
weighted_ITS_mean_norm <- normalize_distance_matrix(weighted_ITS_mean)

unweighted_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__unweighted.csv", fileEncoding = "GBK", row.names = 1)
unweighted_ITS <- unweighted_ITS[row.names(unweighted_ITS) %in% metadata_rs$FileID, names(unweighted_ITS) %in% metadata_rs$FileID]
unweighted_ITS_mean <- calculate_mean_distance(unweighted_ITS, metadata_rs)

### (3) protist 18S距离矩阵
bray_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__bray_curtis.csv", fileEncoding = "GBK", row.names = 1)
bray_18S <- bray_18S[row.names(bray_18S) %in% metadata_rs$FileID, names(bray_18S) %in% metadata_rs$FileID]
bray_18S_mean <- calculate_mean_distance(bray_18S, metadata_rs)

jaccard_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__jaccard.csv", fileEncoding = "GBK", row.names = 1)
jaccard_18S <- jaccard_18S[row.names(jaccard_18S) %in% metadata_rs$FileID, names(jaccard_18S) %in% metadata_rs$FileID]
jaccard_18S_mean <- calculate_mean_distance(jaccard_18S, metadata_rs)

weighted_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__weighted.csv", fileEncoding = "GBK", row.names = 1)
weighted_18S <- weighted_18S[row.names(weighted_18S) %in% metadata_rs$FileID, names(weighted_18S) %in% metadata_rs$FileID]
weighted_18S_mean <- calculate_mean_distance(weighted_18S, metadata_rs)

unweighted_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__unweighted.csv", fileEncoding = "GBK", row.names = 1)
unweighted_18S <- unweighted_18S[row.names(unweighted_18S) %in% metadata_rs$FileID, names(unweighted_18S) %in% metadata_rs$FileID]
unweighted_18S_mean <- calculate_mean_distance(unweighted_18S, metadata_rs)

#######################################################################################
### 距离衰减分析和绘图
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
    } else if (grepl("^weighted", mat_name, ignore.case = TRUE)) {
      dist_metric <- "Weighted UniFrac"
    } else if (grepl("^unweighted", mat_name, ignore.case = TRUE)) {
      dist_metric <- "Unweighted UniFrac"
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

# 准备16S数据
microbe_groups_16S <- list(
  bray_16S = bray_16S_mean,
  jaccard_16S = jaccard_16S_mean,
  weighted_16S = weighted_16S_mean,
  unweighted_16S = unweighted_16S_mean
)

data_16S <- prepare_distance_pairs(tree_dist, microbe_groups_16S)

# 计算16S的Spearman相关系数
spearman_results_16S <- list()
for (dist_name in c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac")) {
  if (dist_name %in% names(data_16S)) {
    x <- data_16S$tree_dist
    y <- data_16S[[dist_name]]
    spearman_results_16S[[dist_name]] <- calculate_spearman_correlation(x, y)
  }
}

# 设置颜色方案
distance_colors <- c("Bray Curtis" = "#609ac6", 
                     "Jaccard" = "#db6786", 
                     "Weighted UniFrac" = "#7bc47f", 
                     "Unweighted UniFrac" = "#d4a356")

distance_fill_colors <- c("Bray Curtis" = "#d6dde2", 
                          "Jaccard" = "#fbf3d5", 
                          "Weighted UniFrac" = "#e1f5e1", 
                          "Unweighted UniFrac" = "#f5e9d1")



# 创建16S的图
plot_data_16S <- data_16S %>%
  pivot_longer(cols = -tree_dist, 
               names_to = "distance_metric", 
               values_to = "microbe_dist")

p_16S <- plot_data_16S %>%
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
       y = "Taxonomic similarity (16S)")

# 保存16S图形
ggsave("./plots/Scatter_tree_vs_16S_distances_spearman.png", p_16S, 
       width = 6, height = 5, dpi = 300, units = "cm")
ggsave("./plots/Scatter_tree_vs_16S_distances_spearman.pdf", p_16S, 
       width = 6, height = 5, units = "cm")

# 准备ITS数据
microbe_groups_ITS <- list(
  bray_ITS = bray_ITS_mean,
  jaccard_ITS = jaccard_ITS_mean,
  weighted_ITS = weighted_ITS_mean_norm,
  unweighted_ITS = unweighted_ITS_mean
)

data_ITS <- prepare_distance_pairs(tree_dist, microbe_groups_ITS)

# 计算ITS的Spearman相关系数
spearman_results_ITS <- list()
for (dist_name in c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac")) {
  if (dist_name %in% names(data_ITS)) {
    x <- data_ITS$tree_dist
    y <- data_ITS[[dist_name]]
    spearman_results_ITS[[dist_name]] <- calculate_spearman_correlation(x, y)
  }
}

# 创建ITS的图
plot_data_ITS <- data_ITS %>%
  pivot_longer(cols = -tree_dist, 
               names_to = "distance_metric", 
               values_to = "microbe_dist")

p_ITS <- plot_data_ITS %>%
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
       y = "Taxonomic similarity (ITS)")

# 保存ITS图形
ggsave("./plots/Scatter_tree_vs_ITS_distances_spearman.png", p_ITS, 
       width = 6, height = 5, dpi = 300, units = "cm")
ggsave("./plots/Scatter_tree_vs_ITS_distances_spearman.pdf", p_ITS, 
       width = 6, height = 5, units = "cm")

# 准备18S数据
microbe_groups_18S <- list(
  bray_18S = bray_18S_mean,
  jaccard_18S = jaccard_18S_mean,
  weighted_18S = weighted_18S_mean,
  unweighted_18S = unweighted_18S_mean
)

data_18S <- prepare_distance_pairs(tree_dist, microbe_groups_18S)

# 计算18S的Spearman相关系数
spearman_results_18S <- list()
for (dist_name in c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac")) {
  if (dist_name %in% names(data_18S)) {
    x <- data_18S$tree_dist
    y <- data_18S[[dist_name]]
    spearman_results_18S[[dist_name]] <- calculate_spearman_correlation(x, y)
  }
}

# 创建18S的图
plot_data_18S <- data_18S %>%
  pivot_longer(cols = -tree_dist, 
               names_to = "distance_metric", 
               values_to = "microbe_dist")

p_18S <- plot_data_18S %>%
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
       y = "Taxonomic similarity (18S)")

# 保存18S图形
ggsave("./plots/Scatter_tree_vs_18S_distances_spearman.png", p_18S, 
       width = 6, height = 5, dpi = 300, units = "cm")
ggsave("./plots/Scatter_tree_vs_18S_distances_spearman.pdf", p_18S, 
       width = 6, height = 5, units = "cm")

# 创建图例
p_legend <- plot_data_16S %>%
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

legend_grob <- get_legend(p_legend)
legend_plot <- as_ggplot(legend_grob)

# 保存图例
ggsave("./plots/Distance_metric_legend.png", legend_plot, 
       width = 4, height = 3, dpi = 300, units = "cm")
ggsave("./plots/Distance_metric_legend.pdf", legend_plot, 
       width = 4, height = 3, units = "cm")

# 组合图形
combined_plot <- p_16S + p_ITS + p_18S + legend_plot +
  plot_layout(widths = c(3, 3, 3, 3))

# 保存组合图形
ggsave("./plots/Combined_scatter_tree_vs_microbial_distances.png", combined_plot, 
       width = 20, height = 6, dpi = 300, units = "cm")
ggsave("./plots/Combined_scatter_tree_vs_microbial_distances.pdf", combined_plot, 
       width = 20, height = 6, units = "cm")

# 整理并保存Spearman相关性结果
results_16S <- data.frame(
  Microbe_Group = "Bacteria (16S)",
  Distance_Metric = names(spearman_results_16S),
  Spearman_rho = sapply(spearman_results_16S, function(x) x$spearman_rho),
  P_value = sapply(spearman_results_16S, function(x) x$spearman_p)
)

results_ITS <- data.frame(
  Microbe_Group = "Fungi (ITS)",
  Distance_Metric = names(spearman_results_ITS),
  Spearman_rho = sapply(spearman_results_ITS, function(x) x$spearman_rho),
  P_value = sapply(spearman_results_ITS, function(x) x$spearman_p)
)

results_18S <- data.frame(
  Microbe_Group = "Protist (18S)",
  Distance_Metric = names(spearman_results_18S),
  Spearman_rho = sapply(spearman_results_18S, function(x) x$spearman_rho),
  P_value = sapply(spearman_results_18S, function(x) x$spearman_p)
)

all_results <- rbind(results_16S, results_ITS, results_18S)
write.csv(all_results, "./stats/ASV_Spearman_correlation_results.csv", 
          row.names = FALSE, fileEncoding = "GBK")

cat("分析完成！结果已保存。\n")