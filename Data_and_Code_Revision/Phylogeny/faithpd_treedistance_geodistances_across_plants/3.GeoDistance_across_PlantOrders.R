### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# 加载包
library(dplyr)
library(ggplot2)
library(rstatix)
library(agricolae)
library(geosphere)

# 读取数据
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]

# 计算地理距离矩阵
tree_coords <- aggregate(cbind(Longitude, Latitude) ~ TreeID, metadata_tree, mean)
row.names(tree_coords) <- tree_coords$TreeID
geo_dist_all <- distm(tree_coords[, c("Longitude", "Latitude")], fun = distGeo) / 1000
row.names(geo_dist_all) <- tree_coords$TreeID
colnames(geo_dist_all) <- tree_coords$TreeID

# 提取每个目的组内地理距离
order_geo_dist_data <- data.frame()
for (group in c(top10orders, "Others")) {
  if (group == "Others") {
    group_ids <- subset(metadata_tree, !(Order %in% top10orders))$TreeID
  } else {
    group_ids <- subset(metadata_tree, Order == group)$TreeID
  }
  group_ids <- intersect(group_ids, row.names(geo_dist_all))
  
  if (length(group_ids) >= 2) {
    dist_values <- geo_dist_all[group_ids, group_ids][lower.tri(diag(length(group_ids)))]
    order_geo_dist_data <- rbind(order_geo_dist_data, 
                                 data.frame(Group = group, Dist = dist_values))
  }
}

# 计算统计
calculate_group_stats <- function(data) {
  summarise_df <- data %>%
    group_by(Group) %>%
    summarise(
      n = n(), mean = mean(Dist), sd = sd(Dist), max = max(Dist),
      se = sd / sqrt(n), ci = qt(0.975, n - 1) * se
    ) %>%
    mutate(allmax = max(max))
  
  kruskal_test <- kruskal.test(Dist ~ Group, data = data)
  dunn_test <- rstatix::dunn_test(data, Dist ~ Group, p.adjust.method = "fdr")
  
  n <- nrow(summarise_df)
  pvalue_matrix <- matrix(1, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      p <- dunn_test$p.adj[dunn_test$group1 == levels(data$Group)[i] & 
                             dunn_test$group2 == levels(data$Group)[j]]
      pvalue_matrix[i, j] <- pvalue_matrix[j, i] <- ifelse(length(p) > 0, p, 1)
    }
  }
  
  letter_df <- agricolae::orderPvalue(summarise_df$Group, summarise_df$mean, 0.05, 
                                      pvalue_matrix, console = FALSE)
  summarise_df$label <- letter_df$groups
  
  list(summary_df = summarise_df, kruskal_test = kruskal_test, dunn_test = dunn_test)
}

# 绘图
order_geo_dist_data$Group <- factor(order_geo_dist_data$Group, 
                                    levels = c(top10orders, "Others"))
stats_result <- calculate_group_stats(order_geo_dist_data)

geo_dist_g <- ggplot(order_geo_dist_data, aes(x = Group, y = Dist, color = Group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
  geom_text(data = stats_result$summary_df,
            aes(y = allmax * 1.1, label = label), size = 8/2.835) +
  labs(y = "Plant geographic distance (km)") +
  scale_color_manual(values = top_tree$Color2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")

# 保存结果
dir.create("./geo_stats", showWarnings = FALSE)
dir.create("./geo_plots", showWarnings = FALSE)
write.csv(order_geo_dist_data, "./geo_stats/plant_geographic_top10orders_dist_data.csv", row.names = FALSE)
ggsave("./geo_plots/Plant_geographic_distance_top10orders.png", geo_dist_g, 
       width = 8, height = 6, dpi = 600, units = "cm")
ggsave("./geo_plots/Plant_geographic_distance_top10orders.pdf", geo_dist_g, 
       width = 8, height = 6, dpi = 600, units = "cm")

# 输出统计摘要
cat("Kruskal-Wallis p-value:", stats_result$kruskal_test$p.value, "\n")