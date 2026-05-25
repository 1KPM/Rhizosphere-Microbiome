### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

library(ggplotify)
library(ggplot2)
library(cowplot)
###########################################
# 定义分析函数
analyze_plant_order <- function(order_name) {
  # 提取当前植物目的ID
  order_id <- subset(metadata_tree, Order == order_name)$TreeID
  
  # 定义颜色方案（在函数外部定义以避免重复）
  if(!exists("distance_colors")) {
    distance_colors <- c("Bray Curtis" = "#609ac6", 
                         "Jaccard" = "#db6786", 
                         "Weighted UniFrac" = "#7bc47f", 
                         "Unweighted UniFrac" = "#d4a356")
    
    distance_fill_colors <- c("Bray Curtis" = "#d6dde2", 
                              "Jaccard" = "#fbf3d5", 
                              "Weighted UniFrac" = "#e1f5e1", 
                              "Unweighted UniFrac" = "#f5e9d1")
  }
  
  # 定义结果存储列表
  all_results <- data.frame()
  plots_list <- list()
  spearman_results_list <- list()
  
  # 定义微生物组类型
  microbe_types <- list(
    "16S" = list(
      bray = bray_16S_mean,
      jaccard = jaccard_16S_mean,
      weighted = weighted_16S_mean,
      unweighted = unweighted_16S_mean,
      suffix = "16S",
      group_name = "Bacteria"
    ),
    "ITS" = list(
      bray = bray_ITS_mean,
      jaccard = jaccard_ITS_mean,
      # 标准化weighted ITS距离矩阵
      weighted = {
        normalize_distance_matrix <- function(dist_matrix) {
          max_val <- max(dist_matrix, na.rm = TRUE)
          if (max_val > 0) {
            return(dist_matrix / max_val)
          } else {
            return(dist_matrix)
          }
        }
        normalize_distance_matrix(weighted_ITS_mean)
      },
      unweighted = unweighted_ITS_mean,
      suffix = "ITS",
      group_name = "Fungi"
    ),
    "18S" = list(
      bray = bray_18S_mean,
      jaccard = jaccard_18S_mean,
      weighted = weighted_18S_mean,
      unweighted = unweighted_18S_mean,
      suffix = "18S",
      group_name = "Protist"
    )
  )
  
  # 对每种微生物组类型进行分析
  for(microbe_name in names(microbe_types)) {
    microbe_data <- microbe_types[[microbe_name]]
    
    # 提取当前植物目在当前微生物组中的样本ID
    microbe_id <- intersect(order_id, row.names(microbe_data$bray))
    
    # 检查是否有足够样本
    if(length(microbe_id) < 3) {
      message(paste("警告:", order_name, microbe_data$group_name, "样本数不足3，跳过分析"))
      next
    }
    
    # 准备距离矩阵
    microbe_groups <- list(
      bray = microbe_data$bray[microbe_id, microbe_id],
      jaccard = microbe_data$jaccard[microbe_id, microbe_id],
      weighted = microbe_data$weighted[microbe_id, microbe_id],
      unweighted = microbe_data$unweighted[microbe_id, microbe_id]
    )
    
    # 准备数据
    data_microbe <- prepare_distance_pairs(tree_dist, microbe_groups)
    
    # 计算Spearman相关系数
    spearman_results <- list()
    for (dist_name in c("Bray Curtis", "Jaccard", "Weighted UniFrac", "Unweighted UniFrac")) {
      if (dist_name %in% names(data_microbe)) {
        x <- data_microbe$tree_dist
        y <- data_microbe[[dist_name]]
        spearman_results[[dist_name]] <- calculate_spearman_correlation(x, y)
      }
    }
    
    # 存储结果
    spearman_results_list[[paste0(order_name, "_", microbe_data$suffix)]] <- spearman_results
    
    # 准备绘图数据
    plot_data <- data_microbe %>%
      pivot_longer(cols = -tree_dist, 
                   names_to = "distance_metric", 
                   values_to = "microbe_dist")
    
    # 创建散点图
    p <- plot_data %>%
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
           y = paste0("Taxonomic similarity (", order_name, "_", microbe_data$suffix, ")"))
    
    # 保存单个图形
    name <- paste0("./plots/Scatter_tree_vs_", order_name, "_", microbe_data$suffix, "_distances_spearman")
    width <- 6
    height <- 5
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 300, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
    
    # 存储图形
    plots_list[[microbe_name]] <- p
    
    # 整理结果
    for (dist_name in names(spearman_results)) {
      result_row <- data.frame(
        Order = order_name,
        Microbe_Group = microbe_data$group_name,
        Distance_Metric = dist_name,
        Spearman_rho = spearman_results[[dist_name]]$spearman_rho,
        P_value = spearman_results[[dist_name]]$spearman_p
      )
      all_results <- rbind(all_results, result_row)
    }
  }
  
  # 如果至少有一种微生物组有数据，创建组合图形
  if (length(plots_list) > 0) {
    # 创建图例
    p_legend <- plots_list[[1]] +
      theme(legend.position = "right") +
      scale_color_manual(values = distance_colors, name = "Distance Metric") +
      scale_fill_manual(values = distance_fill_colors, name = "Distance Metric") +
      theme(legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9, face = "bold"))
    
    legend_grob <- get_legend(p_legend)
    legend_plot <- as.ggplot(legend_grob)
    
    # 创建组合图形
    combined_plot <- Reduce(`+`, c(plots_list, list(legend_plot))) + 
      plot_layout(widths = c(rep(3, length(plots_list)), 3))
    
    # 保存组合图形
    name <- paste0("./plots/", order_name, "_Combined_scatter_tree_vs_microbial_distances")
    width <- 6 * (length(plots_list) + 1)
    height <- 6
    ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 300, units = "cm")
    ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")
  }
  
  # 保存当前目的的结果
  if (nrow(all_results) > 0) {
    write.csv(all_results, 
              paste0("./stats/", order_name, "_ASV_Spearman_correlation_results.csv"), 
              row.names = FALSE, 
              fileEncoding = "GBK")
  }
  
  return(all_results)
}

###########################################
# 循环处理所有植物目
orders_to_analyze <- c("Sapindales", "Arecales", "Asparagales", "Gentianales", "Fabales",
                       "Lamiales", "Malpighiales", "Malvales", "Myrtales", "Rosales")

# 存储所有结果
all_orders_results <- data.frame()

# 创建必要的目录
if (!dir.exists("./plots")) dir.create("./plots")
if (!dir.exists("./stats")) dir.create("./stats")

# 循环分析每个植物目
for (order in orders_to_analyze) {
  message(paste("正在分析:", order))
  order_results <- analyze_plant_order(order)
  all_orders_results <- rbind(all_orders_results, order_results)
}

# 保存所有结果
if (nrow(all_orders_results) > 0) {
  write.csv(all_orders_results, 
            "./stats/All_Orders_ASV_Spearman_correlation_results.csv", 
            row.names = FALSE, 
            fileEncoding = "GBK")
  message("所有植物目分析完成，结果已保存到 ./stats/All_Orders_ASV_Spearman_correlation_results.csv")
} else {
  message("警告：没有生成任何分析结果")
}

# 可选：创建一个汇总所有植物目的图形
if (nrow(all_orders_results) > 0) {
  library(ggplot2)
  
  # 创建汇总图：每个植物目的Spearman相关系数
  summary_plot <- ggplot(all_orders_results, 
                         aes(x = Order, y = Spearman_rho, fill = Microbe_Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ Distance_Metric, ncol = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Plant Order", y = "Spearman Correlation (rho)", 
         title = "Spearman Correlation between Plant Phylogenetic Distance and Microbial Distance",
         fill = "Microbial Group")
  
  ggsave("./plots/Summary_All_Orders_Spearman_Correlations.png", 
         summary_plot, width = 12, height = 8, dpi = 300)
  ggsave("./plots/Summary_All_Orders_Spearman_Correlations.pdf", 
         summary_plot, width = 12, height = 8)
  
  message("汇总图形已保存到 ./plots/Summary_All_Orders_Spearman_Correlations")
}