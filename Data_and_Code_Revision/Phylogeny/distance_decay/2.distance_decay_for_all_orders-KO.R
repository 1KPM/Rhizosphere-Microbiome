########################################################################################
# (12) 绘制多个植物目的KO功能距离衰减曲线
# 定义函数：分析单个植物目的KO功能距离衰减
analyze_KO_for_order <- function(order_name, metadata_tree, tree_dist, 
                                 bray_bac_KO, jaccard_bac_KO,
                                 bray_fun_KO, jaccard_fun_KO,
                                 bray_pro_KO, jaccard_pro_KO,
                                 legend_plot = NULL) {
  
  # 提取当前植物目的ID
  order_id <- subset(metadata_tree, Order == order_name)$TreeID
  
  # 设置颜色方案
  distance_colors <- c("Bray Curtis" = "#609ac6", 
                       "Jaccard" = "#db6786")
  
  distance_fill_colors <- c("Bray Curtis" = "#d6dde2", 
                            "Jaccard" = "#fbf3d5")
  
  # 存储结果和图形
  all_results <- data.frame()
  plots_list <- list()
  
  # 定义微生物组类型及其对应的数据
  microbe_types <- list(
    "Bacteria" = list(
      bray = bray_bac_KO,
      jaccard = jaccard_bac_KO,
      suffix = "bac_KO",
      group_name = "Bacteria_KO"
    ),
    "Fungi" = list(
      bray = bray_fun_KO,
      jaccard = jaccard_fun_KO,
      suffix = "fun_KO",
      group_name = "Fungi_KO"
    ),
    "Protist" = list(
      bray = bray_pro_KO,
      jaccard = jaccard_pro_KO,
      suffix = "pro_KO",
      group_name = "Protist_KO"
    )
  )
  
  # 对每个微生物组进行分析
  for (microbe_name in names(microbe_types)) {
    microbe_data <- microbe_types[[microbe_name]]
    
    # 提取当前植物目在当前微生物组中的样本ID
    microbe_id <- intersect(order_id, row.names(microbe_data$bray))
    
    # 检查是否有足够样本
    if (length(microbe_id) < 3) {
      message(paste("警告:", order_name, microbe_data$group_name, "样本数不足3，跳过分析"))
      next
    }
    
    # 准备距离矩阵
    microbe_groups <- list(
      bray = microbe_data$bray[microbe_id, microbe_id],
      jaccard = microbe_data$jaccard[microbe_id, microbe_id]
    )
    
    # 准备数据
    data_microbe <- prepare_distance_pairs(tree_dist, microbe_groups)
    
    # 计算Spearman相关系数
    spearman_results <- list()
    for (dist_name in c("Bray Curtis", "Jaccard")) {
      if (dist_name %in% names(data_microbe)) {
        x <- data_microbe$tree_dist
        y <- data_microbe[[dist_name]]
        spearman_results[[dist_name]] <- calculate_spearman_correlation(x, y)
      }
    }
    
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
           y = paste0("Functional similarity (", order_name, "_", 
                      microbe_data$group_name, ")"))
    
    # 保存单个图形
    name <- paste0("./meta_plots/Scatter_tree_vs_", order_name, "_", 
                   microbe_data$suffix, "_distances_spearman")
    width <- 6
    height <- 5
    ggsave(paste0(name, ".png"), p, width = width, height = height, 
           dpi = 300, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, 
           units = "cm")
    
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
  
  # 创建图例（如果未提供）
  if (is.null(legend_plot) && length(plots_list) > 0) {
    p_legend <- plots_list[[1]] +
      theme(legend.position = "right") +
      scale_color_manual(values = distance_colors, name = "Distance Metric") +
      scale_fill_manual(values = distance_fill_colors, name = "Distance Metric") +
      theme(legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9, face = "bold"))
    
    legend_grob <- get_legend(p_legend)
    legend_plot <- as.ggplot(legend_grob)
  }
  
  # 创建组合图形
  if (length(plots_list) > 0) {
    combined_plot <- Reduce(`+`, c(plots_list, list(legend_plot))) + 
      plot_layout(widths = c(rep(3, length(plots_list)), 3))
    
    # 保存组合图形
    name <- paste0("./meta_plots/", order_name, 
                   "_Combined_scatter_tree_vs_functional_distances")
    width <- 6 * (length(plots_list) + 1)
    height <- 6
    ggsave(paste0(name, ".png"), combined_plot, width = width, 
           height = height, dpi = 300, units = "cm")
    ggsave(paste0(name, ".pdf"), combined_plot, width = width, 
           height = height, units = "cm")
  }
  
  return(list(results = all_results, plots = plots_list))
}

########################################################################################
# 创建图例（单独创建一次，供所有植物目使用）
create_KO_legend <- function() {
  # 使用一个临时数据框创建图例
  temp_data <- data.frame(
    distance_metric = factor(c("Bray Curtis", "Jaccard"), 
                             levels = c("Bray Curtis", "Jaccard"))
  )
  
  distance_colors <- c("Bray Curtis" = "#609ac6", 
                       "Jaccard" = "#db6786")
  
  distance_fill_colors <- c("Bray Curtis" = "#d6dde2", 
                            "Jaccard" = "#fbf3d5")
  
  p_legend <- ggplot(temp_data, aes(x = 1, y = 1, color = distance_metric, 
                                    fill = distance_metric)) +
    geom_point(shape = 21, size = 3) +
    scale_color_manual(values = distance_colors, name = "Distance Metric") +
    scale_fill_manual(values = distance_fill_colors, name = "Distance Metric") +
    theme_void() +
    theme(legend.position = "right",
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9, face = "bold"))
  
  legend_grob <- get_legend(p_legend)
  legend_plot <- as.ggplot(legend_grob)
  
  # 保存图例
  ggsave("./meta_plots/Distance_metric_legend_KO.png", legend_plot, 
         width = 4, height = 3, dpi = 300, units = "cm")
  ggsave("./meta_plots/Distance_metric_legend_KO.pdf", legend_plot, 
         width = 4, height = 3, units = "cm")
  
  return(legend_plot)
}

########################################################################################
# 主分析流程
# 定义要分析的植物目
orders_to_analyze <- c("Arecales", "Asparagales", "Fabales", "Gentianales", 
                       "Lamiales", "Malpighiales", "Malvales", "Myrtales", 
                       "Rosales", "Sapindales")

# 确保输出目录存在
if (!dir.exists("./meta_plots")) dir.create("./meta_plots", recursive = TRUE)
if (!dir.exists("./meta_stats")) dir.create("./meta_stats", recursive = TRUE)

# 创建图例
KO_legend <- create_KO_legend()

# 存储所有结果
all_orders_KO_results <- data.frame()

# 循环分析每个植物目
for (order in orders_to_analyze) {
  message(paste("正在分析:", order, "的KO功能距离衰减..."))
  
  result <- analyze_KO_for_order(
    order_name = order,
    metadata_tree = metadata_tree,
    tree_dist = tree_dist,
    bray_bac_KO = bray_bac_KO,
    jaccard_bac_KO = jaccard_bac_KO,
    bray_fun_KO = bray_fun_KO,
    jaccard_fun_KO = jaccard_fun_KO,
    bray_pro_KO = bray_pro_KO,
    jaccard_pro_KO = jaccard_pro_KO,
    legend_plot = KO_legend
  )
  
  # 保存当前植物目的结果
  if (nrow(result$results) > 0) {
    all_orders_KO_results <- rbind(all_orders_KO_results, result$results)
    
    # 保存当前植物目的单独结果
    write.csv(result$results, 
              paste0("./meta_stats/", order, "_KO_Spearman_correlation_results.csv"), 
              row.names = FALSE, 
              fileEncoding = "GBK")
    
    message(paste(order, "分析完成，结果已保存"))
  } else {
    message(paste("警告:", order, "没有生成任何分析结果"))
  }
}

# 保存所有植物目的汇总结果
if (nrow(all_orders_KO_results) > 0) {
  write.csv(all_orders_KO_results, 
            "./meta_stats/All_Orders_KO_Spearman_correlation_results.csv", 
            row.names = FALSE, 
            fileEncoding = "GBK")
  message("所有植物目KO功能距离衰减分析完成，汇总结果已保存到 ./meta_stats/All_Orders_KO_Spearman_correlation_results.csv")
} else {
  message("警告：没有生成任何KO功能距离衰减分析结果")
}

########################################################################################
# 创建汇总图形
if (nrow(all_orders_KO_results) > 0) {
  # 创建汇总图：每个植物目的Spearman相关系数
  summary_plot_KO <- all_orders_KO_results %>%
    mutate(Microbe_Group = factor(Microbe_Group, 
                                  levels = c("Bacteria_KO", "Fungi_KO", "Protist_KO"))) %>%
    ggplot(aes(x = Order, y = Spearman_rho, fill = Distance_Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ Microbe_Group, ncol = 1) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    ) +
    scale_fill_manual(values = c("Bray Curtis" = "#609ac6", "Jaccard" = "#db6786")) +
    labs(
      x = "Plant Order", 
      y = "Spearman Correlation Coefficient (ρ)", 
      title = "Spearman Correlation between Plant Phylogenetic Distance and KO Functional Distance",
      fill = "Distance Metric"
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3)
  
  # 保存汇总图形
  ggsave("./meta_plots/Summary_All_Orders_KO_Spearman_Correlations.png", 
         summary_plot_KO, width = 16, height = 12, dpi = 300, units = "cm")
  ggsave("./meta_plots/Summary_All_Orders_KO_Spearman_Correlations.pdf", 
         summary_plot_KO, width = 16, height = 12, units = "cm")
  
  message("汇总图形已保存到 ./meta_plots/Summary_All_Orders_KO_Spearman_Correlations")
  
  # 创建显著性标记的版本
  summary_plot_KO_sig <- all_orders_KO_results %>%
    mutate(
      Microbe_Group = factor(Microbe_Group, 
                             levels = c("Bacteria_KO", "Fungi_KO", "Protist_KO")),
      Significance = ifelse(P_value < 0.001, "***", 
                            ifelse(P_value < 0.01, "**", 
                                   ifelse(P_value < 0.05, "*", ""))),
      y_pos = ifelse(Spearman_rho >= 0, Spearman_rho + 0.02, Spearman_rho - 0.05)
    ) %>%
    ggplot(aes(x = Order, y = Spearman_rho, fill = Distance_Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(y = y_pos, label = Significance), 
              position = position_dodge(width = 0.8), 
              size = 3, vjust = 0) +
    facet_wrap(~ Microbe_Group, ncol = 1) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    ) +
    scale_fill_manual(values = c("Bray Curtis" = "#609ac6", "Jaccard" = "#db6786")) +
    labs(
      x = "Plant Order", 
      y = "Spearman Correlation Coefficient (ρ)", 
      title = "Spearman Correlation between Plant Phylogenetic Distance and KO Functional Distance",
      subtitle = "*P<0.05, **P<0.01, ***P<0.001",
      fill = "Distance Metric"
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3)
  
  ggsave("./meta_plots/Summary_All_Orders_KO_Spearman_Correlations_with_sig.png", 
         summary_plot_KO_sig, width = 16, height = 12, dpi = 300, units = "cm")
  ggsave("./meta_plots/Summary_All_Orders_KO_Spearman_Correlations_with_sig.pdf", 
         summary_plot_KO_sig, width = 16, height = 12, units = "cm")
  
  message("带显著性标记的汇总图形已保存到 ./meta_plots/Summary_All_Orders_KO_Spearman_Correlations_with_sig")
}