### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# 加载必要的包
library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

##Matrix Data Prepare 
kingdom <- c("bac","fun","pro")
for (kin in kingdom){
  tmp_res_all <- read.csv(paste0("./stats/",kin,"_ko_qa_log_tree_metastat_res_diff_KO_core0.2.csv"),row.names = 1)
  diff_sig <- tmp_res_all %>% filter(pvalue < 0.05 & qvalue < 0.2)
  order_list1 <- c("Others","Myrtales","Rosales","Arecales","Sapindales","Lamiales","Malpighiales",
                   "Malvales","Fabales","Gentianales","Asparagales")
  comparison_matrix <- matrix(0, 
                              nrow = length(order_list1),
                              ncol = length(order_list1),
                              dimnames = list(order_list1, order_list1))
  
  # 创建存储每个目富集和排斥KO列表的容器
  enriched_ko_lists <- setNames(replicate(length(order_list1), list(), simplify = FALSE), order_list1)
  depleted_ko_lists <- setNames(replicate(length(order_list1), list(), simplify = FALSE), order_list1)
  
  for (i in 1:(length(order_list1)-1)) {
    for (j in (i+1):length(order_list1)) {
      current_comparison <- paste(order_list1[i], "-", order_list1[j])
      subset_data <- diff_sig %>% filter(Comparison == current_comparison)
      
      # 富集KO
      up_ko <- subset_data %>% filter(Group == order_list1[i]) %>% pull(Taxa)
      down_ko <- subset_data %>% filter(Group == order_list1[j]) %>% pull(Taxa)
      
      comparison_matrix[i, j] <- length(up_ko)
      comparison_matrix[j, i] <- length(down_ko)
      
      # 添加到非冗余列表
      enriched_ko_lists[[order_list1[i]]] <- unique(c(enriched_ko_lists[[order_list1[i]]], up_ko))
      enriched_ko_lists[[order_list1[j]]] <- unique(c(enriched_ko_lists[[order_list1[j]]], down_ko))
      
      depleted_ko_lists[[order_list1[j]]] <- unique(c(depleted_ko_lists[[order_list1[j]]], up_ko))
      depleted_ko_lists[[order_list1[i]]] <- unique(c(depleted_ko_lists[[order_list1[i]]], down_ko))
    }
  }
  diag(comparison_matrix) <- 0
  order_list <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales","Others")
  comparison_matrix <- comparison_matrix[order_list, order_list]
  
  # 计算非冗余总数
  enriched_totals <- sapply(order_list, function(x) length(enriched_ko_lists[[x]]))
  depleted_totals <- sapply(order_list, function(x) length(depleted_ko_lists[[x]]))
  
  # 保存数据
  write.csv(comparison_matrix, file = paste0("./stats/",kin,"_qa_log_metastat_res_diff_KO_matrix(core0.2p0.05q0.2).csv"))
  write.csv(data.frame(Enriched = enriched_totals, Depleted = depleted_totals, row.names = order_list), 
            file = paste0("./stats/",kin,"_qa_log_metastat_res_diff_KO_totals(core0.2p0.05q0.2).csv"))
}


## Plot Matrix with totals
kingdom <- c("bac", "fun", "pro")
subtitle_map <- c(
  "bac" = "Bacteria",
  "fun" = "Fungi", 
  "pro" = "Protists"
)
plot_list <- list()

# 设置更小的图形尺寸
for (kin in kingdom) {
  # 读取数据
  data <- read.csv(paste0("./stats/",kin,"_qa_log_metastat_res_diff_KO_matrix(core0.2p0.05q0.2).csv"), row.names = 1)
  totals <- read.csv(paste0("./stats/",kin,"_qa_log_metastat_res_diff_KO_totals(core0.2p0.05q0.2).csv"), row.names = 1)
  
  # 获取 order 名称
  order_names <- rownames(data)
  n_orders <- length(order_names)
  
  # 1. 上侧条形图
  top_data <- data.frame(
    Order = factor(order_names, levels = order_names),
    Value = totals[order_names, "Depleted"]
  )
  
  p_top <- ggplot(top_data, aes(x = Order, y = Value)) +
    geom_col(fill = "#5BA2D8", width = 0.8) +  # 减小条形宽度
    geom_text(aes(label = Value), vjust = -0.5, size = 2, color = "black") +  # 减小字体
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 6, color = "black"),  # 减小字体
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6, color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),  # 减小边距
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = NULL, y = "Depleted") +
    scale_x_discrete(limits = order_names, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  # 减小扩展
  
  # 2. 右侧条形图
  right_data <- data.frame(
    Order = factor(order_names, levels = rev(order_names)),
    Value = totals[order_names, "Enriched"]
  )
  
  p_right <- ggplot(right_data, aes(x = Order, y = Value)) +
    geom_col(fill = "#EC499A", width = 0.8) +  # 减小条形宽度
    geom_text(aes(label = Value), vjust = 0.5,  hjust = 0, size = 2, color = "black") +  # 减小字体
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, color = "black"),  # 减小字体
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 6, color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),  # 减小边距
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = NULL, y = "Enriched") +
    scale_x_discrete(limits = rev(order_names), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # 减小扩展
    coord_flip()
  
  # 3. 主热图 - 关键优化
  data_matrix <- as.matrix(data)
  for (i in 1:n_orders) {
    for (j in 1:n_orders) {
      if (i > j) {
        data_matrix[i, j] <- -data_matrix[i, j]
      }
    }
  }
  data_long <- melt(data_matrix)
  colnames(data_long) <- c("Row", "Column", "value")
  
  data_long$transformed_value <- apply(data_long, 1, function(x) {
    row_name <- as.character(x["Row"])
    col_name <- as.character(x["Column"])
    
    if (row_name == col_name) {
      return(NA)
    } else if (match(row_name, order_names) < match(col_name, order_names)) {
      return(log1p(as.numeric(x["value"])))
    } else {
      return(-log1p(-as.numeric(x["value"])))
    }
  })
  
  p_main <- ggplot(data_long, aes(x = Column, y = Row, fill = transformed_value)) +
    geom_tile(color = "black", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#5BA2D8", 
      mid = "white", 
      high = "#EC499A", 
      midpoint = 0,
      na.value = "white",
      name = "log(number)",
      breaks = c(min(data_long$transformed_value, na.rm = TRUE), max(data_long$transformed_value, na.rm = TRUE)),
      labels = c("Depleted", "Enriched"),
      limits = c(min(data_long$transformed_value, na.rm = TRUE), max(data_long$transformed_value, na.rm = TRUE))
    ) +
    coord_fixed(ratio = 1, expand = FALSE) +  # 关键修改：添加 expand = FALSE
    theme_void() +  # 使用更简洁的主题
    theme(
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),  # 边距全部设为0
      legend.position = "none",
      panel.border = element_blank(),  # 移除面板边框
      panel.spacing = unit(0, "pt"),  # 移除面板间距
      plot.background = element_blank(),  # 移除背景
      panel.background = element_blank(),  # 移除面板背景
      panel.grid = element_blank(),
      axis.text.x = element_text(
        size = 6, 
        angle = 45, 
        hjust = 1, 
        vjust = 1, 
        color = "black",
        margin = margin(0, 0, 0, 0, unit = "cm")  # 坐标轴文本边距设为0
      ),
      axis.text.y = element_text(
        size = 6, 
        hjust = 1, 
        vjust = 0.5, 
        color = "black",
        margin = margin(0, 0, 0, 0, unit = "cm")  # 坐标轴文本边距设为0
      ),
      axis.title.y = element_text(
        size = 6, 
        angle = 90,
        hjust = 0.5, 
        vjust = 0.5, 
        margin = margin(0, 0, 0, 0, unit = "cm")  # 坐标轴标题边距设为0
      ),
      axis.ticks.length = unit(0, "cm")  # 移除刻度线
    ) +
    labs(title = "", x = "", y = subtitle_map[[kin]]) +
    geom_text(
      aes(label = round(abs(value), 1)),
      color = "black",
      size = 2
    ) +
    scale_x_discrete(limits = order_names, expand = expansion(mult = 0, add = 0)) +  # 明确设为0扩展
    scale_y_discrete(limits = rev(order_names), expand = expansion(mult = 0, add = 0)) +  # 明确设为0扩展
    geom_tile(
      data = subset(data_long, Row == Column),
      fill = "gray",
      color = "black",
      linewidth = 0.3,
      show.legend = FALSE
    )
  
  # 4. 空白图
  p_empty <- ggplot() + 
    theme_void() + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  # 使用cowplot对齐
  aligned_main_top <- align_plots(p_main, p_top, align = "v", axis = "lr")
  aligned_main_right <- align_plots(p_main, p_right, align = "h", axis = "bt")
  
  # 组合图形 - 使用紧凑布局
  row1 <- plot_grid(
    aligned_main_top[[2]],  # 上侧条形图
    p_empty,
    ncol = 2,
    rel_widths = c(3.5, 1)  # 调整比例
  )
  
  row2 <- plot_grid(
    aligned_main_top[[1]],  # 主热图
    aligned_main_right[[2]],  # 右侧条形图
    ncol = 2,
    rel_widths = c(4, 1)  # 调整比例
  )
  
  # 组合两行
  combined <- plot_grid(
    row1, row2,
    ncol = 1,
    rel_heights = c(1, 4)  # 调整高度比例
  )
  
  plot_list[[kin]] <- combined
}

# 最终组合
final_plot <- plot_grid(
  plotlist = plot_list,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1, 1)
)

# 显示图形
print(final_plot)

# 保存图形 - 显著减小宽度
ggsave("./plots/diff_KO_matrix_combined.png", final_plot, 
       width = 8, height = 24, dpi = 300, bg = "white", units = "cm", limitsize = FALSE)
ggsave("./plots/diff_KO_matrix_combined.pdf", final_plot, 
       width = 8, height = 24, bg = "white", units = "cm", limitsize = FALSE)

