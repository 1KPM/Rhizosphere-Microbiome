### Settings -------------------------------------------------------------------
# 设置工作路径
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# 设置随机种子
set.seed(12306)

# 创建结果目录
dir_name <- "results"
if (!file.exists(dir_name)) {
  dir.create(dir_name, recursive = TRUE)
}

# 导入包
library(tidyverse)
library(rstatix)
library(agricolae)
library(RColorBrewer)
library(ggpmisc)
library(ggplot2)
library(cowplot)


### 数据导入函数 ----------------------------------------------------------------
import_data <- function() {
  # 导入元数据
  metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
  metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")
  tree_color_order <- read.csv("../metadata/Tree_top_order_color.csv")
  
  return(list(
    metadata_rs = metadata_rs,
    metadata_tree = metadata_tree,
    tree_color_order = tree_color_order
  ))
}

data_list <- import_data()
metadata_rs <- data_list$metadata_rs
metadata_tree <- data_list$metadata_tree
tree_color_order <- data_list$tree_color_order


### 辅助函数 ----------------------------------------------------------------
# 执行统计检验和标记字母分组
perform_statistical_analysis <- function(fin_df) {
  # 计算每组的最大值
  max_df <- fin_df %>%
    group_by(Group) %>%
    summarise(Value = max(Value)) %>%
    ungroup() %>%
    mutate(GroupMax = max(Value))
  
  # 计算每组的平均值
  mean_df <- fin_df %>%
    group_by(Group) %>%
    summarise(Value = mean(Value))
  
  # 执行成对t检验
  pttest_sig <- pairwise.t.test(fin_df$Value, fin_df$Group, 
                                p.adjust.method = "none", 
                                paired = FALSE)
  
  # 整理p值数据
  pttest_sig_df <- na.omit(reshape2::melt(pttest_sig$p.value)[c(2, 1, 3)])
  names(pttest_sig_df) <- c("group1", "group2", "p.adj")
  
  # 创建p值矩阵
  n <- length(unique(fin_df$Group))
  pvalue_df <- matrix(1, ncol = n, nrow = n)
  
  k <- 0
  for (i in 1:(n - 1)) { 
    for (j in (i + 1):n) { 
      k <- k + 1
      pvalue_df[i, j] <- pttest_sig_df$p.adj[k]
      pvalue_df[j, i] <- pttest_sig_df$p.adj[k]
    }
  }
  
  # 生成字母标记
  label_df <- agricolae::orderPvalue(
    mean_df$Group, 
    mean_df$Value, 
    0.05, 
    pvalue_df, 
    console = TRUE
  )
  label_df <- label_df[levels(fin_df$Group), ]
  max_df$Label <- label_df$groups
  
  return(list(
    max_df = max_df,
    sig_df = pttest_sig_df
  ))
}

# 绘制Alpha多样性图
plot_alpha_diversity <- function(data_df, max_df, tax_color_df, y_label_offset = 0.1) {
  p <- ggplot(data_df, aes(x = Group, y = Value)) + 
    # 数据点
    geom_point(
      aes(color = Group), 
      position = position_jitterdodge(dodge.width = 0.6), 
      alpha = 0.4, 
      size = 3, 
      stroke = 0
    ) +
    # 箱线图
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    # 小提琴图
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    # 统计标记
    geom_text(
      data = max_df, 
      aes(x = Group, y = Value + y_label_offset, label = Label),
      position = position_dodge(0.9), 
      size = 7 / 2.835
    ) + 
    # 标签和颜色
    labs(x = "", y = "Shannon index") + 
    scale_color_manual(values = tax_color_df$Color2) +
    # 分面
    facet_grid(Kingdom ~ ., scales = "free") +
    # 主题设置
    theme_bw() + 
    theme(
      plot.title = element_text(size = 7, color = "black", hjust = 0.5), 
      plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5), 
      axis.title = element_text(size = 7, color = "black"), 
      axis.text = element_text(size = 6, color = "black"), 
      legend.title = element_text(size = 7, color = "black"), 
      legend.text = element_text(size = 6, color = "black"), 
      strip.text = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(), 
      legend.position = "none"
    )
  
  return(p)
}


### 主分析流程 ----------------------------------------------------------------
analyze_alpha_diversity <- function(data_type, kingdoms, metadata) {
  # 初始化存储数据框
  all_data_df <- data.frame()
  all_max_df <- data.frame()
  all_sig_df <- data.frame()
  factor_level <- c()
  
  for (kingdom in kingdoms) {
    cat("Processing:", kingdom, "for", data_type, "\n")
    
    # 1. 数据加载
    if (data_type == "Amplicon") {
      # 扩增子数据
      raw_df <- read.csv(
        paste0("../asv_qa_alpha_beta/01-get_absolute_alpha_diversity/",
               kingdom, "_alpha_diversity_absolute.csv"), 
        fileEncoding = "GBK", 
        row.names = 1
      )
      names(raw_df) <- c("faith_pd", "observed_features", "Chao1", 
                         "shannon_entropy", "Simpson", "pielou_evenness")
      
      merge_df <- merge(raw_df, metadata, by.x = "row.names", by.y = "FileID")
      names(merge_df)[1] <- "Sample"
      
      res_df <- merge_df[c("Order", "shannon_entropy", "Sample")]
      names(res_df) <- c("Group", "Value", "Sample")
      
      kin_text <- ifelse(
        kingdom == "Protist", 
        paste0(data_type, " (18S)"), 
        paste0(data_type, " (", kingdom, ")")
      )
      
    } else {
      # 宏基因组数据
      raw_df <- read.csv(
        paste0("../ko_qa_alpha_beta/diversity/", 
               kingdom, ".quantitative.diversity.index.csv"), 
        row.names = 1
      )
      names(raw_df) <- c("shannon_entropy", "simpson", "observed_KOs", "pielou_evenness")
      row.names(raw_df) <- gsub("RS-", "", row.names(raw_df))
      
      merge_df <- merge(raw_df, metadata, by.x = "row.names", by.y = "TreeID")
      names(merge_df)[1] <- "Tree"
      
      res_df <- merge_df[c("Order", "shannon_entropy", "Tree")]
      names(res_df) <- c("Group", "Value", "Tree")
      
      kin_text <- paste0(data_type, " (", kingdom, ")")
    }
    
    # 2. 数据清洗和处理
    fin_df <- res_df[res_df$Group != "" & !is.na(res_df$Group), ]
    
    # 处理分类组
    fin_df$Group <- ifelse(
      fin_df$Group %in% tree_color_order[, 1], 
      fin_df$Group, 
      tree_color_order[, 1][nrow(tree_color_order)]
    )
    fin_df$Group <- factor(fin_df$Group, levels = tree_color_order[, 1])
    
    # 3. 统计分析
    stats_result <- perform_statistical_analysis(fin_df)
    
    # 4. 添加元信息
    fin_df$Kingdom <- kin_text
    stats_result$max_df$Kingdom <- kin_text
    stats_result$sig_df$Kingdom <- kin_text
    
    # 5. 保存结果
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, stats_result$max_df)
    all_sig_df <- rbind(all_sig_df, stats_result$sig_df)
    factor_level <- c(factor_level, kin_text)
  }
  
  # 6. 因子水平设置
  all_data_df$Kingdom <- factor(all_data_df$Kingdom, levels = factor_level)
  all_max_df$Kingdom <- factor(all_max_df$Kingdom, levels = factor_level)
  
  return(list(
    all_data_df = all_data_df,
    all_max_df = all_max_df,
    all_sig_df = all_sig_df
  ))
}


### 执行分析 ----------------------------------------------------------------
# 1. 扩增子数据分析
amplicon_kingdoms <- c("16S", "ITS", "Protist")
amplicon_results <- analyze_alpha_diversity(
  data_type = "Amplicon",
  kingdoms = amplicon_kingdoms,
  metadata = metadata_rs
)

# 保存结果
write.csv(
  amplicon_results$all_sig_df,
  paste0(dir_name, "/alpha_diversity_Order_Amplicon_sig.csv"),
  row.names = FALSE
)
write.csv(
  amplicon_results$all_data_df,
  paste0(dir_name, "/alpha_diversity_Order_Amplicon_data.csv"),
  row.names = FALSE
)
write.csv(
  amplicon_results$all_max_df,
  paste0(dir_name, "/alpha_diversity_Order_Amplicon_mean.csv"),
  row.names = FALSE
)

# 绘制扩增子图
p1 <- plot_alpha_diversity(
  data_df = amplicon_results$all_data_df,
  max_df = amplicon_results$all_max_df,
  tax_color_df = tree_color_order,
  y_label_offset = 0.5
)


# 2. 宏基因组数据分析
metagenome_kingdoms <- c("Bacteria", "Fungi", "Protist")
metagenome_results <- analyze_alpha_diversity(
  data_type = "Metagenome",
  kingdoms = metagenome_kingdoms,
  metadata = metadata_tree
)

# 保存结果
write.csv(
  metagenome_results$all_sig_df,
  paste0(dir_name, "/alpha_diversity_Order_Metagenome_sig.csv"),
  row.names = FALSE
)
write.csv(
  metagenome_results$all_data_df,
  paste0(dir_name, "/alpha_diversity_Order_Metagenome_data.csv"),
  row.names = FALSE
)
write.csv(
  metagenome_results$all_max_df,
  paste0(dir_name, "/alpha_diversity_Order_Metagenome_mean.csv"),
  row.names = FALSE
)

# 绘制宏基因组图
p2 <- plot_alpha_diversity(
  data_df = metagenome_results$all_data_df,
  max_df = metagenome_results$all_max_df,
  tax_color_df = tree_color_order,
  y_label_offset = 0.1
)


# 3. 合并图形并保存
combined_plot_width <- 14
combined_plot_height <- 21

combined_plot <- plot_grid(
  p1 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
  p2 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
  ncol = 2, 
  rel_widths = c(1, 1)
)

# 保存组合图
ggsave(
  paste0(dir_name, "/shannon_entropy_across_plant_orders.png"), 
  combined_plot, 
  width = combined_plot_width, 
  height = combined_plot_height, 
  dpi = 600, 
  units = "cm"
)
ggsave(
  paste0(dir_name, "/shannon_entropy_across_plant_orders.pdf"), 
  combined_plot, 
  width = combined_plot_width, 
  height = combined_plot_height, 
  units = "cm"
)


### 分类与功能多样性相关性分析 ------------------------------------------------
analyze_correlation <- function() {
  # 1. 处理扩增子数据
  amplicon_raw_df <- read.csv("results/alpha_diversity_Order_Amplicon_data.csv")
  amplicon_raw_df$Type <- "Amplicon"
  amplicon_raw_df$Kingdom <- ifelse(
    grepl("16S", amplicon_raw_df$Kingdom), "Bacteria", 
    ifelse(grepl("ITS", amplicon_raw_df$Kingdom), "Fungi", "Protists")
  )
  
  res_df <- merge(
    amplicon_raw_df, 
    metadata_rs[c("FileID", "TreeID")], 
    by.x = "Sample", 
    by.y = "FileID"
  )
  
  # 聚合扩增子数据
  aggregated_df <- res_df %>%
    group_by(Group, Kingdom, Type, TreeID) %>%
    summarise(
      Value = mean(Value),
      .groups = "drop"
    ) %>%
    rename(Tree = TreeID, Amplicon = Value)
  
  # 2. 处理宏基因组数据
  meta_df <- read.csv("results/alpha_diversity_Order_Metagenome_data.csv")
  meta_df$Type <- "Metagenome"
  meta_df$Kingdom <- ifelse(
    grepl("Bacteria", meta_df$Kingdom), "Bacteria", 
    ifelse(grepl("Fungi", meta_df$Kingdom), "Fungi", "Protists")
  )
  
  # 重命名宏基因组数据
  metagenome_df <- meta_df %>%
    select(Group, Kingdom, Tree, Value) %>%
    rename(Metagenome = Value)
  
  # 3. 合并数据
  common_trees <- intersect(aggregated_df$Tree, metagenome_df$Tree)
  
  correlation_data <- merge(
    aggregated_df %>% filter(Tree %in% common_trees),
    metagenome_df %>% filter(Tree %in% common_trees),
    by = c("Group", "Kingdom", "Tree")
  )
  
  # 设置因子水平
  correlation_data$Group <- factor(
    correlation_data$Group, 
    levels = c(tree_color_order$Order, "Other order")
  )
  
  return(correlation_data)
}

# 执行相关性分析
correlation_data <- analyze_correlation()

# 绘制相关性图
p3 <- ggplot(correlation_data, aes(x = Metagenome, y = Amplicon)) +
  geom_point(aes(color = Group), size = 1) +
  geom_smooth(
    method = "lm", 
    formula = y ~ x, 
    color = "#F8766D", 
    linewidth = 0.7
  ) +
  labs(
    x = "Shannon index of metagenomic KO functions",
    y = "Shannon index of microbial species"
  ) +
  facet_wrap(~ Kingdom, ncol = 1, scales = "free", strip.position = "right") +
  scale_color_manual(
    values = tree_color_order$Color2, 
    name = "Plant",
    guide = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 7, color = "black", hjust = 0.5), 
    plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5), 
    axis.title = element_text(size = 7, color = "black"), 
    axis.text = element_text(size = 6, color = "black"), 
    strip.text = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    legend.margin = margin(-0.1, 0, 0, 0, "cm"),
    legend.title = element_text(size = 7, color = "black"), 
    legend.text = element_text(size = 6, color = "black"), 
    legend.key.size = unit(0.25, "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom"
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~")),
    parse = TRUE,
    label.x = 0.95,
    label.y = 0.95,
    size = 7 / 2.835,
    p.digits = 3,
    small.p = TRUE
  )


# 4. 创建最终组合图
final_plot <- plot_grid(
  p1,
  p2, 
  p3,
  ncol = 3, 
  rel_widths = c(1, 1, 1)
)

final_width <- 17.5
final_height <- 15
figure_name <- paste0(dir_name, "/Figure6")

# 保存最终图
ggsave(
  paste0(dir_name, "/shannon_entropy_across_plant_orders_and_taxonomy_vs_function.png"), 
  final_plot, 
  width = final_width, 
  height = final_height, 
  dpi = 600, 
  units = "cm"
)
ggsave(
  paste0(dir_name, "/shannon_entropy_across_plant_orders_and_taxonomy_vs_function.pdf"), 
  final_plot, 
  width = final_width, 
  height = final_height, 
  units = "cm"
)




############################################################################################################
### 分类与功能多样性相关性分析（Spearman）-----------------------------------------------
analyze_correlation <- function() {
  # 1. 处理扩增子数据
  amplicon_raw_df <- read.csv("results/alpha_diversity_Order_Amplicon_data.csv")
  amplicon_raw_df$Type <- "Amplicon"
  amplicon_raw_df$Kingdom <- ifelse(
    grepl("16S", amplicon_raw_df$Kingdom), "Bacteria", 
    ifelse(grepl("ITS", amplicon_raw_df$Kingdom), "Fungi", "Protists")
  )
  
  res_df <- merge(
    amplicon_raw_df, 
    metadata_rs[c("FileID", "TreeID")], 
    by.x = "Sample", 
    by.y = "FileID"
  )
  
  # 聚合扩增子数据
  aggregated_df <- res_df %>%
    group_by(Group, Kingdom, Type, TreeID) %>%
    summarise(
      Value = mean(Value),
      .groups = "drop"
    ) %>%
    rename(Tree = TreeID, Amplicon = Value)
  
  # 2. 处理宏基因组数据
  meta_df <- read.csv("results/alpha_diversity_Order_Metagenome_data.csv")
  meta_df$Type <- "Metagenome"
  meta_df$Kingdom <- ifelse(
    grepl("Bacteria", meta_df$Kingdom), "Bacteria", 
    ifelse(grepl("Fungi", meta_df$Kingdom), "Fungi", "Protists")
  )
  
  # 重命名宏基因组数据
  metagenome_df <- meta_df %>%
    select(Group, Kingdom, Tree, Value) %>%
    rename(Metagenome = Value)
  
  # 3. 合并数据
  common_trees <- intersect(aggregated_df$Tree, metagenome_df$Tree)
  
  correlation_data <- merge(
    aggregated_df %>% filter(Tree %in% common_trees),
    metagenome_df %>% filter(Tree %in% common_trees),
    by = c("Group", "Kingdom", "Tree")
  )
  
  # 设置因子水平
  correlation_data$Group <- factor(
    correlation_data$Group, 
    levels = c(tree_color_order$Order, "Other order")
  )
  
  return(correlation_data)
}

# 执行相关性分析
correlation_data <- analyze_correlation()

# 计算每个Kingdom的Spearman相关系数和p值
spearman_stats <- correlation_data %>%
  group_by(Kingdom) %>%
  summarise(
    rho = cor(Amplicon, Metagenome, method = "spearman", use = "complete.obs"),
    p.value = cor.test(Amplicon, Metagenome, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    rho_label = sprintf("rho == %.3f", rho),
    p_label = ifelse(p.value < 0.001, "p < 0.001", 
                     sprintf("p == %.3f", p.value)),
    label = paste(rho_label, p_label, sep = "~~")
  )

# 添加统计标签到数据
correlation_data <- correlation_data %>%
  left_join(spearman_stats %>% select(Kingdom, label), by = "Kingdom")

# 创建自定义的标签位置（每个分面的右上角）
label_positions <- correlation_data %>%
  group_by(Kingdom) %>%
  summarise(
    x = max(Metagenome, na.rm = TRUE) * 0.95,
    y = max(Amplicon, na.rm = TRUE) * 0.95
  ) %>%
  left_join(spearman_stats %>% select(Kingdom, label), by = "Kingdom")

# 绘制相关性图（使用Spearman相关性）
p3 <- ggplot(correlation_data, aes(x = Metagenome, y = Amplicon)) +
  geom_point(aes(color = Group), size = 1) +
  geom_smooth(
    method = "lm", 
    formula = y ~ x, 
    color = "#F8766D", 
    linewidth = 0.7,
    se = TRUE
  ) +
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label),
    parse = TRUE,
    hjust = 1,
    vjust = 1,
    size = 6 / 2.835,
    color = "black"
  ) +
  labs(
    x = "Shannon index of metagenomic KO functions",
    y = "Shannon index of microbial species"
  ) +
  facet_wrap(~ Kingdom, ncol = 1, scales = "free", strip.position = "right") +
  scale_color_manual(
    values = tree_color_order$Color2, 
    name = "Plant",
    guide = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size = 7, color = "black", hjust = 0.5), 
    plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5), 
    axis.title = element_text(size = 7, color = "black"), 
    axis.text = element_text(size = 6, color = "black"), 
    strip.text = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    legend.margin = margin(-0.1, 0, 0, 0, "cm"),
    legend.title = element_text(size = 7, color = "black"), 
    legend.text = element_text(size = 6, color = "black"), 
    legend.key.size = unit(0.25, "cm"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom"
  )

# 保存Spearman相关性结果
write.csv(
  spearman_stats,
  paste0(dir_name, "/spearman_correlation_results.csv"),
  row.names = FALSE
)


# 4. 创建最终组合图
final_plot <- plot_grid(
  p1,
  p2, 
  p3,
  ncol = 3, 
  rel_widths = c(1, 1, 1)
)

final_width <- 17.5
final_height <- 15
figure_name <- paste0(dir_name, "/Figure6")

# 保存最终图
ggsave(
  paste0(dir_name, "/shannon_entropy_across_plant_orders_and_taxonomy_vs_function-spearman.png"), 
  final_plot, 
  width = final_width, 
  height = final_height, 
  dpi = 600, 
  units = "cm"
)
ggsave(
  paste0(dir_name, "/shannon_entropy_across_plant_orders_and_taxonomy_vs_function-spearman.pdf"), 
  final_plot, 
  width = final_width, 
  height = final_height, 
  units = "cm"
)

cat("Spearman相关性分析完成！\n")
cat("Spearman相关性结果已保存到：", paste0(dir_name, "/spearman_correlation_results.csv\n"))

