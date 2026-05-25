# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# 加载必要包
library(dplyr)
library(tidyr)
library(rstatix)
library(agricolae)
library(stats)
library(ggplot2)
library(patchwork)
library(ape)

# 设置主题
ggtheme2 <- 
  theme(
    text = element_text(color = "black", size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, color = "black"),
    strip.text = element_text(color = "black", size = 8, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
  )

# 读取数据
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")

# 统计函数
calculate_group_stats <- function(data, value_col, group_col, alpha = 0.05, p_adjust_method = "fdr") {
  names(data)[names(data) == value_col] <- "Dist"
  names(data)[names(data) == group_col] <- "Group"
  
  if (!is.factor(data$Group)) {
    data$Group <- factor(data$Group)
  }
  
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
  
  kruskal_test <- kruskal.test(Dist ~ Group, data = data)
  sig_value <- kruskal_test$p.value
  dunn_test <- rstatix::dunn_test(data, Dist ~ Group, p.adjust.method = p_adjust_method)
  
  pvalue_df <- data.frame(
    KruskalWallis = sig_value, 
    dunn_test[c("group1", "group2", "p.adj")]
  )
  pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                            ifelse(pvalue_df$p.adj < 0.01, "**", 
                                   ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
  
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
  
  letter_df <- agricolae::orderPvalue(
    summarise_df$Group, 
    summarise_df$mean, 
    alpha, 
    pvalue_matrix, 
    console = FALSE
  )
  
  letter_df <- letter_df[levels(data$Group), ]
  summarise_df$label <- letter_df$groups
  
  return(list(summary_df = summarise_df))
}

# 16S数据
alpha_16S <- read.csv("../asv_qa_alpha_beta/01-get_absolute_alpha_diversity/16S_alpha_diversity_absolute.csv", fileEncoding = "GBK")
alpha_16S <- merge.data.frame(metadata_rs[,c(1,2,10)], alpha_16S, by.x = "FileID", by.y = "X")
names(alpha_16S)[4:9] <- c("faith_pd", "observed_features", "Chao1", "shannon_entropy", "Simpson", "pielou_evenness")
alpha_16S$Group <- ifelse(alpha_16S$Order %in% top10orders, alpha_16S$Order, "Others")
alpha_16S$Group <- factor(alpha_16S$Group, levels = c(top10orders, "Others"))

# ITS数据
alpha_ITS <- read.csv("../asv_qa_alpha_beta/01-get_absolute_alpha_diversity/ITS_alpha_diversity_absolute.csv", fileEncoding = "GBK")
alpha_ITS <- merge.data.frame(metadata_rs[,c(1,2,10)], alpha_ITS, by.x = "FileID", by.y = "X")
names(alpha_ITS)[4:9] <- c("faith_pd", "observed_features", "Chao1", "shannon_entropy", "Simpson", "pielou_evenness")
alpha_ITS$Group <- ifelse(alpha_ITS$Order %in% top10orders, alpha_ITS$Order, "Others")
alpha_ITS$Group <- factor(alpha_ITS$Group, levels = c(top10orders, "Others"))

# 18S数据
alpha_18S <- read.csv("../asv_qa_alpha_beta/01-get_absolute_alpha_diversity/Protist_alpha_diversity_absolute.csv", fileEncoding = "GBK")
alpha_18S <- merge.data.frame(metadata_rs[,c(1,2,10)], alpha_18S, by.x = "FileID", by.y = "X")
names(alpha_18S)[4:9] <- c("faith_pd", "observed_features", "Chao1", "shannon_entropy", "Simpson", "pielou_evenness")
alpha_18S$Group <- ifelse(alpha_18S$Order %in% top10orders, alpha_18S$Order, "Others")
alpha_18S$Group <- factor(alpha_18S$Group, levels = c(top10orders, "Others"))

# 统计计算
alpha_16S_faith_pd_summarise <- calculate_group_stats(data = alpha_16S, value_col = "faith_pd", group_col = "Group")
alpha_ITS_faith_pd_summarise <- calculate_group_stats(data = alpha_ITS, value_col = "faith_pd", group_col = "Group")
alpha_18S_faith_pd_summarise <- calculate_group_stats(data = alpha_18S, value_col = "faith_pd", group_col = "Group")

# 生成组合图
regions <- c("16S", "ITS", "18S")
plot_list <- list()

for (i in seq_along(regions)) {
  region <- regions[i]
  alpha_data <- get(paste0("alpha_", region))
  summary_data <- get(paste0("alpha_", region, "_faith_pd_summarise"))$summary_df
  
  alpha_data$Group <- factor(alpha_data$Group, levels = rev(levels(alpha_data$Group)))
  summary_data$Group <- factor(summary_data$Group, levels = rev(levels(summary_data$Group)))
  
  if (region == "16S") {
    y_label <- "Bacteria"
    y_limit <- 400
    text_y_position <- 400
  } else if (region == "ITS") {
    y_label <- "Fungi"
    y_limit <- 600
    text_y_position <- 600
  } else {
    y_label <- "Protists"
    y_limit <- 300
    text_y_position <- 300
  }
  
  p <- ggplot(data = alpha_data, mapping = aes(x = Group, y = faith_pd)) +
    geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
    geom_text(
      data = summary_data,
      mapping = aes(x = Group, y = text_y_position, label = label),
      position = position_dodge(0.9),
      size = 8 / 2.835
    ) +
    coord_flip() +
    labs(x = NULL, y = y_label) + 
    scale_color_manual(values = rev(top_tree$Color2)) +
    scale_y_continuous(limits = c(0, y_limit)) +
    theme_classic() + 
    ggtheme2
  
  if (region != "16S") {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  plot_list[[region]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(widths = c(1, 1, 1), heights = 1, guides = "collect") &
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"), panel.spacing = unit(2, "mm"))

# 保存图形
plot_dir <- "./plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
name <- file.path(plot_dir, "combined_Faiths_PD_top10orders")
ggsave(paste0(name, ".png"), combined_plot, width = 10, height = 9, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = 10, height = 9, units = "cm")
cat("组合图形生成完成！保存位置：", plot_dir, "\n")
print(combined_plot)