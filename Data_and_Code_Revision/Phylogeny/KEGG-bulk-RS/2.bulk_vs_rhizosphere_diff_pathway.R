# 设置工作路径
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# 加载包
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

# 定义分组
kingdoms <- c("Bacteria", "Fungi", "Protist")

# 为每个kingdom执行差异分析
for (kingdom in kingdoms) {
  cat("正在处理", kingdom, "...\n")
  
  # 导入数据
  bs_df_path <- paste0("./data/", kingdom, "_bulk_tpm_pathway_abundance.csv")
  rs_df_path <- paste0("./data/", kingdom, "_RS_tpm_pathway_abundance.csv")
  
  bs_df <- read.csv(bs_df_path, check.names = F, row.names = 1)
  rs_df <- read.csv(rs_df_path, check.names = F, row.names = 1)
  
  # 移除Unassign
  bs_df <- bs_df[row.names(bs_df) != "Unassign", ]
  rs_df <- rs_df[row.names(rs_df) != "Unassign", ]
  
  # 数据预处理 - 匹配样本
  bs_core <- gsub("^S-", "", colnames(bs_df))
  names(bs_core) <- colnames(bs_df)
  rs_core <- gsub("^RS-", "", colnames(rs_df))
  names(rs_core) <- colnames(rs_df)
  common_ids <- intersect(bs_core, rs_core)
  
  bs_filtered <- bs_df[, names(bs_core)[bs_core %in% common_ids], drop = FALSE]
  rs_filtered <- rs_df[, names(rs_core)[rs_core %in% common_ids], drop = FALSE]
  
  # 确保行名一致
  all_rows <- union(rownames(bs_filtered), rownames(rs_filtered))
  bs_complete <- matrix(0, nrow = length(all_rows), ncol = ncol(bs_filtered),
                        dimnames = list(all_rows, colnames(bs_filtered)))
  bs_complete[rownames(bs_filtered), ] <- as.matrix(bs_filtered)
  bs_complete <- as.data.frame(bs_complete)
  
  rs_complete <- matrix(0, nrow = length(all_rows), ncol = ncol(rs_filtered),
                        dimnames = list(all_rows, colnames(rs_filtered)))
  rs_complete[rownames(rs_filtered), ] <- as.matrix(rs_filtered)
  rs_complete <- as.data.frame(rs_complete)
  
  # 进行配对Wilcoxon秩和检验
  cat("  进行配对Wilcoxon检验...\n")
  wilcoxon_results <- data.frame(
    Pathway = rownames(bs_complete),
    p_value = NA,
    adjusted_p = NA,
    log2FC = NA,
    mean_bs = NA,
    mean_rs = NA
  )
  
  for (i in 1:nrow(bs_complete)) {
    pathway <- rownames(bs_complete)[i]
    bs_values <- as.numeric(bs_complete[i, ])
    rs_values <- as.numeric(rs_complete[i, ])
    
    # 计算log2 Fold Change
    mean_bs <- mean(bs_values, na.rm = TRUE)
    mean_rs <- mean(rs_values, na.rm = TRUE)
    
    # 避免除0错误
    if (mean_bs == 0) mean_bs <- 0.001
    if (mean_rs == 0) mean_rs <- 0.001
    
    log2fc <- log2(mean_rs / mean_bs)
    
    # 执行配对Wilcoxon检验
    if (all(bs_values == 0) && all(rs_values == 0)) {
      p_val <- 1
    } else {
      test_result <- wilcox.test(bs_values, rs_values, paired = TRUE, exact = FALSE)
      p_val <- test_result$p.value
    }
    
    wilcoxon_results[i, "p_value"] <- p_val
    wilcoxon_results[i, "log2FC"] <- log2fc
    wilcoxon_results[i, "mean_bs"] <- mean_bs
    wilcoxon_results[i, "mean_rs"] <- mean_rs
  }
  
  # 校正p值
  wilcoxon_results$adjusted_p <- p.adjust(wilcoxon_results$p_value, method = "fdr")
  
  # 定义显著性阈值
  wilcoxon_results$Significance <- "Not significant"
  wilcoxon_results$Significance[wilcoxon_results$adjusted_p < 0.05 & 
                                  abs(wilcoxon_results$log2FC) > 0] <- "Significant"
  wilcoxon_results$Significance[wilcoxon_results$adjusted_p < 0.05 & 
                                  wilcoxon_results$log2FC > 0] <- "Up in RS"
  wilcoxon_results$Significance[wilcoxon_results$adjusted_p < 0.05 & 
                                  wilcoxon_results$log2FC < 0] <- "Up in BS"
  
  # 标记top显著通路用于标注
  sig_results <- wilcoxon_results[wilcoxon_results$Significance != "Not significant", ]
  if (nrow(sig_results) > 0) {
    # 按p值和fold change排序
    sig_results$score <- -log10(sig_results$adjusted_p) * abs(sig_results$log2FC)
    sig_results <- sig_results[order(sig_results$score, decreasing = TRUE), ]
    
    # 选择top 20标注
    top_n <- min(20, nrow(sig_results))
    top_pathways <- head(sig_results$Pathway, top_n)
    wilcoxon_results$Label <- ifelse(wilcoxon_results$Pathway %in% top_pathways, 
                                     wilcoxon_results$Pathway, "")
  } else {
    wilcoxon_results$Label <- ""
  }
  
  # 保存结果
  write.csv(wilcoxon_results, 
            file = paste0(kingdom, "_BS_RS_pathway_wilcoxon_results.csv"), 
            row.names = FALSE)
  
  # 创建火山图
  cat("  绘制火山图...\n")
  volcano_plot <- ggplot(wilcoxon_results, aes(x = log2FC, y = -log10(adjusted_p))) +
    geom_point(aes(color = Significance, size = abs(log2FC)), alpha = 0.7) +
    scale_color_manual(values = c(
      "Not significant" = "gray70",
      "Up in RS" = "red2",
      "Up in BS" = "blue2",
      "Significant" = "orange"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.7) +
    labs(
      x = paste(kingdom, "log2 Fold Change (RS/BS)"),
      y = "-log10(Adjusted p-value)",
      color = "Significance"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.box = "vertical"
    ) +
    guides(size = "none")
  
  # 添加标签
  if (any(wilcoxon_results$Label != "")) {
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = wilcoxon_results[wilcoxon_results$Label != "", ],
        aes(label = Label),
        size = 1,
        max.overlaps = 20,
        box.padding = 0.5,
        segment.color = "gray50",
        segment.alpha = 0.5
      )
  }
  
  # 保存火山图
  ggsave(paste0(kingdom, "_volcano_plot.png"), 
         volcano_plot, 
         width = 7, 
         height = 8, 
         dpi = 600,
         units = "cm")
  ggsave(paste0(kingdom, "_volcano_plot.pdf"), 
         volcano_plot, 
         width = 7, 
         height = 8,
         units = "cm")
  
  # 创建结果摘要
  cat("\n=== ", kingdom, " 结果摘要 ===\n")
  cat("总通路数:", nrow(wilcoxon_results), "\n")
  cat("显著差异通路 (adjusted p < 0.05 & |log2FC| > 1):", 
      sum(wilcoxon_results$Significance %in% c("Up in RS", "Up in BS")), "\n")
  cat("  - RS中上调通路:", sum(wilcoxon_results$Significance == "Up in RS"), "\n")
  cat("  - BS中上调通路:", sum(wilcoxon_results$Significance == "Up in BS"), "\n")
  
  # 显示top显著通路
  sig_up_rs <- wilcoxon_results[wilcoxon_results$Significance == "Up in RS", ]
  sig_up_bs <- wilcoxon_results[wilcoxon_results$Significance == "Up in BS", ]
  
  if (nrow(sig_up_rs) > 0) {
    cat("\nTop 5 RS中上调通路:\n")
    print(head(sig_up_rs[order(sig_up_rs$log2FC, decreasing = TRUE), 
                         c("Pathway", "log2FC", "adjusted_p")], 5))
  }
  
  if (nrow(sig_up_bs) > 0) {
    cat("\nTop 5 BS中上调通路:\n")
    print(head(sig_up_bs[order(sig_up_bs$log2FC), 
                         c("Pathway", "log2FC", "adjusted_p")], 5))
  }
  cat("\n")
}

cat("分析完成！\n")