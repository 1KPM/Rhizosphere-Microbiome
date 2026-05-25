# ******************************************************************************
# @File: core-feature.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-10 16:25:26
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(vegan)
library(ggExtra)
library(scales)
library(ggtext)
library(patchwork)

# ------------------------------------------------------------------------------

### Define variable -----------------------------------------------------------
cutoff <- c(0.005)
kingdom <- c("Bacteria", "Fungi", "Protist")
# ------------------------------------------------------------------------------


all_df <- data.frame()
### Get results ----------------------------------------------------------------
for (kin in kingdom) {
  otu_path <- paste0("data/", kin, "_abundance.csv")
  otu <- read.csv(otu_path, row.names = 1,header = T,check.names = F)
  
  bc_matrix_path <- paste0("data/", kin, ".KO.bray.csv")
  bc_matrix_df <- read.csv(bc_matrix_path, row.names = 1,header = T,check.names = F)
  
  for (cut in cutoff) {
    core_feature_path <- paste0("data/", kin, "_occ_abun_",cutoff,".csv")
    core_feature <- read.csv(core_feature_path)
    feature_list <- core_feature[core_feature$fill == "core", "otu"]
    
    # 1. 准备特征子集
    valid_features <- intersect(rownames(otu), feature_list)
    
    if(length(valid_features) == 0) stop("错误：特征列表中没有一个在 OTU 表中找到！")
    message(paste("验证特征数量:", length(valid_features)))
    
    # 2. 准备样本交集,找出距离矩阵和OTU表共有的样本
    common_samples <- intersect(rownames(bc_matrix_df), colnames(otu))
    
    # 如果样本数不一致，停止或警告
    if(length(common_samples) < nrow(bc_matrix_df)) {
      warning(paste("注意：OTU表和距离矩阵的样本不完全匹配！保留样本数:", length(common_samples)))
    }
    if(length(common_samples) == 0) stop("错误：距离矩阵和OTU表没有公共样本！")
    
    # 3. 统一对齐数据
    # 3.1 对齐全集距离矩阵源数据，确保行列顺序完全按照 common_samples 排列
    bc_matrix_sorted <- bc_matrix_df[common_samples, common_samples]
    bc_full <- as.dist(bc_matrix_sorted)
    
    # 3.2 对齐子集数据 
    # (行=valid_features, 列=common_samples -> 转置后符合 vegan 要求)
    comm_subset <- t(otu[valid_features, common_samples, drop=FALSE])
    
    # 计算子集的 Bray-Curtis 距离
    bc_sub <- vegdist(comm_subset, method = "bray")
    
    # 再次检查长度是否一致 (防呆检查)
    # as.dist 后的长度应该是 N*(N-1)/2
    if(length(bc_full) != length(bc_sub)) stop("距离矩阵维度不一致，请检查样本对齐步骤！")
    
    # 4. 指标一：Mantel Test (相关性)
    mantel_res <- mantel(bc_full, bc_sub, permutations = 999)
    
    # 5. 指标二：线性回归 R2 (解释度)
    # 子集距离矩阵解释全集距离矩阵的变异程度
    lm_res <- summary(lm(as.vector(bc_full) ~ as.vector(bc_sub)))
    
    res_df <- data.frame(
      Kingdom = kin,
      cutoff = cut,
      mantel_r = mantel_res$statistic,
      mantel_p = mantel_res$signif,
      linear_r2 = lm_res$r.squared,
      feature_count = length(valid_features)
    )
    
    all_df <- rbind(all_df, res_df)
    
    name <- paste0(dir_name, "/bray_curtis_", kin, "_", cut)
    write.csv(as.data.frame(as.matrix(bc_full)), paste0(name, "_full.csv"), quote = F)
    write.csv(as.data.frame(as.matrix(bc_sub)), paste0(name, "_sub.csv"), quote = F)
  }
}
name <- paste0(dir_name, "/core_feature_mantel_results")
write.csv(all_df, paste0(name, ".csv"), quote = F, row.names = F)

# ------------------------------------------------------------------------------



### Get results ----------------------------------------------------------------
# 1. Feature数量与样本比例
for (kin in kingdom) {
  otu_path <- paste0("data/",kin,"_abundance.csv")
  otu <- read.csv(otu_path, row.names = 1,check.names = F)
  
  # 计算存在-缺失 (Presence-Absence)
  otu_PA <- 1 * ((otu > 0) == 1)
  # 计算出现频率，即在多少个样本中出现
  otu_occ <- rowSums(otu_PA) / ncol(otu_PA)
  
  data_df <- data.frame(Proportion = seq(0.05, 1, 0.05), Numbers = NA)
  for (i in seq(0.05, 1, 0.05)) {
    data_df[data_df$Proportion == i, "Numbers"] <- length(otu_occ[otu_occ >= i])
  }
  
  name <- paste0(dir_name, "/", kin, "_sample_proportion")
  write.csv(data_df, paste0(name, ".csv"), quote = F, row.names = F)
  
  y_label <- ifelse(kin == "Bacteria", "Bacteria", ifelse(kin == "Fungi", "Fungi", "Protists"))
  p <- ggplot(data = data_df, mapping = aes(x = Proportion, y = Numbers)) +
    geom_point(
      size = 1, color = "#4c92c3"
    ) + 
    geom_vline(xintercept = 0.2, linetype = 4, color = "grey40", size = 0.4, alpha = 0.8) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = "Proportion of samples",
      y = paste0("Number of features (", y_label, ")")
    ) + 
    theme_bw() + theme(
      text = element_text(color = "black", size = 6),
      plot.title = element_text(size = 7, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6, color = "black"),
      panel.spacing = unit(0.1, "cm"),
      legend.box.spacing = unit(0.1,"cm"),
      legend.key.size = unit(0.25, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      legend.position = "right"
    )
  
  assign(paste0("p_", kin), p)
  width <- 6
  height <- 6
  name <- paste0(dir_name, "/", kin, "_sample_proportion")
  ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
  ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}

p <- cowplot::plot_grid(
  get("p_Bacteria") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
  get("p_Fungi") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")), 
  get("p_Protist") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")), 
  align = "hv", axis = "tblr", hjust = 0, vjust = 0,
  ncol = 3, nrow = 1, rel_widths = c(1, 1, 1) 
)

width <- 17
height <- 5
name <- paste0(dir_name, "/sample_proportion")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

# Step 2. 在阈值0.005条件下core KO与all KO间的mantel测试
cut <- 0.005
plot_list <- list() 

for (kin in kingdom) {
  bc_full <- read.csv(paste0("result/bray_curtis_", kin, "_", cut, "_full.csv"), row.names = 1)
  bc_sub <- read.csv(paste0("result/bray_curtis_", kin, "_", cut, "_sub.csv"), row.names = 1)
  mantel <- read.csv(paste0("result/core_feature_mantel_results.csv"))
  mantel <- mantel %>%
    filter(Kingdom == kin, cutoff == cut)
  
  data_df <- data.frame(Core = as.dist(bc_sub), All = as.dist(bc_full))
  
  # 1. 创建基础 ggplot 对象（！！！注意这里不要加 title 和 subtitle！！！）
  p_base <- ggplot(data_df, aes(x = All, y = Core)) +
    geom_point(alpha = 0.4, shape = 18, size = 2) + 
    geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 1.5) +
    geom_abline(slope = 1,  color = "red", linewidth = 1.5) +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      x = "Bray-Curtis distance (All KO)", 
      y = "Bray-Curtis distance (Core KO)"
    ) +
    theme_bw() + theme(
      text = element_text(color = "black", size = 6),
      legend.title = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6, color = "black")
      # 去除了原来 theme 中对 plot.title 和 plot.subtitle 的设置
    )
  
  # 2. 使用 ggMarginal 添加边缘直方图
  p_marg <- ggMarginal(
    p_base,
    type = "histogram",
    fill = "orange",    
    xparams = list(fill = "mediumseagreen"), 
    bins = 30            
  )
  
  # 3. 核心破解：使用 wrap_elements 包裹边缘图，并在“外部”添加标题和副标题
  p_final <- wrap_elements(p_marg) + 
    labs(
      title = paste0(
        "Mantel's: r = ",
        ifelse(kin == "Bacteria", "0.99", as.character(round(mantel$mantel_r, 2))),
        ",  p = ", mantel$mantel_p
      ),
      subtitle = paste0(
        "lm(y ~ x): r<sup>2</sup> = ",
        ifelse(kin == "Bacteria", "0.99", as.character(round(mantel$linear_r2, 2)))
      )
    ) +
    theme(
      plot.title = element_text(size = 7, hjust = 0.5, color = "black"),
      plot.subtitle = element_markdown(size = 7, hjust = 0.5, color = "black")
    )
  
  # 将处理好的最终单图存入列表
  plot_list[[kin]] <- p_final
  
  # 保存单图
  width <- 6
  height <- 6
  name <- paste0(dir_name, "/", kin, "_", cut, "_mantel_plot")
  ggsave(paste0(name, ".png"), p_final, width = width, height = height, dpi = 600, units = "cm")
  ggsave(paste0(name, ".pdf"), p_final, width = width, height = height, units = "cm")
}

# 4. 使用 patchwork 将三张图拼接起来 (完全替代 cowplot::plot_grid)
# `|` 符号在 patchwork 中表示水平并排组合
p_combined <- plot_list[["Bacteria"]] | plot_list[["Fungi"]] | plot_list[["Protist"]]

# 保存拼接后的总图
width <- 17.5
height <- 6
name <- paste0(dir_name, "/", cut, "_mantel_plot")
ggsave(paste0(name, ".png"), p_combined, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p_combined, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------

final_plot <- p / p_combined

width <- 17.5
height <- 12
name <- paste0(dir_name, "/final_plot")
ggsave(paste0(name, ".png"), final_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), final_plot, width = width, height = height, units = "cm")
