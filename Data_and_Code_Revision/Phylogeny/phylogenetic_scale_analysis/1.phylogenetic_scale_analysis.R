### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

### 补充分析3：系统发育尺度分析 - 检验系统发育保守性在不同植物类群中的变化 （大于10的样本）
### 回答"系统发育信号是否在所有分类水平一致"

library(picante)
library(ggplot2)

### Import data ----------------------------------------------------------------
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

### Keep sample with tree_ITS_reads
metadata_rs <- metadata_rs[metadata_rs$TreeID %in% tree_file$tip.label, ]
match(tree_file$tip.label, unique(metadata_rs$TreeID)) ###说明有的树木样本没有扩增子测序数据

tree_dist <- as.data.frame(cophenetic.phylo(tree_file)) 
###cophenetic.phylo computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.

######################################################################################################
# 读取top10植物目信息
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]

ggtheme <- 
  theme(
    text = element_text(color = "black", size = 8),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
  )


# 对每个微生物组执行分析（需要先准备好距离矩阵）
# 以Bray-Curtis距离为例
bray_16S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_16S_bray_curtis_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
jaccard_16S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_16S_jaccard_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
weighted_16S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_16S_weighted_unifrac_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
unweighted_16S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_16S_unweighted_unifrac_tree_mean.csv", row.names = 1, fileEncoding = "GBK")

bray_ITS_mean <-read.csv(file = "../asv_qa_alpha_beta/rhizo_ITS_bray_curtis_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
jaccard_ITS_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_ITS_jaccard_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
weighted_ITS_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_ITS_weighted_unifrac_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
unweighted_ITS_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_ITS_unweighted_unifrac_tree_mean.csv", row.names = 1, fileEncoding = "GBK")

bray_18S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_18S_bray_curtis_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
jaccard_18S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_18S_jaccard_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
weighted_18S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_18S_weighted_unifrac_tree_mean.csv", row.names = 1, fileEncoding = "GBK")
unweighted_18S_mean <- read.csv(file = "../asv_qa_alpha_beta/rhizo_18S_unweighted_unifrac_tree_mean.csv", row.names = 1, fileEncoding = "GBK")


# 函数：在不同分类水平进行Mantel检验
perform_phylogenetic_scale_analysis <- function(microbe_dist_matrix, microbe_name) {
  # 获取共有的样本
  common_samples <- intersect(rownames(microbe_dist_matrix), intersect(metadata_tree$TreeID, tree_file$tip.label))
  
  # 准备结果数据框
  scale_results <- data.frame()
  
  # 在不同分类水平进行分析
  taxonomic_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
  
  for(level in taxonomic_levels) {
    # 获取该水平的分类信息
    taxon_groups <- split(metadata_tree$TreeID, metadata_tree[[level]])
    
    # 只分析包含足够样本的组（至少10个样本）
    taxon_groups <- taxon_groups[sapply(taxon_groups, length) >= 10]
    
    for(group_name in names(taxon_groups)) {
      group_samples <- taxon_groups[[group_name]]
      group_samples <- intersect(group_samples, common_samples)
      
      if(length(group_samples) >= 10) {
        # 提取子距离矩阵
        microbe_subset <- as.dist(microbe_dist_matrix[group_samples, group_samples])
        tree_subset <- as.dist(cophenetic.phylo(tree_file)[group_samples, group_samples])
        
        # Mantel检验
        mantel_result <- mantel(tree_subset, microbe_subset, method = "spearman", permutations = 999)
        
        # 保存结果
        result_row <- data.frame(
          Taxonomic_Level = level,
          Group = group_name,
          Microbe_Group = microbe_name,
          Sample_Size = length(group_samples),
          Mantel_r = mantel_result$statistic,
          P_value = mantel_result$signif
        )
        
        scale_results <- rbind(scale_results, result_row)
      }
    }
  }
  
  return(scale_results)
}


scale_16S <- perform_phylogenetic_scale_analysis(bray_16S_mean, "Bacteria")
scale_ITS <- perform_phylogenetic_scale_analysis(bray_ITS_mean, "Fungi")
scale_18S <- perform_phylogenetic_scale_analysis(bray_18S_mean, "Protists")

write.csv(rbind.data.frame(scale_16S, scale_ITS, scale_18S), file = "./stats/Mantel_tests_between_microbial_bray_curtis_and_plant_phylogenetic_relatedness_all_clade_above_10_trees_amplicons.csv")

# 可视化结果
# 绘制图钉图
scale_16S_plot <- ggplot(scale_16S, aes(x = Mantel_r, y = reorder(Group, Mantel_r))) +
  geom_segment(aes(x = 0, xend = Mantel_r, yend = Group), 
               color = "gray70", linewidth = 0.8) +
  geom_point(aes(color = Mantel_r), size = 4, alpha = 0.8) +
  scale_color_gradient(low = "steelblue", high = "firebrick") +
  labs(x = "Mantel's r", y = "Bacteria",
       title = "",
       color = "Mantel's r") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


name <- "./plots/Mantel_tests_between_microbial_community_similarity_and_plant_phylogenetic_relatedness_all_clade_above_10_trees-16S"
width <- 12
height <- 18
ggsave(paste0(name, ".png"), scale_16S_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), scale_16S_plot, width = width, height = height, units = "cm")


###############################################################################合并
library(ggplot2)
library(patchwork)
library(dplyr)

# 将三个数据集的Group转换为有序因子，使用16S的顺序
# 计算三个数据集的最大Mantel_r，用于统一x轴范围
max_r <- max(c(scale_16S$Mantel_r, scale_ITS$Mantel_r, scale_18S$Mantel_r), na.rm = TRUE)

# 创建自定义主题函数，避免重复代码
create_lollipop_plot <- function(data, title, show_y_axis = TRUE) {
  p <- ggplot(data, aes(x = Mantel_r, y = Group)) +
    geom_segment(aes(x = 0, xend = Mantel_r, yend = Group), 
                 color = "gray70", linewidth = 0.6) +
    geom_point(aes(color = Mantel_r), size = 2.5, alpha = 0.8) +
    scale_color_gradient(low = "#1E88E5", high = "#D81B60", 
                         limits = c(-0.3, max_r)) +
    labs(x = "Mantel's r", y = "", title = title) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      axis.text.x = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = c(-0.3, max_r * 1.1))
  
  # 设置y轴标签
  if (show_y_axis) {
    p <- p + theme(axis.text.y = element_text(size = 9))
  } else {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  return(p)
}

# 创建三个子图
plot_16S <- create_lollipop_plot(scale_16S, "Bacterial ASV", show_y_axis = TRUE)
plot_ITS <- create_lollipop_plot(scale_ITS, "Fungal ASV", show_y_axis = FALSE)
plot_18S <- create_lollipop_plot(scale_18S, "Protistan ASV", show_y_axis = FALSE)


# 使用patchwork组合三个图
combined_plot <- plot_16S + plot_ITS + plot_18S +
  plot_layout(nrow = 1, widths = c(1, 1, 1)) +  # 调整宽度比例
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  )

# 显示图形
print(combined_plot)

# 保存合并后的图
name <- "./plots/Mantel_tests_between_microbial_community_similarity_and_plant_phylogenetic_relatedness_all_clade_above_10_trees-combined"
width <- 24  # 增加宽度以容纳三个子图
height <- 16
ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")


########################################################################
###宏基因组数据
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


####
scale_bray_bac_KO <- perform_phylogenetic_scale_analysis(bray_bac_KO, "Bacteria")
scale_bray_fun_KO <- perform_phylogenetic_scale_analysis(bray_fun_KO, "Fungi")
scale_bray_pro_KO <- perform_phylogenetic_scale_analysis(bray_pro_KO, "Protists")

write.csv(rbind.data.frame(scale_bray_bac_KO, scale_bray_fun_KO, scale_bray_pro_KO), 
          file = "./stats/Mantel_tests_between_microbial_bray_curtis_and_plant_phylogenetic_relatedness_all_clade_above_10_trees_metagenomes.csv")

# 创建三个子图
plot_bray_bac_KO <- create_lollipop_plot(scale_bray_bac_KO, "Bacterial KO", show_y_axis = TRUE)
plot_bray_fun_KO <- create_lollipop_plot(scale_bray_fun_KO, "Fungal KO", show_y_axis = FALSE)
plot_bray_pro_KO <- create_lollipop_plot(scale_bray_pro_KO, "Protistan KO", show_y_axis = FALSE)


# 使用patchwork组合三个图
combined_plot <- plot_bray_bac_KO + plot_bray_fun_KO + plot_bray_pro_KO +
  plot_layout(nrow = 1, widths = c(1, 1, 1)) +  # 调整宽度比例
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  )

# 显示图形
print(combined_plot)

# 保存合并后的图
name <- "./plots/Mantel_tests_between_microbial_community_similarity_and_plant_phylogenetic_relatedness_all_clade_above_10_trees-combined-KO"
width <- 24  # 增加宽度以容纳三个子图
height <- 16
ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")









