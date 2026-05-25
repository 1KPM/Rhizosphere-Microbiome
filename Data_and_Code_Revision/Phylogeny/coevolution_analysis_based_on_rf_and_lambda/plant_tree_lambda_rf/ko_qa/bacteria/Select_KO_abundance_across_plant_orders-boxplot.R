### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
# 加载必要的包
library(dplyr)
library(tidyr)
library(rstatix)
library(agricolae)
library(stats)
library(ggplot2)
library(ggtree)
library(ape)
library(geosphere)

#########################################################
#### 在Plant tree sig KO barplot的基础上进一步统计可视化
####循环可视化每个KO的组间差异boxplot

###在这基础上进一步可视化7个植物目的代表性KO
# KOID	Kingdom	lambda	logL	logL0	P	padj	KOLabel	gene	Description	Tree_Class	rf_FDR_Class	Tree_Order	rf_FDR_Order	Tree_Family	rf_FDR_Family
# K14491	Bacteria	0.300388696	-1082.955038	-1111.902425	2.77E-14	2.15E-11	K14491|ARR-B	ARR-B	two-component response regulator ARR-B family	NA	NA	Rosales	3.37E-07	Moraceae	3.12E-11
# K14658	Bacteria	0.143856173	-937.1654584	-954.2102675	5.26E-09	4.85E-07	K14658|nodA	nodA	nodulation protein A [EC:2.3.1.-]	NA	NA	Fabales	0.009238176	Fabaceae	0.006479555
# K02640	Bacteria	0.529077505	-1164.188035	-1209.644791	1.50E-21	1.87E-17	K02640|petG	petG	cytochrome b6-f complex subunit 5	Magnoliopsida	2.95E-06	Gentianales	0.000164682	Apocynaceae	2.95E-12
# K06972	Bacteria	0.103627772	-904.046762	-911.1944249	0.000156252	0.001057571	K06972|PITRM1, PreP, CYM1	PITRM1, PreP, CYM1	presequence protease [EC:3.4.24.-]	NA	NA	Lamiales	0.001620885	Acanthaceae	0.005973177
# K16261	Bacteria	0.216381059	-793.0485293	-804.3985435	1.89E-06	4.69E-05	K16261|YAT	YAT	yeast amino acid transporter	NA	NA	Asparagales	0.000795494	Asparagaceae	0.000243796
# K19110	Bacteria	0.17238518	-761.1655003	-774.682993	2.00E-07	8.86E-06	K19110|tblC	tblC	putative 2-oxoglutarate oxygenase	Others	0.027753053	Myrtales	0.000186137	Myrtaceae	9.36E-05
# K17413	Bacteria	0.409139082	-1447.813161	-1475.343746	1.17E-13	6.93E-11	K17413|MRPS35	MRPS35	small subunit ribosomal protein S35	Liliopsida	0.017674094	Arecales	2.71E-11	Arecaceae	5.30E-10


### Import data ----------------------------------------------------------------
tree_file <- read.tree("../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../../../../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../../../../metadata/tree_metadata_merge_info.csv")

### Keep sample with tree_ITS_reads
metadata_rs <- metadata_rs[metadata_rs$TreeID %in% tree_file$tip.label, ]
match(tree_file$tip.label, unique(metadata_rs$TreeID)) ###说明有的树木样本没有扩增子测序数据

tree_dist <- as.data.frame(cophenetic.phylo(tree_file)) 
###cophenetic.phylo computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.

######################################################################################################
# 读取top10植物目信息
top_tree <- read.csv("../../../../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]


####(1) select KO
select_KO <- read.csv(file = "./stats/bac_order_coevoled_top42_ko_qa.csv", fileEncoding = "GBK")
row.names(select_KO) <- select_KO$KOID
select_KO <- as.data.frame(t(select_KO[,-c(1:3)]))
select_KO$TreeID <- row.names(select_KO)

select_KO <- merge.data.frame(metadata_tree[,c(1,8)], select_KO, by = "TreeID")
select_KO$Group <- ifelse(select_KO$Order %in% top10orders, select_KO$Order, "Others")
select_KO$Group <- factor(select_KO$Group, levels = c(top10orders, "Others"))

select_KO_info <- read.csv(file = "./stats/bac_order_coevoled_top42_ko_lambda_rf.csv", fileEncoding = "GBK")

###########################################
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

### 计算组间丰度差异统计和多重比较结果
# 测试 data = select_KO; value_col = "K01369"; group_col = "Group"; alpha = 0.05; p_adjust_method = "fdr"

calculate_group_stats <- function(data, value_col, group_col, alpha = 0.05, p_adjust_method = "fdr") {
  
  # 重命名列以便于处理
  names(data)[names(data) == value_col] <- "Abundance"
  names(data)[names(data) == group_col] <- "Group"
  
  # 确保Group是因子
  if (!is.factor(data$Group)) {
    data$Group <- factor(data$Group)
  }
  
  ###########################################
  # 1. 计算每组的基本统计量
  ###########################################
  summarise_df <- data %>%
    group_by(Group) %>%
    summarise(
      n = length(Abundance),
      mean = mean(Abundance, na.rm = TRUE),
      min = min(Abundance, na.rm = TRUE),
      max = max(Abundance, na.rm = TRUE),
      sd = sd(Abundance, na.rm = TRUE),
      se = sd / sqrt(n),
      ci = qt(1 - alpha/2, df = n - 1) * se
    ) %>%
    mutate(allmax = max(max))
  
  #########################################
  # 2. 非参数检验
  ###########################################
  # Kruskal-Wallis检验
  kruskal_test <- kruskal.test(Abundance ~ Group, data = data)
  sig_value <- kruskal_test$p.value
  
  # Dunn事后检验
  dunn_test <- rstatix::dunn_test(data, Abundance ~ Group, p.adjust.method = p_adjust_method)
  
  # 添加显著性标记
  pvalue_df <- data.frame(
    KruskalWallis = sig_value, 
    dunn_test[c("group1", "group2", "p.adj")]
  )
  pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                            ifelse(pvalue_df$p.adj < 0.01, "**", 
                                   ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
  
  ###########################################
  # 3. 生成p值矩阵
  ###########################################
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
  
  ###########################################
  # 4. 生成显著性字母标记
  ###########################################
  # 使用agricolae包
  letter_df <- agricolae::orderPvalue(
    summarise_df$Group, 
    summarise_df$mean, 
    alpha, 
    pvalue_matrix, 
    console = FALSE
  )
  
  # 确保字母标记的顺序与summarise_df一致
  letter_df <- letter_df[levels(data$Group), ]
  summarise_df$label <- letter_df$groups
  
  ###########################################
  # 5. 返回结果列表
  ###########################################
  result_list <- list(
    summary_df = summarise_df,
    kruskal_test = kruskal_test,
    dunn_test = dunn_test,
    letter_df = letter_df,
    pvalue_matrix = pvalue_matrix
  )
  
  return(result_list)
}



###统计abundance差异
names(select_KO)
select_KO_abundance_summarise <- calculate_group_stats(data = select_KO, value_col = "K01369", group_col = "Group")
select_KO_abundance_summarise$summary_df
#####################################
###循环出图
# 定义要循环的因子
ko_ids <- names(select_KO)[3:44]

# 设置保存目录
plot_dir <- "./boxplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# 循环生成所有图形

for (ko in ko_ids) {
    
    # 获取对应的数据框
    select_KO_abundance_summarise <- calculate_group_stats(data = select_KO, value_col = ko, group_col = "Group")
    summary_data <- select_KO_abundance_summarise$summary_df
    
    ko_info <- select_KO_info[select_KO_info$KOID == ko, c("KOID", "Description")]
    ko_lambda <- select_KO_info[select_KO_info$KOID == ko, c("lambda", "padj")]
    
  
    # 创建图形
    p <- ggplot(data = select_KO, mapping = aes(x = Group, y = .data[[ko]])) +
      geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
      geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
      geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
      geom_text(
        data = summary_data,
        mapping = aes(x = Group, y = max + 10, label = label),
        position = position_dodge(0.9),
        size = 8 / 2.835
      ) +
      
      labs(
        x = NULL,
        y = "TPM",
        title = paste0(ko_info$KOID[1], " | ", ko_info$Description[1]),
        subtitle = paste0(
          "Pagel's Lambda = ", format(ko_lambda$lambda[1], scientific = F, digits = 3),
          ", FDR = ", format(ko_lambda$padj[1], scientific = TRUE, digits = 3), "\n",
          "Kruskal-Wallis P = ", format(select_KO_abundance_summarise$kruskal_test$p.value, 
                                        scientific = TRUE, digits = 3)
        )
      ) + 
      scale_color_manual(values = top_tree$Color2) +
      theme_bw() + 
      ggtheme
    
    # 设置保存参数
    name <- file.path(plot_dir, paste0(ko, "_abundance_across_top10orders"))
    width <- 8
    height <- 6
    
    # 保存图形
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

}



###############################################################################
###可视化7个目的代表性KO
# K14491
# K14658
# K02640
# K06972
# K16261
# K19110
# K17413

# 定义要循环的因子
library("patchwork")

ggtheme2 <- 
  theme(
    text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 6, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6, color = "black"),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(color = "black", size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
  )


###############################################################################
### 可视化7个目的代表性KO
# 定义要循环的KO列表
rep_kos <- c("K14658", "K14491", "K06972", "K02640",  
             "K16261", "K19110", "K17413")

# 设置保存目录
plot_dir <- "./boxplots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# 创建存储图形的列表
plot_list <- list()

# 定义格式化函数
format_p_value <- function(p) {
  if (p < 0.001) {
    return("< 0.001")
  } else {
    return(format(p, scientific = FALSE, digits = 3))
  }
}

# 循环生成所有图形
for (i in seq_along(rep_kos)) {
  ko <- rep_kos[i]
  
  # 获取对应的数据框
  select_KO_abundance_summarise <- calculate_group_stats(data = select_KO, value_col = ko, group_col = "Group")
  summary_data <- select_KO_abundance_summarise$summary_df
  
  # 获取KO信息
  ko_info <- select_KO_info[select_KO_info$KOID == ko, c("KOID", "KOLabel", "Description")]
  ko_lambda <- select_KO_info[select_KO_info$KOID == ko, c("lambda", "padj")]
  
  # 创建y轴标签（包含描述信息）
  y_label <- if (nrow(ko_info) > 0) {
    paste0(ko_info$KOID[1], ": ", ko_info$Description[1])
  } else {
    ko
  }
  
  # 创建图形（使用coord_flip旋转坐标轴）
  p <- ggplot(data = select_KO, mapping = aes(x = Group, y = log(.data[[ko]]+1))) +
    geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
    geom_text(
      data = summary_data,
      mapping = aes(x = Group, y = log(max*2), label = label),  # 使用固定的偏移量
      position = position_dodge(0.9),
      size = 6 / 2.835
    ) +
    coord_flip() +  # 旋转坐标轴
    labs(
      x = NULL,
      y = paste0(ko_info$KOLabel[1], "\n",
                 "(λ=", round(ko_lambda$lambda[1], 2),
                 ", P=", format(ko_lambda$padj[1], scientific = TRUE, digits = 2), ")"),
      title = NULL
    ) + 
    scale_color_manual(values = top_tree$Color2) +
    theme_bw() + 
    ggtheme2
  
  # 只对第一个图显示Y轴标签
  if (i != 1) {
    p <- p + theme(axis.text.y = element_blank(),  # 隐藏y轴文本
                   axis.ticks.y = element_blank())  # 隐藏y轴刻度
  } else {
    # 第一个图保持完整的Y轴标签
    p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 1))
  }
  
  # 将图存入列表
  plot_list[[i]] <- p
  
  cat(paste0("已生成: ", ko, "\n"))
}

# 设置图标题（可选）
# 如果需要为整个组合图添加标题，可以使用patchwork的plot_annotation
library(patchwork)

# 组合图形（7个图，1行布局）
combined_plot <- wrap_plots(plot_list, ncol = 7) +
  plot_layout(widths = c(1, 1, 1, 1, 1, 1, 1),  # 稍微调整宽度比例
              heights = 1,
              guides = "collect") &
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),  # 减小图形边距
    panel.spacing = unit(2, "mm")  # 减小子图之间的间距
  )

# 设置保存参数
name <- file.path(plot_dir, "representative_KOs_top10orders-wide-version2")
width <- 16  # 增加宽度以适应3列
height <- 6  # 增加高度以适应更多行

# 保存组合图形
ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")

# 显示组合图
print(combined_plot)

