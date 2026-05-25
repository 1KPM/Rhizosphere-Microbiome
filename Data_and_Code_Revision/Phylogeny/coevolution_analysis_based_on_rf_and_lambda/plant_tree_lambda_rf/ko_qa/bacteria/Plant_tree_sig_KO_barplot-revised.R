# 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# 加载必要的R包
library(dplyr)
library(tibble)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(ape)
library(tidyr)

# 1. 加载必要的元数据和颜色信息
plant_colors <- read.csv("../../../../metadata/Tree_top_order_color.csv")

# 2. 加载细菌KO基因的统计结果和定量数据
bac_lambda_rf_ko <- read.csv(file = "./stats/TableS7D-bac_KO_qa_log_lamda_indicators.csv", row.names = 1) 
table(bac_lambda_rf_ko$padj < 0.05)
# 58.86% of core bacterial KOs (7326/12447) exhibited significant phylogenetic conservation (Pagel’s λ, FDR < 0.05) 

# 虽然显著性检验使用的是绝对丰度，但使用TPM数据展示更好看（避免绝对丰度异常值）
bac_ko_qa <- read.csv(file = "../../../../kegg_abundance/relative/Bacteria_core0.2_tpm.csv", row.names = 1, header = T)
names(bac_ko_qa) <- gsub("RS.", "", names(bac_ko_qa))
bac_ko_qa <- rownames_to_column(bac_ko_qa)
names(bac_ko_qa)[1] <- "KOID"

# 3. 筛选出与植物目共进化且显著的KO
bac_order_coevoled_ko <- bac_lambda_rf_ko[!is.na(bac_lambda_rf_ko$Tree_Order) & bac_lambda_rf_ko$padj < 0.05,]
# 2779个KO，还不是最终的筛选

# 4. 进一步筛选：KO对应的植物科属于同一目，且满足多重检验 （和丰度阈值）
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
bac_order_coevoled_top_ko <- merge(bac_order_coevoled_ko, unique(metadata_tree[,c("Family","Order")]),
                                   by.x = "Tree_Family", by.y = "Family", all.x = T)
bac_order_coevoled_top_ko <- bac_order_coevoled_top_ko[bac_order_coevoled_top_ko$Tree_Order == bac_order_coevoled_top_ko$Order,]
names(bac_order_coevoled_top_ko)
bac_order_coevoled_top_ko <- bac_order_coevoled_top_ko[,c(2:15,1,16)]
bac_order_coevoled_top_ko <- subset(bac_order_coevoled_top_ko, Tree_Order != "Others")
# 508个有发育信号的KO对应的植物科属于同一目 (这是很重要排除假阳性的筛选)
# 保存这个筛选结果
bac_lambda_rf_ko$Coevolved <- NA  # 初始化新列
matching <- bac_lambda_rf_ko$KOID %in% bac_order_coevoled_top_ko$KOID
bac_lambda_rf_ko$Coevolved[matching] <- bac_lambda_rf_ko$Tree_Order[matching]
write.csv(bac_lambda_rf_ko, file = "./stats/TableS7D-revised-bac_KO_qa_log_lamda_indicators.csv") 

# 丰度阈值筛选，可选
bac_ko_qa_mean <- data.frame(KOID = bac_ko_qa$KOID, TPM_mean = rowMeans(bac_ko_qa[,-1]))
bac_order_coevoled_top_ko <- merge.data.frame(bac_ko_qa_mean, bac_order_coevoled_top_ko, by = "KOID")

# bac_order_coevoled_top_ko <- subset(bac_order_coevoled_top_ko, padj < 0.01 & rf_FDR_Order < 0.01 & TPM_mean > 1)
# 42个KO
bac_order_coevoled_top_ko <- subset(bac_order_coevoled_top_ko, padj < 0.01 & rf_FDR_Order < 0.01)
# 67个KO

# 5. 准备用于热图绘制的KO丰度数据
bac_order_coevoled_top_ko_qa <- merge.data.frame(bac_order_coevoled_top_ko[,c("KOID", "Tree_Order")], bac_ko_qa, by = "KOID")
common_samples <- intersect(read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')$tip.label,
                            colnames(bac_order_coevoled_top_ko_qa)[!colnames(bac_order_coevoled_top_ko_qa) %in% c("KOID", "Tree_Order")])
ko_tpm_data_filtered <- bac_order_coevoled_top_ko_qa[, c("KOID", "Tree_Order", common_samples)]

# 数据标准化 (按行Z-score)
ko_tpm_data_filtered[,-c(1:2)] <- log(ko_tpm_data_filtered[,-c(1:2)]+1)
ko_tpm_data_filtered[,-c(1:2)] <- t(scale(t(ko_tpm_data_filtered[,-c(1:2)])))

# 转换成长格式
ko_tpm_long_scaled <- ko_tpm_data_filtered %>%
  pivot_longer(
    cols = all_of(common_samples),
    names_to = "PlantSample",
    values_to = "TPM_scaled"
  )

# 设置因子水平，保证与树的顺序一致
ko_order_info <- unique(ko_tpm_long_scaled[, c("KOID", "Tree_Order")])
ko_order_info <- ko_order_info[order(ko_order_info$Tree_Order, ko_order_info$KOID), ]
ko_tpm_long_scaled$KOID <- factor(ko_tpm_long_scaled$KOID, levels = ko_order_info$KOID)
ko_tpm_long_scaled$Tree_Order <- factor(ko_tpm_long_scaled$Tree_Order,
                                        levels = plant_colors[plant_colors$Order %in% names(table(ko_tpm_long_scaled$Tree_Order)),]$Order)

# 6. 准备系统发育树
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
fin_tree <- drop.tip(tree_file, tree_file$tip.label[!tree_file$tip.label %in% common_samples])
ko_tpm_long_scaled$PlantSample <- factor(ko_tpm_long_scaled$PlantSample, levels = fin_tree$tip.label)

metadata_tree <- subset(metadata_tree, metadata_tree$TreeID %in% fin_tree$tip.label)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")
metadata_tree$Order <- factor(metadata_tree$Order, levels = plant_colors$Order)

plant_groupInfo <- split(metadata_tree$TreeID, metadata_tree$Order)
plant_tree_grouped <- groupOTU(fin_tree, plant_groupInfo, group_name = "Lineage")

# 7. 创建分面热图
plant_tree <- ggtree(plant_tree_grouped, aes(color = Lineage), size = 1,
                     ladderize = T, linewidth = 0.2, alpha = 0.6,
                     layout = 'fan', branch.length = "none", open.angle = 15) +
  scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) +
  new_scale_color()

# 按植物目拆分数据并计算每个分面的宽度
ko_order_list <- split(ko_tpm_long_scaled, ko_tpm_long_scaled$Tree_Order)
pwidth_each <- table(ko_tpm_long_scaled$Tree_Order)/dim(ko_tpm_long_scaled)[1]
line.colors <- plant_colors[plant_colors$Order %in% names(table(ko_tpm_long_scaled$Tree_Order)),]$Color2

min_tpm <- min(ko_tpm_long_scaled$TPM_scaled)
max_tpm <- max(ko_tpm_long_scaled$TPM_scaled)
plant_final_plot_scaled <- plant_tree
n_orders <- length(names(ko_order_list))

# 循环添加每个目的热图分面
for(i in 1:n_orders) {
  order_name <- names(ko_order_list)[i]
  order_data <- ko_order_list[[i]]
  
  plant_final_plot_scaled <- plant_final_plot_scaled +
    geom_fruit(
      data = order_data,
      geom = geom_tile,
      aes(y = PlantSample, x = KOID, fill = TPM_scaled),
      color = NA,
      size = 2,
      linewidth = 1,
      offset = 0.03,
      pwidth = pwidth_each[i],
      axis.params = list(
        axis = "x",
        line.color = line.colors[i],
        line.size = 1,
        text.size = 0.5,
        nbreak = 2,
        text.angle = 45,
        vjust = 1,
        hjust = 1
      )
    )
}

# 设置颜色梯度
plant_final_plot_scaled <- plant_final_plot_scaled +
  scale_fill_gradientn(
    colours = c("#596A98", "#596A98", "white", "#E41A1C", "#E41A1C"),
    values = scales::rescale(
      c(min_tpm, min_tpm/2, 0, max_tpm/5, max_tpm),
      from = c(min_tpm, max_tpm)
    ),
    name = "Scaled TPM",
    limits = c(min_tpm, max_tpm)
  ) +
  theme(legend.position = "none")

# 8. 保存最终图形
ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-Facet-Final-used-revised.pdf",
       width = 8, height = 8, units = "cm", dpi = 900, bg = "white")

ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-Facet-Final-used-revised.png",
       width = 8, height = 8, units = "cm", dpi = 900, bg = "white")
