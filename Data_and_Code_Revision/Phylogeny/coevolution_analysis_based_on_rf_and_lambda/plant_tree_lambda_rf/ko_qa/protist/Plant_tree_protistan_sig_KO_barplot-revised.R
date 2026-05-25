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
pro_lambda_rf_ko <- read.csv(file = "./stats/TableS7F-pro_KO_qa_log_lamda_indicators.csv", row.names = 1) 
table(pro_lambda_rf_ko$padj < 0.05)
# functional coevolution was limited in fungi ((160/2458 = 6.51%) of KOs) 

# 虽然显著性检验使用的是绝对丰度，但使用TPM数据展示更好看（避免绝对丰度异常值）
pro_ko_qa <- read.csv(file = "../../../../kegg_abundance/relative/Protists_core0.2_tpm.csv", row.names = 1, header = T)
names(pro_ko_qa) <- gsub("RS.", "", names(pro_ko_qa))
pro_ko_qa <- rownames_to_column(pro_ko_qa)
names(pro_ko_qa)[1] <- "KOID"

# 3. 筛选出与植物目共进化且显著的KO
pro_order_coevoled_ko <- pro_lambda_rf_ko[!is.na(pro_lambda_rf_ko$Tree_Order) & pro_lambda_rf_ko$padj < 0.05,]
# 100个KO，还不是最终的筛选

# 4. 进一步筛选：KO对应的植物科属于同一目，且满足多重检验 （和丰度阈值）
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
pro_order_coevoled_top_ko <- merge(pro_order_coevoled_ko, unique(metadata_tree[,c("Family","Order")]),
                                   by.x = "Tree_Family", by.y = "Family", all.x = T)
pro_order_coevoled_top_ko <- pro_order_coevoled_top_ko[pro_order_coevoled_top_ko$Tree_Order == pro_order_coevoled_top_ko$Order,]
names(pro_order_coevoled_top_ko)
pro_order_coevoled_top_ko <- pro_order_coevoled_top_ko[,c(2:15,1,16)]
pro_order_coevoled_top_ko <- subset(pro_order_coevoled_top_ko, Tree_Order != "Others")

# 只有18个有发育信号的KO对应的植物科属于同一目 (这是很重要排除假阳性的筛选)
# 保存这个筛选结果
pro_lambda_rf_ko$Coevolved <- NA  # 初始化新列
matching <- pro_lambda_rf_ko$KOID %in% pro_order_coevoled_top_ko$KOID
pro_lambda_rf_ko$Coevolved[matching] <- pro_lambda_rf_ko$Tree_Order[matching]
write.csv(pro_lambda_rf_ko, file = "./stats/TableS7F-revised-pro_KO_qa_log_lamda_indicators.csv") 

# 丰度阈值筛选，可选
pro_ko_qa_mean <- data.frame(KOID = pro_ko_qa$KOID, TPM_mean = rowMeans(pro_ko_qa[,-1]))
pro_order_coevoled_top_ko <- merge.data.frame(pro_ko_qa_mean, pro_order_coevoled_top_ko, by = "KOID")

# 继续筛选
# pro_order_coevoled_top_ko1 <- subset(pro_order_coevoled_top_ko, padj < 0.01 & rf_FDR_Order < 0.01)
# 8个KO

# 5. 准备用于热图绘制的KO丰度数据
pro_order_coevoled_top_ko_qa <- merge.data.frame(pro_order_coevoled_top_ko[,c("KOID", "Tree_Order")], pro_ko_qa, by = "KOID")

# write.csv(pro_order_coevoled_top_ko_qa, file = "./stats/pro_order_coevoled_top18_ko_qa.csv")
# write.csv(pro_order_coevoled_top_ko, file = "./stats/pro_order_coevoled_top18_ko_lambda_rf.csv")

common_samples <- intersect(read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')$tip.label,
                            colnames(pro_order_coevoled_top_ko_qa)[!colnames(pro_order_coevoled_top_ko_qa) %in% c("KOID", "Tree_Order")])
ko_tpm_data_filtered <- pro_order_coevoled_top_ko_qa[, c("KOID", "Tree_Order", common_samples)]

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
# pwidth_each <- table(ko_tpm_long_scaled$Tree_Order)/dim(ko_tpm_long_scaled)[1]
pwidth_each <- c(1.5, 0.3, 1.5, 1.5, 0.1)  #pwidth_each这个参数始终没有搞明白
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
        text.size = 1,
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
       filename = "./plant_tree_barplot/plant_order_level_coevolved_Protistan_top_ko-Facet-Final-used-revised.pdf",
       width = 8, height = 8, units = "cm", dpi = 900, bg = "white")

ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_Protistan_top_ko-Facet-Final-used-revised.png",
       width = 8, height = 8, units = "cm", dpi = 900, bg = "white")
