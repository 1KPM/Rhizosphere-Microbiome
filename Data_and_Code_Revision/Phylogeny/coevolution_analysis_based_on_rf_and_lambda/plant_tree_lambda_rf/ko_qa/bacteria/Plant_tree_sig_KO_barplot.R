### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

######################################################################
######################################################################
library(dplyr)
library(tibble)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(reshape2)
library(scales)
library(Rmisc)
library(tidyr)
library(ape)

### 绘制ko水平的系统发育树+barplot图

plant_colors <- read.csv("../../../../metadata/Tree_top_order_color.csv")

bac_lambda_rf_ko <- read.csv(file = "./stats/bac_KO_qa_log_lamda_indicators.csv", row.names = 1)

m <- unique(bac_lambda_rf_ko$KOID) ###12447 KO

n <- unique(subset(bac_lambda_rf_ko, padj < 0.05)$KOID) ###7326 KO


table(!is.na(bac_lambda_rf_ko$Tree_Order) & bac_lambda_rf_ko$padj < 0.05) 
table(subset(bac_lambda_rf_ko, !is.na(Tree_Order) & padj < 0.05)$Tree_Order)

###
bac_ko_ra <- read.csv(file = "../../../../kegg_abundance/relative/Bacteria_core0.2_tpm.csv", row.names = 1, header = T)
names(bac_ko_ra) <- gsub("RS.", "", names(bac_ko_ra)) #693个tree
bac_ko_qa <- read.csv(file = "../../../../kegg_abundance/absolute/Bacteria_core0.2_quantitative.csv", row.names = 1, header = T)
names(bac_ko_qa) <- gsub("RS.", "", names(bac_ko_qa)) #660个tree


####因为有的样本绝对丰度太高，导致所有KO丰度都很高，
####不适合绘制绝对丰度，使用qa样本的相对丰度TPM绘制（coevolved KO还是是基于logqa鉴定的）
bac_ko_qa <- bac_ko_ra[,names(bac_ko_ra) %in% names(bac_ko_qa)]
bac_ko_qa <- rownames_to_column(bac_ko_qa)
names(bac_ko_qa)[1] <- "KOID"

######################################################################
######################################################################

### 在以上的基础上，（1）选择CO-EVOLVED ko绘制植物水平的系统发育树；取其总丰度呈现
bac_order_coevoled_ko <- bac_lambda_rf_ko[!is.na(bac_lambda_rf_ko$Tree_Order) & bac_lambda_rf_ko$padj < 0.05,]
bac_order_coevoled_ko_qa <- merge.data.frame(bac_order_coevoled_ko[,c("KOID", "Tree_Order")], bac_ko_qa, by = "KOID")
bac_order_coevoled_ko_qa <- unique.data.frame(bac_order_coevoled_ko_qa)
bac_order_coevoled_ko_qa_total <- aggregate(bac_order_coevoled_ko_qa[,-c(1:2)], by = list(bac_order_coevoled_ko_qa$Tree_Order), FUN = sum)
bac_order_coevoled_ko_qa_total <- bac_order_coevoled_ko_qa_total %>% column_to_rownames("Group.1")
bac_order_coevoled_ko_qa_total <- as.data.frame(t(bac_order_coevoled_ko_qa_total))
bac_order_coevoled_ko_qa_total <- bac_order_coevoled_ko_qa_total[, plant_colors$Order]
names(bac_order_coevoled_ko_qa_total) <- paste(names(bac_order_coevoled_ko_qa_total), "_coevolved_ko", sep = "")
  
#######################################tree
### 
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% row.names(bac_order_coevoled_ko_qa_total)]
fin_tree <- drop.tip(tree_file, drop_list)
fin_table <- bac_order_coevoled_ko_qa_total[row.names(bac_order_coevoled_ko_qa_total) %in% fin_tree$tip.label,]

###
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
metadata_tree <- subset(metadata_tree, metadata_tree$TreeID %in% fin_tree$tip.label)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")
table(metadata_tree$Order)
metadata_tree$Order <- factor(metadata_tree$Order, levels = plant_colors$Order)

fin_table <- fin_table[fin_tree$tip.label, ] ###重新排序与tip.label一致
fin_table <- fin_table %>% rownames_to_column("label")  ###ggtree的名称是label, 不能改为其它的

### 
plant_groupInfo <- split(
  metadata_tree$TreeID, 
  metadata_tree$Order  # 使用Order作为分组依据
)
# 使用groupOTU整合分组信息
plant_tree_grouped <- groupOTU(
  fin_tree,
  plant_groupInfo,
  group_name = "Lineage"  # 自定义分组变量名
)


### 完整可视化流程 ----------------------------------------------------------

plant_tree <- ggtree(plant_tree_grouped, aes(color=Lineage), size = 1, ladderize = T, linewidth = 0.1, 
                      layout = 'fan', branch.length = "none", right = TRUE, open.angle = 180) +
  scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) + 
  new_scale_color()

# 定义需要绘制的列名
coevoled_ko_columns <- names(fin_table)[-1]
# 创建命名颜色向量（确保顺序匹配）
color_vector <- setNames(plant_colors$Color2, coevoled_ko_columns)

###
plant_final_plot1 <- plant_tree
for (col in coevoled_ko_columns) {
  # 为当前列获取对应的颜色
  col_color <- color_vector[col]
  
  plant_final_plot1 <- plant_final_plot1 + 
    geom_fruit(
      data = fin_table,
      geom = geom_bar,
      stat = 'identity',
      width = 2,
      aes_string(y = "label", x = col),  # 使用 aes_string 简化语法
      pwidth = 0.1,
      fill = col_color
    )
}
plant_final_plot1

# 保存结果
ggsave(plant_final_plot1,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_ko_loop-TPM.pdf", 
       width = 12, height = 12, units = "cm", dpi = 900, bg = "white")

###
final_plot1 <- plant_final_plot1 + theme(legend.position = "none")
ggsave(final_plot1,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_ko_loop_nolegend-TPM.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot1,
       filename = "plant_tree_barplot/plant_order_level_coevolved_bacterial_ko_loop_nolegend-TPM.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot1,
       filename = "plant_tree_barplot/plant_order_level_coevolved_bacterial_ko_loop_nolegend-TPM.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")



#####(2) 可视化特别的一些pathway
######################################################################
######################################################################
### 在以上的基础上，（1）选择CO-EVOLVED ko绘制植物水平的系统发育树；取其总丰度呈现
bac_order_coevoled_ko <- bac_lambda_rf_ko[!is.na(bac_lambda_rf_ko$Tree_Order) & bac_lambda_rf_ko$padj < 0.05,]
select_pathways <- c("Sporulation and spore development",
                     "Symbiotic nitrogen fixation",
                     "Secretion systems",
                     "Quorum sensing",
                     "Methane metabolism")

# 使用正则表达式模糊匹配包含"spore"或"sporulation"的描述
Sporulation_ko <- bac_order_coevoled_ko[
  grepl("spore|sporulation", bac_order_coevoled_ko$Description, ignore.case = TRUE),
]
Sporulation_ko$Pathway <- "Sporulation_and_spore_development"

Symbiotic_ko <- bac_order_coevoled_ko[
  grepl("nod|nif|vnf|anf", bac_order_coevoled_ko$gene, ignore.case = FALSE),
]
Symbiotic_ko$Pathway <- "Symbiotic_nitrogen_fixation"

Secretion_ko <- bac_order_coevoled_ko[
  grepl("type III secretion|type IV secretion", bac_order_coevoled_ko$Description, ignore.case = TRUE),
]
Secretion_ko$Pathway <- "Secretion_systems"

Quorum_ko <- bac_order_coevoled_ko[
  grepl("quorum-sensing system", bac_order_coevoled_ko$Description, ignore.case = TRUE),
]
Quorum_ko$Pathway <- "Quorum_sensing_system"

Methane_ko <- bac_order_coevoled_ko[
  grepl("methane/ammonia monooxygenase", bac_order_coevoled_ko$Description, ignore.case = FALSE),
]
Methane_ko$Pathway <- "Methane_metabolism"

pathways_list <- rbind(Sporulation_ko, Symbiotic_ko, Secretion_ko, Quorum_ko, Methane_ko)
pathways_list <- unique(pathways_list[,c("KOID", "Pathway")])
write.csv(pathways_list, file = "./stats/pathways_list.csv")



bac_order_coevoled_module_qa <- merge.data.frame(pathways_list, bac_ko_qa, by = "KOID")
bac_order_coevoled_module_qa_total <- aggregate(bac_order_coevoled_module_qa[,-c(1:2)], by = list(bac_order_coevoled_module_qa$Pathway), FUN = sum)
bac_order_coevoled_module_qa_total <- bac_order_coevoled_module_qa_total %>% column_to_rownames("Group.1")
bac_order_coevoled_module_qa_total <- as.data.frame(t(bac_order_coevoled_module_qa_total))

#######################################tree
### 
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% row.names(bac_order_coevoled_module_qa_total)]
fin_tree <- drop.tip(tree_file, drop_list)
fin_table_module <- bac_order_coevoled_module_qa_total[row.names(bac_order_coevoled_module_qa_total) %in% fin_tree$tip.label,]

###
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
metadata_tree <- subset(metadata_tree, metadata_tree$TreeID %in% fin_tree$tip.label)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")
table(metadata_tree$Order)
metadata_tree$Order <- factor(metadata_tree$Order, levels = plant_colors$Order)

fin_table_module <- fin_table_module[fin_tree$tip.label, ] ###重新排序与tip.label一致
fin_table_module <- fin_table_module %>% rownames_to_column("label")  ###ggtree的名称是label, 不能改为其它的

### 
plant_groupInfo <- split(
  metadata_tree$TreeID, 
  metadata_tree$Order  # 使用Order作为分组依据
)
# 使用groupOTU整合分组信息
plant_tree_grouped <- groupOTU(
  fin_tree,
  plant_groupInfo,
  group_name = "Lineage"  # 自定义分组变量名
)


### 完整可视化流程 ----------------------------------------------------------

plant_tree <- ggtree(plant_tree_grouped, aes(color=Lineage), size = 1, ladderize = T, linewidth = 0.1, 
                     layout = 'fan', branch.length = "none", right = TRUE, open.angle = 180) +
  scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) + 
  new_scale_color()

# 定义需要绘制的列名
coevoled_module_columns <- names(fin_table_module)[-1][c(4,5,1,3,2)]
# 创建命名颜色向量（确保顺序匹配）
color_vector <- setNames(plant_colors$Color2[1:5], coevoled_module_columns)

###
plant_final_plot2 <- plant_tree
for (col in coevoled_module_columns) {
  # 为当前列获取对应的颜色
  col_color <- color_vector[col]
  
  plant_final_plot2 <- plant_final_plot2 + 
    geom_fruit(
      data = fin_table_module,
      geom = geom_bar,
      stat = 'identity',
      width = 2,
      aes_string(y = "label", x = col),  # 使用 aes_string 简化语法
      pwidth = 0.2,
      fill = col_color
    )
}
plant_final_plot2

# 保存结果
ggsave(plant_final_plot2,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_module_loop-TPM.pdf", 
       width = 12, height = 12,  dpi = 900, bg = "white")

###
final_plot2 <- plant_final_plot2 + theme(legend.position = "none")
ggsave(final_plot2,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_module_loop_nolegend-TPM.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot2,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_module_loop_nolegend-TPM.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot2,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_module_loop_nolegend-TPM.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")


###################################################################################
########### 可视化每个目最显著的 KO，要求：
###########（1）该KO对应的植物科也属于这个植物目
###########(2) padj < 0.01 & rf_FDR_Order < 0.01
###########（3）TPM > 1, 不是随机造成的

bac_order_coevoled_top_ko <- merge(bac_order_coevoled_ko, unique(metadata_tree[,c("Family","Order")]), 
                                     by.x = "Tree_Family", by.y = "Family", all.x = T)


bac_order_coevoled_top_ko <- bac_order_coevoled_top_ko[bac_order_coevoled_top_ko$Tree_Order == bac_order_coevoled_top_ko$Order,]
names(bac_order_coevoled_top_ko)
bac_order_coevoled_top_ko <- unique(bac_order_coevoled_top_ko[,c(2:11,18:21,1,22)])
bac_order_coevoled_top_ko <- subset(bac_order_coevoled_top_ko, Tree_Order != "Others")

###进一步选择lambda比较高且丰度也比较高的
bac_ko_qa_mean <- data.frame(KOID = bac_ko_qa$KOID, TPM_mean = rowMeans(bac_ko_qa[,-1]))
bac_order_coevoled_top_ko <- merge.data.frame(bac_ko_qa_mean, bac_order_coevoled_top_ko, by = "KOID")
bac_order_coevoled_top_ko <- subset(bac_order_coevoled_top_ko, padj < 0.01 & rf_FDR_Order < 0.01 & TPM_mean > 1) #7个目42个KO

# 系统发育信号具体的强弱划分（经验法则）
# λ > 0.8：通常被认为存在很强的系统发育信号。这表明性状在亲缘关系较近的物种间高度相似，系统发育关系在解释性状变异时至关重要。在分析时必须考虑系统发育结构（使用系统发育比较方法）。
# 0.5 < λ < 0.8：存在中等至较强的系统发育信号。系统发育关系有重要影响，但可能也存在其他因素（如适应性趋同）的作用。分析时强烈建议考虑系统发育。
# 0.2 < λ < 0.5：存在较弱的系统发育信号。性状有一定程度的系统发育依赖性，但影响不大。
# λ < 0.2：系统发育信号非常弱或可忽略不计。基本可以认为物种性状是独立演化的。
# λ ≈ 1：如前所述，这是理论上的最强信号。

###（1）可视化7个目42个KO的系统发育树+barplot
###（2）可视化7个目42个KO的相对丰度TPM热图

bac_order_coevoled_top_ko_qa <- merge.data.frame(bac_order_coevoled_top_ko[,c("KOID", "Tree_Order")], bac_ko_qa, by = "KOID")
# 
# write.csv(bac_order_coevoled_top_ko_qa, file = "./stats/bac_order_coevoled_top42_ko_qa.csv")
# write.csv(bac_order_coevoled_top_ko, file = "./stats/bac_order_coevoled_top42_ko_lambda_rf.csv")
# 
bac_order_coevoled_top_ko_qa_total <- aggregate(bac_order_coevoled_top_ko_qa[,-c(1:2)], by = list(bac_order_coevoled_top_ko_qa$Tree_Order), FUN = sum)
bac_order_coevoled_top_ko_qa_total <- bac_order_coevoled_top_ko_qa_total %>% column_to_rownames("Group.1")
bac_order_coevoled_top_ko_qa_total <- as.data.frame(t(bac_order_coevoled_top_ko_qa_total))

#######################################tree
### 
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% row.names(bac_order_coevoled_top_ko_qa_total)]
fin_tree <- drop.tip(tree_file, drop_list)
fin_table_top_ko <- bac_order_coevoled_top_ko_qa_total[row.names(bac_order_coevoled_top_ko_qa_total) %in% fin_tree$tip.label,]

###
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
metadata_tree <- subset(metadata_tree, metadata_tree$TreeID %in% fin_tree$tip.label)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")
table(metadata_tree$Order)
metadata_tree$Order <- factor(metadata_tree$Order, levels = plant_colors$Order)

fin_table_top_ko <- fin_table_top_ko[fin_tree$tip.label, ] ###重新排序与tip.label一致
fin_table_top_ko <- fin_table_top_ko %>% rownames_to_column("label")  ###ggtree的名称是label, 不能改为其它的

### 
plant_groupInfo <- split(
  metadata_tree$TreeID, 
  metadata_tree$Order  # 使用Order作为分组依据
)
# 使用groupOTU整合分组信息
plant_tree_grouped <- groupOTU(
  fin_tree,
  plant_groupInfo,
  group_name = "Lineage"  # 自定义分组变量名
)


### 完整可视化流程 ----------------------------------------------------------

plant_tree <- ggtree(plant_tree_grouped, aes(color=Lineage), size = 1, ladderize = T, linewidth = 0.1, 
                     layout = 'fan', branch.length = "none", right = TRUE, open.angle = 15) +
  scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) + 
  new_scale_color()

# 定义需要绘制的列名
match(plant_colors$Order, names(fin_table_top_ko))
coevoled_top_ko_columns <- names(fin_table_top_ko)[rev(c(4,8,6,5,3,7,2))]
# 创建命名颜色向量（确保顺序匹配）
color_vector <- setNames(plant_colors$Color2[rev(c(1,2,3,6,7,9,10))], coevoled_top_ko_columns)

###
plant_final_plot3 <- plant_tree
for (col in coevoled_top_ko_columns) {
  # 为当前列获取对应的颜色
  col_color <- color_vector[col]
  
  plant_final_plot3 <- plant_final_plot3 + 
    geom_fruit(
      data = fin_table_top_ko,
      geom = geom_bar,
      stat = 'identity',
      width = 2,
      aes_string(y = "label", x = col),  # 使用 aes_string 简化语法
      pwidth = 0.2,
      fill = col_color
    )
}
plant_final_plot3

# 保存结果
ggsave(plant_final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko_loop-TPM.pdf", 
       width = 12, height = 12,  dpi = 900, bg = "white")

###
final_plot3 <- plant_final_plot3 + theme(legend.position = "none")
ggsave(final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko_loop_nolegend-TPM.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko_loop_nolegend-TPM.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_module_loop_nolegend-TPM.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")


######################################################################
###（2）可视化7个目42个KO的相对丰度TPM热图
# 读取植物系统发育树
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')

# 读取植物元数据
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
metadata_tree <- subset(metadata_tree, metadata_tree$TreeID %in% tree_file$tip.label)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")
metadata_tree$Order <- factor(metadata_tree$Order, levels = plant_colors$Order)

# 准备热图数据
# bac_order_coevoled_top_ko_qa 应该包含：KOID, Tree_Order, 以及每个植物样本的TPM值
# 我们需要将其转换为长格式，行为植物样本，列为KO

# 确保树中的植物样本在数据中存在
plant_samples_in_tree <- tree_file$tip.label
sample_cols <- colnames(bac_order_coevoled_top_ko_qa)[!colnames(bac_order_coevoled_top_ko_qa) %in% c("KOID", "Tree_Order")]

# 找出树和数据中共同存在的样本
common_samples <- intersect(plant_samples_in_tree, sample_cols)
cat("树中有", length(plant_samples_in_tree), "个植物样本\n")
cat("数据中有", length(sample_cols), "个植物样本\n")
cat("共同存在的样本有", length(common_samples), "个\n")

# 过滤树，只保留共同样本
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% common_samples]
fin_tree <- drop.tip(tree_file, drop_list)

# 过滤TPM数据，只保留共同样本
ko_tpm_data_filtered <- bac_order_coevoled_top_ko_qa[, c("KOID", "Tree_Order", common_samples)]
ko_tpm_data_filtered[,-c(1:2)] <- log(ko_tpm_data_filtered[,-c(1:2)]+1)
ko_tpm_data_filtered[,-c(1:2)] <- t(scale(t(ko_tpm_data_filtered[,-c(1:2)])))   #按行（样本）标准化

# 准备热图数据（宽格式 -> 长格式）
ko_tpm_long_scaled <- ko_tpm_data_filtered %>%
  pivot_longer(
    cols = all_of(common_samples),
    names_to = "PlantSample",
    values_to = "TPM_scaled"
  )

# 确保植物样本顺序与树一致
ko_tpm_long_scaled$PlantSample <- factor(ko_tpm_long_scaled$PlantSample, levels = fin_tree$tip.label)

# 确保KO按Order排序
ko_order_info <- unique(ko_tpm_long_scaled[, c("KOID", "Tree_Order")])
ko_order_info <- ko_order_info[order(ko_order_info$Tree_Order, ko_order_info$KOID), ]
ko_tpm_long_scaled$KOID <- factor(ko_tpm_long_scaled$KOID, levels = ko_order_info$KOID)

# 为树添加分组信息
plant_groupInfo <- split(metadata_tree$TreeID, metadata_tree$Order)
plant_tree_grouped <- groupOTU(fin_tree, plant_groupInfo, group_name = "Lineage")

# 绘制系统发育树（矩形布局）
plant_tree <- ggtree(plant_tree_grouped, aes(color = Lineage), size = 1, ladderize = T, linewidth = 0.2, alpha = 0.6,
                     layout = 'fan', branch.length = "none", open.angle = 15) +
  scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) + 
  new_scale_color()

###
min <- min(ko_tpm_long_scaled$TPM_scaled)
max <- max(ko_tpm_long_scaled$TPM_scaled)

plant_final_plot_scaled <- plant_tree +
  geom_fruit(
    data = ko_tpm_long_scaled,
    geom = geom_tile,
    aes(y = PlantSample, x = KOID, fill = TPM_scaled),
    color = NA,
    size = 2,
    linewidth = 1,
    offset = 0.01,
    pwidth = 1,
    #给热图标注行名
    axis.params = list(
      axis="x",
      line.color="grey",
      text.size = 1,
      nbreak = 2, 
      text.angle = 45, 
      vjust = 1, 
      hjust = 1
    )
  ) +
  scale_fill_gradientn(
    colours = c("#596A98", "#596A98", "white", "#E41A1C", "#E41A1C"),
    values = scales::rescale(
      c(min, min/2, 0, max/2, max),  # 这些位置对应的颜色
      from = c(min, max)  # 原始数据范围
    ),
    name = "Scaled TPM",
    limits = c(min, max)  # 确保颜色范围与数据范围一致
  )

plant_colors$Color2
# "#E41A1C" "#596A98" "#449B75" "#6B886D" "#AC5782" "#FF7F00" "#FFE528" "#C9992C" "#C66764" "#E485B7" "#999999"
# 保存结果
ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-TPM.pdf", 
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


plant_final_plot4 <- plant_final_plot_scaled + theme(legend.position = "none")

ggsave(plant_final_plot4,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-TPM-nolegend.pdf", 
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


ggsave(plant_final_plot4,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko_nolegend-TPM.png",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")

ggsave(plant_final_plot4,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko_nolegend-TPM.jpg",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


###############

# 将数据按Tree_Order分组
ko_tpm_long_scaled$Tree_Order <- factor(ko_tpm_long_scaled$Tree_Order, 
                                        levels = plant_colors[plant_colors$Order %in% names(table(ko_tpm_long_scaled$Tree_Order)),]$Order)
ko_order_list <- split(ko_tpm_long_scaled, ko_tpm_long_scaled$Tree_Order)

# 计算每个分组的偏移量和宽度
pwidth_each <- table(ko_tpm_long_scaled$Tree_Order)/dim(ko_tpm_long_scaled)[1]
line.colors <- plant_colors[plant_colors$Order %in% names(table(ko_tpm_long_scaled$Tree_Order)),]$Color2

# 创建基础树图
plant_tree <- ggtree(plant_tree_grouped, aes(color = Lineage), size = 1, 
                     ladderize = T, linewidth = 0.2, alpha = 0.6,
                     layout = 'fan', branch.length = "none", open.angle = 15) +
  scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) + 
  new_scale_color()

# 循环添加每个Order的热图
plant_final_plot_scaled <- plant_tree
n_orders <- length(names(ko_order_list))

for(i in 1:n_orders) {
  order_name <- names(ko_order_list)[i]
  order_data <- ko_order_list[[i]]
  
  # 为每个分组创建热图
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

# 添加颜色标度
plant_final_plot_scaled <- plant_final_plot_scaled +
  scale_fill_gradientn(
    colours = c("#596A98", "#596A98", "white", "#E41A1C", "#E41A1C"),
    values = scales::rescale(
      c(min, min/2, 0, max/2, max),
      from = c(min, max)
    ),
    name = "Scaled TPM",
    limits = c(min, max)
  )

plant_final_plot_scaled <- plant_final_plot_scaled + theme(legend.position = "none")

ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-Facet-Final-used.pdf", 
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-Facet-Final-used.png",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")

ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_bacterial_top_ko-Facet-Final-used.jpg",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")











