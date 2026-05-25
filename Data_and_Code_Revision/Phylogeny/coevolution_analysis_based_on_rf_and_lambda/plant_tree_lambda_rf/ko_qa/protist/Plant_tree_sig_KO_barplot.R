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

pro_lambda_rf_ko <- read.csv(file = "./stats/pro_KO_qa_log_lamda_indicators.csv", row.names = 1)

m <- unique(pro_lambda_rf_ko$KOID) ###2458 KO

n <- unique(subset(pro_lambda_rf_ko, padj < 0.05)$KOID) ###160 KO


table(!is.na(pro_lambda_rf_ko$Tree_Order) & pro_lambda_rf_ko$padj < 0.05) 
table(subset(pro_lambda_rf_ko, !is.na(Tree_Order) & padj < 0.05)$Tree_Order)

###
pro_ko_ra <- read.csv(file = "../../../../kegg_abundance/relative/Protists_core0.2_tpm.csv", row.names = 1, header = T)
names(pro_ko_ra) <- gsub("RS.", "", names(pro_ko_ra)) #693个tree
pro_ko_qa <- read.csv(file = "../../../../kegg_abundance/absolute/Protists_core0.2_quantitative.csv", row.names = 1, header = T)
names(pro_ko_qa) <- gsub("RS.", "", names(pro_ko_qa)) #660个tree


####因为有的样本绝对丰度太高，导致所有KO丰度都很高，
####不适合绘制绝对丰度，使用qa样本的相对丰度TPM绘制（coevolved KO还是是基于logqa鉴定的）
pro_ko_qa <- pro_ko_ra[,names(pro_ko_ra) %in% names(pro_ko_qa)]
pro_ko_qa <- rownames_to_column(pro_ko_qa)
names(pro_ko_qa)[1] <- "KOID"

######################################################################
######################################################################

### 在以上的基础上，（1）选择CO-EVOLVED ko绘制植物水平的系统发育树；取其总丰度呈现
pro_order_coevoled_ko <- pro_lambda_rf_ko[!is.na(pro_lambda_rf_ko$Tree_Order) & pro_lambda_rf_ko$padj < 0.05,]
pro_order_coevoled_ko_qa <- merge.data.frame(pro_order_coevoled_ko[,c("KOID", "Tree_Order")], pro_ko_qa, by = "KOID")
pro_order_coevoled_ko_qa <- unique.data.frame(pro_order_coevoled_ko_qa)
pro_order_coevoled_ko_qa_total <- aggregate(pro_order_coevoled_ko_qa[,-c(1:2)], by = list(pro_order_coevoled_ko_qa$Tree_Order), FUN = sum)
pro_order_coevoled_ko_qa_total <- pro_order_coevoled_ko_qa_total %>% column_to_rownames("Group.1")
pro_order_coevoled_ko_qa_total <- as.data.frame(t(pro_order_coevoled_ko_qa_total))

pro_co_orders <- intersect(plant_colors$Order, names(pro_order_coevoled_ko_qa_total)) 
pro_order_coevoled_ko_qa_total <- pro_order_coevoled_ko_qa_total[, pro_co_orders]
names(pro_order_coevoled_ko_qa_total) <- paste(names(pro_order_coevoled_ko_qa_total), "_coevolved_ko", sep = "")
  
#######################################tree
### 
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% row.names(pro_order_coevoled_ko_qa_total)]
fin_tree <- drop.tip(tree_file, drop_list)
fin_table <- pro_order_coevoled_ko_qa_total[row.names(pro_order_coevoled_ko_qa_total) %in% fin_tree$tip.label,]

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
pro_co_colors <- plant_colors[plant_colors$Order %in% pro_co_orders,]$Color2
color_vector <- setNames(pro_co_colors, coevoled_ko_columns)

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
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_ko_loop-TPM.pdf", 
       width = 12, height = 12, units = "cm", dpi = 900, bg = "white")

###
final_plot1 <- plant_final_plot1 + theme(legend.position = "none")
ggsave(final_plot1,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_ko_loop_nolegend-TPM.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot1,
       filename = "plant_tree_barplot/plant_order_level_coevolved_protistan_ko_loop_nolegend-TPM.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot1,
       filename = "plant_tree_barplot/plant_order_level_coevolved_protistan_ko_loop_nolegend-TPM.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")


##############################################################################################
###以下结果不用了
############################################################################################


#####(2) 可视化每个目最显著的 KO，要求：
###########（1）该KO对应的植物科也属于这个植物目
###########(2) padj < 0.01 & rf_FDR_Order < 0.01
###########（3）TPM > 1, 不是随机造成的

pro_order_coevoled_top_ko <- merge(pro_order_coevoled_ko, unique(metadata_tree[,c("Family","Order")]), 
                                     by.x = "Tree_Family", by.y = "Family", all.x = T)

###也不要求
# pro_order_coevoled_top_ko <- pro_order_coevoled_top_ko[pro_order_coevoled_top_ko$Tree_Order == pro_order_coevoled_top_ko$Order,]
names(pro_order_coevoled_top_ko)
pro_order_coevoled_top_ko <- unique(pro_order_coevoled_top_ko[,c(2:11,18:21,1,22)])
pro_order_coevoled_top_ko <- subset(pro_order_coevoled_top_ko, Tree_Order != "Others")

###进一步选择lambda比较高且丰度也比较高的
pro_ko_qa_mean <- data.frame(KOID = pro_ko_qa$KOID, TPM_mean = rowMeans(pro_ko_qa[,-1]))
pro_order_coevoled_top_ko <- merge.data.frame(pro_ko_qa_mean, pro_order_coevoled_top_ko, by = "KOID")

##不p值过滤还有18个KO
# pro_order_coevoled_top_ko <- subset(pro_order_coevoled_top_ko, padj < 0.01 & rf_FDR_Order < 0.01 & TPM_mean > 1) 
###过滤没了

#####无法满足以上条件，我们直接可视化lambda最高的前42个KO
pro_order_coevoled_top_ko <- pro_order_coevoled_top_ko[order(pro_order_coevoled_top_ko$lambda, decreasing = T), ] 
pro_order_coevoled_top_ko <- pro_order_coevoled_top_ko[1:42, ]

#####
pro_order_coevoled_top_ko_qa <- merge.data.frame(pro_order_coevoled_top_ko[,c("KOID", "Tree_Order")], pro_ko_qa, by = "KOID")
pro_order_coevoled_top_ko_qa_total <- aggregate(pro_order_coevoled_top_ko_qa[,-c(1:2)], by = list(pro_order_coevoled_top_ko_qa$Tree_Order), FUN = sum)
pro_order_coevoled_top_ko_qa_total <- pro_order_coevoled_top_ko_qa_total %>% column_to_rownames("Group.1")
pro_order_coevoled_top_ko_qa_total <- as.data.frame(t(pro_order_coevoled_top_ko_qa_total))

#######################################tree
### 
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% row.names(pro_order_coevoled_top_ko_qa_total)]
fin_tree <- drop.tip(tree_file, drop_list)
fin_table_top_ko <- pro_order_coevoled_top_ko_qa_total[row.names(pro_order_coevoled_top_ko_qa_total) %in% fin_tree$tip.label,]

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
coevoled_top_ko_columns <- names(fin_table_top_ko)[rev(c(6,4,7,3,2,5))]
# 创建命名颜色向量（确保顺序匹配）
color_vector <- setNames(plant_colors$Color2[rev(c(2,3,5,6,7,8))], coevoled_top_ko_columns)

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
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko_loop-TPM.pdf", 
       width = 12, height = 12,  dpi = 900, bg = "white")

###
final_plot3 <- plant_final_plot3 + theme(legend.position = "none")
ggsave(final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko_loop_nolegend-TPM.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko_loop_nolegend-TPM.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot3,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko_loop_nolegend-TPM.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")


######################################################################
###（2）可视化这些top
# 读取植物系统发育树
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')

# 读取植物元数据
metadata_tree <- read.csv('../../../../metadata/tree_metadata_merge_info.csv')
metadata_tree <- subset(metadata_tree, metadata_tree$TreeID %in% tree_file$tip.label)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")
metadata_tree$Order <- factor(metadata_tree$Order, levels = plant_colors$Order)

# 准备热图数据
# pro_order_coevoled_top_ko_qa 应该包含：KOID, Tree_Order, 以及每个植物样本的TPM值
# 我们需要将其转换为长格式，行为植物样本，列为KO

# 确保树中的植物样本在数据中存在
plant_samples_in_tree <- tree_file$tip.label
sample_cols <- colnames(pro_order_coevoled_top_ko_qa)[!colnames(pro_order_coevoled_top_ko_qa) %in% c("KOID", "Tree_Order")]

# 找出树和数据中共同存在的样本
common_samples <- intersect(plant_samples_in_tree, sample_cols)
cat("树中有", length(plant_samples_in_tree), "个植物样本\n")
cat("数据中有", length(sample_cols), "个植物样本\n")
cat("共同存在的样本有", length(common_samples), "个\n")

# 过滤树，只保留共同样本
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% common_samples]
fin_tree <- drop.tip(tree_file, drop_list)

# 过滤TPM数据，只保留共同样本
ko_tpm_data_filtered <- pro_order_coevoled_top_ko_qa[, c("KOID", "Tree_Order", common_samples)]
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
    offset = 0.03,
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
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko-TPM.pdf", 
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


plant_final_plot4 <- plant_final_plot_scaled + theme(legend.position = "none")

ggsave(plant_final_plot4,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko-TPM-nolegend.pdf", 
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


ggsave(plant_final_plot4,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko_nolegend-TPM.png",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")

ggsave(plant_final_plot4,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko_nolegend-TPM.jpg",
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

order_data <- ko_order_list[[1]]

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
      pwidth = pwidth_each[i],    #pwidth_each这个参数没有调整好，第一个pwidth应该是与进化树直接的相对距离，
      # 第二个pwidth应该是与进化树+第一个热图的相对距离。反复调整都不行，是不是只有一行没法绘制geom_tile?
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
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko-Facet-Final-used.pdf", 
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")


ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko-Facet-Final-used.png",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")

ggsave(plant_final_plot_scaled,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_protistan_top_ko-Facet-Final-used.jpg",
       width = 7, height = 7, units = "cm", dpi = 900, bg = "white")

