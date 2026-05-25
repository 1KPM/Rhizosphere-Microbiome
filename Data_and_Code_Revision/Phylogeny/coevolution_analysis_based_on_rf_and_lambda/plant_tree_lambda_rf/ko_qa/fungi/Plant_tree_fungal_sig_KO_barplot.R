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

fun_lambda_rf_ko <- read.csv(file = "./stats/fun_KO_qa_log_lamda_indicators.csv", row.names = 1)

m <- unique(fun_lambda_rf_ko$KOID) ###5160 KO

n <- unique(subset(fun_lambda_rf_ko, padj < 0.05)$KOID) ###261 KO


table(!is.na(fun_lambda_rf_ko$Tree_Order) & fun_lambda_rf_ko$padj < 0.05) 
table(subset(fun_lambda_rf_ko, !is.na(Tree_Order) & padj < 0.05)$Tree_Order)

###
fun_ko_ra <- read.csv(file = "../../../../kegg_abundance/relative/Fungi_core0.2_tpm.csv", row.names = 1, header = T)
names(fun_ko_ra) <- gsub("RS.", "", names(fun_ko_ra)) #693个tree
fun_ko_qa <- read.csv(file = "../../../../kegg_abundance/absolute/Fungi_core0.2_quantitative.csv", row.names = 1, header = T)
names(fun_ko_qa) <- gsub("RS.", "", names(fun_ko_qa)) #660个tree


####因为有的样本绝对丰度太高，导致所有KO丰度都很高，
####不适合绘制绝对丰度，使用qa样本的相对丰度TPM绘制（coevolved KO还是是基于logqa鉴定的）
fun_ko_qa <- fun_ko_ra[,names(fun_ko_ra) %in% names(fun_ko_qa)]
fun_ko_qa <- rownames_to_column(fun_ko_qa)
names(fun_ko_qa)[1] <- "KOID"

######################################################################
######################################################################

### 在以上的基础上，（1）选择CO-EVOLVED ko绘制植物水平的系统发育树；取其总丰度呈现
fun_order_coevoled_ko <- fun_lambda_rf_ko[!is.na(fun_lambda_rf_ko$Tree_Order) & fun_lambda_rf_ko$padj < 0.05,]
fun_order_coevoled_ko_qa <- merge.data.frame(fun_order_coevoled_ko[,c("KOID", "Tree_Order")], fun_ko_qa, by = "KOID")
fun_order_coevoled_ko_qa <- unique.data.frame(fun_order_coevoled_ko_qa)
fun_order_coevoled_ko_qa_total <- aggregate(fun_order_coevoled_ko_qa[,-c(1:2)], by = list(fun_order_coevoled_ko_qa$Tree_Order), FUN = sum)
fun_order_coevoled_ko_qa_total <- fun_order_coevoled_ko_qa_total %>% column_to_rownames("Group.1")
fun_order_coevoled_ko_qa_total <- as.data.frame(t(fun_order_coevoled_ko_qa_total))

fun_co_orders <- intersect(plant_colors$Order, names(fun_order_coevoled_ko_qa_total)) 
fun_order_coevoled_ko_qa_total <- fun_order_coevoled_ko_qa_total[, fun_co_orders]
names(fun_order_coevoled_ko_qa_total) <- paste(names(fun_order_coevoled_ko_qa_total), "_coevolved_ko", sep = "")
  
#######################################tree
### 
tree_file <- read.tree('../../../../metadata/tree_metadata_merge_info_final_align_tree.nwk')
drop_list <- tree_file$tip.label[!tree_file$tip.label %in% row.names(fun_order_coevoled_ko_qa_total)]
fin_tree <- drop.tip(tree_file, drop_list)
fin_table <- fun_order_coevoled_ko_qa_total[row.names(fun_order_coevoled_ko_qa_total) %in% fin_tree$tip.label,]

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
fun_co_colors <- plant_colors[plant_colors$Order %in% fun_co_orders,]$Color2
color_vector <- setNames(fun_co_colors, coevoled_ko_columns)

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
       filename = "./plant_tree_barplot/plant_order_level_coevolved_fungal_ko_loop-TPM.pdf", 
       width = 12, height = 12, units = "cm", dpi = 900, bg = "white")

###
final_plot1 <- plant_final_plot1 + theme(legend.position = "none")
ggsave(final_plot1,
       filename = "./plant_tree_barplot/plant_order_level_coevolved_fungal_ko_loop_nolegend-TPM.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot1,
       filename = "plant_tree_barplot/plant_order_level_coevolved_fungal_ko_loop_nolegend-TPM.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot1,
       filename = "plant_tree_barplot/plant_order_level_coevolved_fungal_ko_loop_nolegend-TPM.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")



#####(2) 可视化每个目最显著的 KO，要求：
###########（1）该KO对应的植物科也属于这个植物目
###########(2) padj < 0.01 & rf_FDR_Order < 0.01
###########（3）TPM > 1, 不是随机造成的

fun_order_coevoled_top_ko <- merge(fun_order_coevoled_ko, unique(metadata_tree[,c("Family","Order")]), 
                                     by.x = "Tree_Family", by.y = "Family", all.x = T)


fun_order_coevoled_top_ko <- fun_order_coevoled_top_ko[fun_order_coevoled_top_ko$Tree_Order == fun_order_coevoled_top_ko$Order,]
names(fun_order_coevoled_top_ko)
fun_order_coevoled_top_ko <- unique(fun_order_coevoled_top_ko[,c(2:11,18:21,1,22)])
fun_order_coevoled_top_ko <- subset(fun_order_coevoled_top_ko, Tree_Order != "Others")

###进一步选择lambda比较高且丰度也比较高的
fun_ko_qa_mean <- data.frame(KOID = fun_ko_qa$KOID, TPM_mean = rowMeans(fun_ko_qa[,-1]))
fun_order_coevoled_top_ko <- merge.data.frame(fun_ko_qa_mean, fun_order_coevoled_top_ko, by = "KOID") #4个KO

# fun_order_coevoled_top_ko <- subset(fun_order_coevoled_top_ko, padj < 0.01 & rf_FDR_Order < 0.01 & TPM_mean > 1) 

###过滤没了

#############################################################################################
###以下代码暂时就没有跑了
#############################################################################################
#############################################################################################
###以下代码暂时就没有跑了
#############################################################################################


