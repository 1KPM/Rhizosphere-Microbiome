### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

### 加载必要包 ----------------------------------------------------------------
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(Vennerable)
library(ggplot2)
library(scales)
library(tibble)

### 数据准备与合并 (为后续统计做准备) ------------------------------------------
# 导入核心ASV分类信息
core_taxonomy <- read.csv('../../../feature_table/All_core_ASV_taxonomy.csv', row.names = 1)

# 导入系统发育信号分析结果
lamda_file <- read.csv("../../../12.phylogenetic_signal_λ/amplicon/core_asv/Fungi_absolute_format_core_ASV_phylogenetic_signal.csv", row.names = 1) %>%
  rownames_to_column("FeatureID")

# 导入随机森林（Class, Order, Family水平）鉴定的差异ASV结果
rf_class <- read.csv("../../../13.random_forest/amplicon/fun_qa_log_Class_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
rf_order <- read.csv("../../../13.random_forest/amplicon/fun_qa_log_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
rf_family <- read.csv("../../../13.random_forest/amplicon/fun_qa_log_Family_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)

# 统一处理ASVID
row.names(rf_class) <- gsub("\\|.*", "", row.names(rf_class))
rf_class$ASVID <- row.names(rf_class)
row.names(rf_order) <- gsub("\\|.*", "", row.names(rf_order))
rf_order$ASVID <- row.names(rf_order)
row.names(rf_family) <- gsub("\\|.*", "", row.names(rf_family))
rf_family$ASVID <- row.names(rf_family)

# 合并数据
lambda_rf_file <- lamda_file %>%
  left_join(core_taxonomy[,c("FeatureID", "ASVID", "ASVLabel", "Class", "Order", "Family", "Genus")], by = "FeatureID") %>%
  left_join(rf_class[,c("ASVID", "Group", "P.adj")], by = "ASVID") %>%
  left_join(rf_order[,c("ASVID", "Group", "P.adj")], by = "ASVID") %>%
  left_join(rf_family[,c("ASVID", "Group", "P.adj")], by = "ASVID") %>%
  column_to_rownames("ASVID")

# 重命名列以区分分类水平
names(lambda_rf_file)[13:18] <- c("Tree_Class","rf_FDR_Class", "Tree_Order", "rf_FDR_Order", "Tree_Family", "rf_FDR_Family")

write.csv(lambda_rf_file, file = "./stats/fun_core_asv_qa_log_lambda_indicators.csv")

###################################################################################################
###################################################################################################

# 进一步筛选：ASV随机森林鉴定的植物科属随机森林鉴定的同一目，定义为更严格的特定植物目共进化的ASV

metadata_tree <- read.csv('../../../metadata/tree_metadata_merge_info.csv')
fun_order_coevoled_asv <- lambda_rf_file[!is.na(lambda_rf_file$Tree_Order) & lambda_rf_file$padj < 0.05,]
fun_order_coevoled_top_asv <- merge(fun_order_coevoled_asv, unique(metadata_tree[,c("Family","Order")]),
                                    by.x = "Tree_Family", by.y = "Family", all.x = T)
fun_order_coevoled_top_asv <- fun_order_coevoled_top_asv[fun_order_coevoled_top_asv$Tree_Order == fun_order_coevoled_top_asv$Order.y,]
names(fun_order_coevoled_top_asv)
fun_order_coevoled_top_asv <- fun_order_coevoled_top_asv[,c(2:17,1,18)]
fun_order_coevoled_top_asv <- subset(fun_order_coevoled_top_asv, Tree_Order != "Others")
# 71个有发育信号的ASV对应的植物科属于同一目 (这是很重要排除假阳性的筛选)

# 保存这个筛选结果
fun_lambda_rf_asv <- lambda_rf_file
fun_lambda_rf_asv$Coevolved <- NA  # 初始化新列
matching <- fun_lambda_rf_asv$FeatureID %in% fun_order_coevoled_top_asv$FeatureID
fun_lambda_rf_asv$Coevolved[matching] <- fun_lambda_rf_asv$Tree_Order[matching]
write.csv(fun_lambda_rf_asv, file = "./stats/TableS7B-revised-fun_ASV_qa_log_lamda_indicators.csv") 

###################################################################################################
###################################################################################################


### 绘制韦恩图 --------------------------------------------------------------
# 获取有显著系统发育信号的ASV列表 (padj < 0.05)
lamda_sig_list <- row.names(lambda_rf_file[lambda_rf_file$padj < 0.05, ])
# 166个ASV; 15.30% (166/1085) of fungal

# Class水平
class_sig_list <- rf_class$ASVID
class_venn <- Venn(list('Phylogenetically conserved ASVs' = lamda_sig_list,
                        'Class indicators' = class_sig_list))

win.metafile("./plots/Fungal_class_venn.emf")
plot(class_venn, doWeights = TRUE, type="circles",
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

# Order水平
order_sig_list <- rf_order$ASVID
order_venn <- Venn(list('Phylogenetically conserved ASVs' = lamda_sig_list,
                        'Order indicators' = order_sig_list))

win.metafile("./plots/Fungal_order_venn.emf")
plot(order_venn, doWeights = TRUE, type="circles",
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

# Family水平
family_sig_list <- rf_family$ASVID
family_venn <- Venn(list('Phylogenetically conserved ASVs' = lamda_sig_list,
                         'Family indicators' = family_sig_list))

win.metafile("./plots/Fungal_family_venn.emf")
plot(family_venn, doWeights = TRUE, type="circles",
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()



### 统计并绘制Class水平的分布条形图 ------------------------------------------
lambda_rf_file <- read.csv(file = "./stats/fun_core_asv_qa_log_lambda_indicators.csv")
top_class_bac <- read.csv('../../../feature_table/ITS_top10Class_color_meanValue.csv')

raw_df <- lambda_rf_file
# 处理Class分类，将不在前十的归为Others
raw_df$Taxa <- raw_df$Class
raw_df$Taxa <- ifelse(raw_df$Taxa == "", "Unassigned", raw_df$Taxa)
raw_df$Taxa <- ifelse(raw_df$Taxa %in% top_class_bac$Taxa, raw_df$Taxa, 'Others')

# 统计各组的ASV数量
all_df <- data.frame(t(table(raw_df$Taxa)))
lambda_df <- data.frame(t(table(raw_df[raw_df$padj < 0.05, 'Taxa'])))
rf_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]), 'Taxa'])))
co_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]) & raw_df$padj < 0.05, 'Taxa'])))
rf_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]), 'Taxa'])))
co_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]) & raw_df$padj < 0.05, 'Taxa'])))
rf_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]), 'Taxa'])))
co_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]) & raw_df$padj < 0.05, 'Taxa'])))

# 合并统计数据
res_df <- all_df %>%
  left_join(lambda_df, by = "Var2") %>%
  left_join(rf_df_class, by = "Var2") %>%
  left_join(co_df_class, by = "Var2") %>%
  left_join(rf_df_order, by = "Var2") %>%
  left_join(co_df_order, by = "Var2") %>%
  left_join(rf_df_family, by = "Var2") %>%
  left_join(co_df_family, by = "Var2")

res_df <- res_df[,c(2,3,5,7,9,11,13,15,17)]
res_df[is.na(res_df)] <- 0
names(res_df) <- c("Taxonomy", "All", "Phylogenetically conserved",
                   "Class indicators", "Class co-evolved",
                   "Order indicators", "Order co-evolved",
                   "Family indicators", "Family co-evolved")

# 转换数据为长格式用于绘图
fin_df <- melt(res_df, id.vars = 'Taxonomy')

# 设置颜色
color_df <- data.frame(
  Taxonomy = top_class_bac$Taxa,
  Color = top_class_bac$Color
)
color_df <- color_df[color_df$Taxonomy %in% fin_df$Taxonomy,]
fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)
color_manual <- color_df$Color

# 设置图形主题
theme_used <-
  theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 7, color = 'black'),
        axis.line = element_blank(),
        axis.text = element_text(size = 6, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 7, color = 'black'),
        legend.text = element_text(size = 6,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        legend.position = 'right')

# 绘制百分比堆叠条形图
p <- ggplot(fin_df, aes(x = variable, y = value, fill = Taxonomy)) +
  geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
  labs(x = '', y = 'Proportion') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_manual, name = "Class") +
  theme_bw() +
  theme_used

# 保存图形
width = 5.5 * 2; height = 6.5
name = './plots/Class distribution co-evolved Fungi ASVs'
ggsave(paste0(name, '.pdf'), p, width = width, height = height, units = 'cm')
ggsave(paste0(name, '.png'), p, width = width, height = height, units = 'cm')
