
### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

####randomForest—biomarkers

library(varSelRF)
library(pROC)
library(randomForest)
library(vegan)
library(rbiom)
library(microeco)
library(magrittr)
library(ggplot2)
library(ggpubr)

###################（1）ASV水平随机森林分析
####数据准备####
plant_colors <- read.csv("../rawdata/amplicon/Tree_top_order_color.csv")

core_taxonomy <- read.csv('../rawdata/amplicon/taxonomy/All_core_ASV_taxonomy.csv', row.names = 2)
names(core_taxonomy)
core_taxonomy <- core_taxonomy[,c(15,2,4)]

metadata_rs <- read.csv('../rawdata/amplicon/metadata/rhizosphere_metadata_merge_info.csv')
row.names(metadata_rs) <- metadata_rs$FileID
metadata_rs$Order <- ifelse(metadata_rs$Order %in% plant_colors$Order, metadata_rs$Order, "Others")

###
bac_ra <- read.table('../rawdata/amplicon/feature_table/16S/core_feature_table_relative.tsv',header = T, row.names = 1)
fun_ra <- read.table('../rawdata/amplicon/feature_table/ITS/core_feature_table_relative.tsv',header = T, row.names = 1)
pro_ra <- read.table('../rawdata/amplicon/feature_table/Protist/core_feature_table_relative.tsv',header = T, row.names = 1)

bac_qa <- read.table('../rawdata/amplicon/feature_table/16S/core_feature_table_absolute.tsv',header = T, row.names = 1)
fun_qa <- read.table('../rawdata/amplicon/feature_table/ITS/core_feature_table_absolute.tsv',header = T, row.names = 1)
pro_qa <- read.table('../rawdata/amplicon/feature_table/Protist/core_feature_table_absolute.tsv',header = T, row.names = 1)

bac_qa_log <- read.table('../rawdata/amplicon/feature_table/16S/core_feature_table_absolute.tsv',header = T, row.names = 1)
bac_qa_log <- log(bac_qa_log+1)
fun_qa_log <- read.table('../rawdata/amplicon/feature_table/ITS/core_feature_table_absolute.tsv',header = T, row.names = 1)
fun_qa_log <- log(fun_qa_log+1)
pro_qa_log <- read.table('../rawdata/amplicon/feature_table/Protist/core_feature_table_absolute.tsv',header = T, row.names = 1)
pro_qa_log <- log(pro_qa_log+1)

##测试： otu_table = pro_qa[1:100, ]; data_name = "pro_qa"
###循环
# 定义要处理的数据集列表
datasets <- list(
  list(name = "bac_ra", otu_table = bac_ra),
  list(name = "fun_ra", otu_table = fun_ra),
  list(name = "pro_ra", otu_table = pro_ra),
  list(name = "bac_qa", otu_table = bac_qa),
  list(name = "fun_qa", otu_table = fun_qa),
  list(name = "pro_qa", otu_table = pro_qa),
  list(name = "bac_qa_log", otu_table = bac_qa_log),
  list(name = "fun_qa_log", otu_table = fun_qa_log),
  list(name = "pro_qa_log", otu_table = pro_qa_log)
)

# 循环处理每个数据集
for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  
  # 创建microtable对象
  dataset <- microtable$new(
    sample_table = metadata_rs,
    otu_table = otu_table,
    tax_table = core_taxonomy
  )
  
  # 整理数据集
  dataset$tidy_dataset()
  
  # 计算丰度（不进行标准化）
  dataset$cal_abund(rel = FALSE) #该函数会计算每个colname的丰度，所以taxonomy没必要放那么多变量（colname）
  
  # 执行随机森林差异分析
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "rf",
    group = "Order",
    taxa_level = "ASVLabel" #经常重新命名为各colname的拼接，暂时找不到解决方案, 使用taxonomy的第一行貌似不会拼接
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0(data_name, "_rf_res_diff_ASV.csv"))
  write.csv(t1$abund_table, file = paste0(data_name, "_rf_abund_table_ASV.csv"))
  write.csv(t1$res_abund, file = paste0(data_name, "_rf_res_abund_ASV.csv"))
  
  # 打印进度信息
  message("已完成处理: ", data_name)
}

######################################
###（2）taxa水平随机森林分析
bac_ra_taxa <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/16S_barplot_data_relative_sorted_final.csv',header = T)
row.names(bac_ra_taxa) <- paste(tolower(substr(bac_ra_taxa$Level,1,1)), bac_ra_taxa$Taxa, sep = "_")
bac_ra_taxa$Taxa <- row.names(bac_ra_taxa)
fun_ra_taxa <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/ITS_barplot_data_relative_sorted_final.csv',header = T)
row.names(fun_ra_taxa) <- paste(tolower(substr(fun_ra_taxa$Level,1,1)), fun_ra_taxa$Taxa, sep = "_")
fun_ra_taxa$Taxa <- row.names(fun_ra_taxa)
pro_ra_taxa <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/Protist_barplot_data_relative_sorted_final.csv',header = T)
row.names(pro_ra_taxa) <- paste(tolower(substr(pro_ra_taxa$Level,1,1)), pro_ra_taxa$Taxa, sep = "_")
pro_ra_taxa$Taxa <- row.names(pro_ra_taxa)

###
bac_qa_taxa <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/16S_barplot_data_absolute_sorted_final.csv', header = T)
row.names(bac_qa_taxa) <- paste(tolower(substr(bac_qa_taxa$Level,1,1)), bac_qa_taxa$Taxa, sep = "_")
bac_qa_taxa$Taxa <- row.names(bac_qa_taxa)
fun_qa_taxa <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/ITS_barplot_data_absolute_sorted_final.csv',header = T)
row.names(fun_qa_taxa) <- paste(tolower(substr(fun_qa_taxa$Level,1,1)), fun_qa_taxa$Taxa, sep = "_")
fun_qa_taxa$Taxa <- row.names(fun_qa_taxa)
pro_qa_taxa <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/Protist_barplot_data_absolute_sorted_final.csv',header = T)
row.names(pro_qa_taxa) <- paste(tolower(substr(pro_qa_taxa$Level,1,1)), pro_qa_taxa$Taxa, sep = "_")
pro_qa_taxa$Taxa <- row.names(pro_qa_taxa)

###
bac_qa_taxa_log <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/16S_barplot_data_absolute_sorted_final.csv', header = T)
row.names(bac_qa_taxa_log) <- paste(tolower(substr(bac_qa_taxa_log$Level,1,1)), bac_qa_taxa_log$Taxa, sep = "_")
bac_qa_taxa_log_taxa <- bac_qa_taxa_log[,1:2]
bac_qa_taxa_log_taxa$Taxa <- row.names(bac_qa_taxa_log_taxa)
bac_qa_taxa_log <- log(bac_qa_taxa_log[,-c(1:2)] +1)

fun_qa_taxa_log <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/ITS_barplot_data_absolute_sorted_final.csv',header = T)
row.names(fun_qa_taxa_log) <- paste(tolower(substr(fun_qa_taxa_log$Level,1,1)), fun_qa_taxa_log$Taxa, sep = "_")
fun_qa_taxa_log_taxa <- fun_qa_taxa_log[,1:2]
fun_qa_taxa_log_taxa$Taxa <- row.names(fun_qa_taxa_log_taxa)
fun_qa_taxa_log <- log(fun_qa_taxa_log[,-c(1:2)]+1)

pro_qa_taxa_log <- read.csv('../rawdata/amplicon/each_Taxa/sorted_data/Protist_barplot_data_absolute_sorted_final.csv',header = T)
row.names(pro_qa_taxa_log) <- paste(tolower(substr(pro_qa_taxa_log$Level,1,1)), pro_qa_taxa_log$Taxa, sep = "_")
pro_qa_taxa_log_taxa <- pro_qa_taxa_log[,1:2]
pro_qa_taxa_log_taxa$Taxa <- row.names(pro_qa_taxa_log_taxa)
pro_qa_taxa_log <- log(pro_qa_taxa_log[,-c(1:2)]+1)


# 定义要处理的数据集列表
datasets <- list(
  list(name = "bac_ra_taxa", otu_table = bac_ra_taxa[,-c(1:2)], tax_table = bac_ra_taxa[,1:2]),
  list(name = "fun_ra_taxa", otu_table = fun_ra_taxa[,-c(1:2)], tax_table = fun_ra_taxa[,1:2]),
  list(name = "pro_ra_taxa", otu_table = pro_ra_taxa[,-c(1:2)], tax_table = pro_ra_taxa[,1:2]),
  list(name = "bac_qa_taxa", otu_table = bac_qa_taxa[,-c(1:2)], tax_table = bac_qa_taxa[,1:2]),
  list(name = "fun_qa_taxa", otu_table = fun_qa_taxa[,-c(1:2)], tax_table = fun_qa_taxa[,1:2]),
  list(name = "pro_qa_taxa", otu_table = pro_qa_taxa[,-c(1:2)], tax_table = pro_qa_taxa[,1:2]),
  list(name = "bac_qa_taxa_log", otu_table = bac_qa_taxa_log, tax_table = bac_qa_taxa_log_taxa),
  list(name = "fun_qa_taxa_log", otu_table = fun_qa_taxa_log, tax_table = fun_qa_taxa_log_taxa),
  list(name = "pro_qa_taxa_log", otu_table = pro_qa_taxa_log, tax_table = pro_qa_taxa_log_taxa)
)

# 循环处理每个数据集
for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  tax_table <- data_info$tax_table
  
  # 创建microtable对象
  dataset <- microtable$new(
    sample_table = metadata_rs,
    otu_table = otu_table,
    tax_table = tax_table
  )
  
  # 整理数据集
  dataset$tidy_dataset()
  
  # 计算丰度（不进行标准化）
  dataset$cal_abund(rel = FALSE)
  
  # 执行随机森林差异分析
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "rf",
    group = "Order",
    taxa_level = "Taxa"
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0(data_name, "_rf_res_diff_taxa.csv"))
  write.csv(t1$abund_table, file = paste0(data_name, "_rf_abund_table_taxa.csv"))
  write.csv(t1$res_abund, file = paste0(data_name, "_rf_res_abund_taxa.csv"))
  
  # 绘制MeanDecreaseGini条形图
  g1 <- t1$plot_diff_bar(
    use_number = 1:50,
    group_order = plant_colors$Order
  ) + 
    theme(legend.position = "none")
  
  # 绘制对应分类单元的丰度图
  g2 <- t1$plot_diff_abund(
    group_order = plant_colors$Order,
    select_taxa = t1$plot_diff_bar_taxa
  ) + 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  
  # 合并图形并保存
  g12 <- gridExtra::grid.arrange(g1, g2, ncol = 2, widths = c(2, 1.7))
  ggsave(
    filename = paste0(data_name, "_MeanDecreaseGini_abundance_barplot_Taxa.pdf"),
    plot = g12,
    width = 18,
    height = 9
  )
  
  # 打印进度信息
  message("已完成处理: ", data_name)
}
