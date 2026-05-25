
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

####数据准备####
plant_colors <- read.csv("./data/Tree_top10_order_color.csv")

core_taxonomy <- read.csv('./data/All_core_ASV_taxonomy.csv', row.names = 2)
names(core_taxonomy)
core_taxonomy <- core_taxonomy[,c(15,2,4)]

metadata <- read.csv('./data/rhizosphere_metadata_merge_info.csv')
row.names(metadata) <- metadata$FileID
metadata$Order <- ifelse(metadata$Order %in% plant_colors$Order, metadata$Order, "Others")

###
bac_qa <- read.table('../../feature_table/16S/core_feature_table_absolute.tsv',header = T, row.names = 1)
fun_qa <- read.table('../../feature_table/ITS/core_feature_table_absolute.tsv',header = T, row.names = 1)
pro_qa <- read.table('../../feature_table/Protist/core_feature_table_absolute.tsv',header = T, row.names = 1)

bac_qa_log <- log(bac_qa+1)
fun_qa_log <- log(fun_qa+1)
pro_qa_log <- log(pro_qa+1)

# 定义要处理的数据集列表
datasets <- list(
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
    sample_table = metadata,
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
    method = "metastat",
    group = "Order",
    taxa_level = "ASVLabel" #经常重新命名为各colname的拼接，暂时找不到解决方案, 使用taxonomy的第一行貌似不会拼接
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0("./stats/", data_name, "_metastat_res_diff_ASV_core0.2.csv"))
  write.csv(t1$abund_table, file = paste0("./stats/", data_name, "_metastat_abund_table_ASV_core0.2.csv"))
  write.csv(t1$res_abund, file = paste0("./stats/", data_name, "_metastat_res_abund_ASV_core0.2.csv"))
  
  # 打印进度信息
  message("已完成处理: ", data_name)
}

