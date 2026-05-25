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

#### 下面分别使用Tree ID和Sample ID进行了随机森林分析，基于分析结果，最终选择Sample ID的结果进行展示

#### (1.1) 基于sample id，没有对树id进行平均

####数据准备####
core_taxonomy <- read.csv('../rawdata/amplicon/taxonomy/All_core_ASV_taxonomy.csv', row.names = 2)
names(core_taxonomy)
core_taxonomy <- core_taxonomy[,c(15,2,4)]

metadata_rs <- read.csv('../rawdata/amplicon/metadata/rhizosphere_metadata_merge_info.csv') #与nifH不是一个metadata
row.names(metadata_rs) <- metadata_rs$FileID
metadata_rs$Order <- ifelse(metadata_rs$Order %in% plant_colors$Order, metadata_rs$Order, "Others")
###其它分类水平
class_used <- read.csv("../rawdata/amplicon/Tree_top_class_color.csv")
order_used <- read.csv("../rawdata/amplicon/Tree_top_order_color.csv")
family_used <- read.csv("../rawdata/amplicon/Tree_top_family_color.csv")
genus_used <- read.csv("../rawdata/amplicon/Tree_top_genus_color.csv")

metadata_rs$Class <- ifelse(metadata_rs$Class %in% class_used$Class, metadata_rs$Class, "Others")
metadata_rs$Family <- ifelse(metadata_rs$Family %in% family_used$Family, metadata_rs$Family, "Others")
metadata_rs$Genus <- ifelse(metadata_rs$Genus %in% genus_used$Genus, metadata_rs$Genus, "Others")


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

# 定义要处理的数据集列表, qa数据
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
    sample_table = metadata_rs,
    otu_table = otu_table,
    tax_table = core_taxonomy
  )
  
  # 整理数据集
  dataset$tidy_dataset()
  
  # 计算丰度（不进行标准化）
  dataset$cal_abund(rel = FALSE) 
  
  # 循环处理每个分组变量
  for (group_var in c("Class", "Family", "Genus")) {
    # 执行随机森林差异分析
    t1 <- trans_diff$new(
      dataset = dataset,
      method = "rf",
      group = group_var,
      taxa_level = "ASVLabel" 
    )
    
    # 保存结果文件（添加分组变量到文件名）
    write.csv(t1$res_diff, file = paste0(data_name, "_", group_var, "_rf_res_diff_ASV.csv"))
    write.csv(t1$abund_table, file = paste0(data_name, "_", group_var, "_rf_abund_table_ASV.csv"))
    write.csv(t1$res_abund, file = paste0(data_name, "_", group_var, "_rf_res_abund_ASV.csv"))
    
    # 打印详细的进度信息
    message("已完成处理: 数据集 = ", data_name, " | 分组变量 = ", group_var)
  }
}


# 定义要处理的数据集列表, ra数据
datasets <- list(
  list(name = "bac_ra", otu_table = bac_ra),
  list(name = "fun_ra", otu_table = fun_ra),
  list(name = "pro_ra", otu_table = pro_ra)
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
  dataset$cal_abund(rel = FALSE) 
  
  # 循环处理每个分组变量
  for (group_var in c("Class", "Family", "Genus")) {
    # 执行随机森林差异分析
    t1 <- trans_diff$new(
      dataset = dataset,
      method = "rf",
      group = group_var,
      taxa_level = "ASVLabel" 
    )
    
    # 保存结果文件（添加分组变量到文件名）
    write.csv(t1$res_diff, file = paste0(data_name, "_", group_var, "_rf_res_diff_ASV.csv"))
    write.csv(t1$abund_table, file = paste0(data_name, "_", group_var, "_rf_abund_table_ASV.csv"))
    write.csv(t1$res_abund, file = paste0(data_name, "_", group_var, "_rf_res_abund_ASV.csv"))
    
    # 打印详细的进度信息
    message("已完成处理: 数据集 = ", data_name, " | 分组变量 = ", group_var)
  }
}
##################################

#########################################taxa水平分析
### （2.2）使用全部的taxa进行分析, 基于sample id
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


# 定义要处理的数据集列表, qa数据
datasets <- list(
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
  
  # 循环处理每个分组变量
  for (group_var in c("Class", "Family", "Genus")) {
  
  # 执行随机森林差异分析
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "rf",
    group = group_var,
    taxa_level = "Taxa"
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0(data_name, "_", group_var, "_rf_res_diff_taxa.csv"))
  write.csv(t1$abund_table, file = paste0(data_name, "_", group_var, "_rf_abund_table_taxa.csv"))
  write.csv(t1$res_abund, file = paste0(data_name, "_", group_var, "_rf_res_abund_taxa.csv"))
  
  # 打印进度信息
  message("已完成处理: 数据集 = ", data_name, " | 分组变量 = ", group_var)
  }
}



# 定义要处理的数据集列表, ra数据
datasets <- list(
  list(name = "bac_ra_taxa", otu_table = bac_ra_taxa[,-c(1:2)], tax_table = bac_ra_taxa[,1:2]),
  list(name = "fun_ra_taxa", otu_table = fun_ra_taxa[,-c(1:2)], tax_table = fun_ra_taxa[,1:2]),
  list(name = "pro_ra_taxa", otu_table = pro_ra_taxa[,-c(1:2)], tax_table = pro_ra_taxa[,1:2])
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
  
  # 循环处理每个分组变量
  for (group_var in c("Class", "Family", "Genus")) {
    
    # 执行随机森林差异分析
    t1 <- trans_diff$new(
      dataset = dataset,
      method = "rf",
      group = group_var,
      taxa_level = "Taxa"
    )
    
    # 保存结果文件
    write.csv(t1$res_diff, file = paste0(data_name, "_", group_var, "_rf_res_diff_taxa.csv"))
    write.csv(t1$abund_table, file = paste0(data_name, "_", group_var, "_rf_abund_table_taxa.csv"))
    write.csv(t1$res_abund, file = paste0(data_name, "_", group_var, "_rf_res_abund_taxa.csv"))
    
    # 打印进度信息
    message("已完成处理: 数据集 = ", data_name, " | 分组变量 = ", group_var)
  }
}

