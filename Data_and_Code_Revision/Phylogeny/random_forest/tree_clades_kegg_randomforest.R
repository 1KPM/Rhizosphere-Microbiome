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


###分类水平
class_used <- read.csv("../metadata/Tree_top_class_color.csv")
order_used <- read.csv("../metadata/Tree_top_order_color.csv")
family_used <- read.csv("../metadata/Tree_top_family_color.csv")
genus_used <- read.csv("../metadata/Tree_top_genus_color.csv")


metadata_tree <- read.csv(file = "../metadata/tree_metadata_merge_info.csv", row.names = 1)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% order_used$Order, metadata_tree$Order, "Others")
metadata_tree$Class <- ifelse(metadata_tree$Class %in% class_used$Class, metadata_tree$Class, "Others")
metadata_tree$Family <- ifelse(metadata_tree$Family %in% family_used$Family, metadata_tree$Family, "Others")
metadata_tree$Genus <- ifelse(metadata_tree$Genus %in% genus_used$Genus, metadata_tree$Genus, "Others")
table(metadata_tree$Genus)
table(metadata_tree$Family)

core_annotation <- read.csv('../kegg_abundance/K_gene_name.csv', row.names = 2)
names(core_annotation)
core_annotation <- core_annotation[,-1]
core_annotation$KOLabel <- paste(row.names(core_annotation), core_annotation$gene, sep = "|")
core_annotation <- core_annotation[,c(3,1,2)]

bac_ko_core0.2_tpm <- read.csv(file = "../kegg_abundance/relative/Bacteria_core0.2_tpm.csv", row.names = 1, header = T)
names(bac_ko_core0.2_tpm) <- gsub("RS.", "", names(bac_ko_core0.2_tpm))
bac_ko_core0.2_qa <- read.csv(file = "../kegg_abundance/absolute/Bacteria_core0.2_quantitative.csv", row.names = 1, header = T)
names(bac_ko_core0.2_qa) <- gsub("RS.", "", names(bac_ko_core0.2_qa))

fun_ko_core0.2_tpm <- read.csv(file = "../kegg_abundance/relative/Fungi_core0.2_tpm.csv", row.names = 1, header = T)
names(fun_ko_core0.2_tpm) <- gsub("RS.", "", names(fun_ko_core0.2_tpm))
fun_ko_core0.2_qa <- read.csv(file = "../kegg_abundance/absolute/Fungi_core0.2_quantitative.csv", row.names = 1, header = T)
names(fun_ko_core0.2_qa) <- gsub("RS.", "", names(fun_ko_core0.2_qa))

pro_ko_core0.2_tpm <- read.csv(file = "../kegg_abundance/relative/Protists_core0.2_tpm.csv", row.names = 1, header = T)
names(pro_ko_core0.2_tpm) <- gsub("RS.", "", names(pro_ko_core0.2_tpm))
pro_ko_core0.2_qa <- read.csv(file = "../kegg_abundance/absolute/Protists_core0.2_quantitative.csv", row.names = 1, header = T)
names(pro_ko_core0.2_qa) <- gsub("RS.", "", names(pro_ko_core0.2_qa))


#############################################
###
datasets <- list(
  list(name = "pro_ko_qa_log", otu_table = log(pro_ko_core0.2_qa+1)),
  list(name = "fun_ko_qa_log", otu_table = log(fun_ko_core0.2_qa+1)),
  list(name = "bac_ko_qa_log", otu_table = log(bac_ko_core0.2_qa+1)),
  list(name = "pro_ko_ra", otu_table = pro_ko_core0.2_tpm)
  list(name = "fun_ko_ra", otu_table = fun_ko_core0.2_tpm),
  list(name = "bac_ko_ra", otu_table = bac_ko_core0.2_tpm)
)

for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  
  # 创建microtable对象
  dataset <- microtable$new(
    sample_table = metadata_tree,
    otu_table = otu_table,
    tax_table = core_annotation
  )
  
  # 整理数据集
  dataset$tidy_dataset()
  
  # 计算丰度（不进行标准化）
  dataset$cal_abund(rel = FALSE)
  
  # 执行随机森林差异分析
  for (group_var in c("Class", "Family")) { #genus 样本太少，结果出问题
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "rf",
    alpha = 0.5,
    group = group_var,
    taxa_level = "KOLabel"
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0(data_name, "_", group_var, "_tree_rf_res_diff_KO.csv"))
  write.csv(t1$abund_table, file = paste0(data_name, "_", group_var, "_tree_rf_abund_table_KO.csv"))
  write.csv(t1$res_abund, file = paste0(data_name, "_", group_var, "_tree_rf_res_abund_KO.csv"))
 
  message("已完成处理: ", data_name)
  }
}

















