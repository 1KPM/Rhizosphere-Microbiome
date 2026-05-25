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


plant_colors <- read.csv("./data/Tree_top10_order_color.csv")
core_annotation <- read.csv('./data/K_gene_name.csv', row.names = 2)
names(core_annotation)
core_annotation <- core_annotation[,-1]
core_annotation$KOLabel <- paste(row.names(core_annotation), core_annotation$gene, sep = "|")
core_annotation <- core_annotation[,c(3,1,2)]
metadata_tree <- read.csv(file = "./data/metadata.csv", row.names = 1)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")


set.seed(2024)
ko_qa <- read.csv(file = "../../kegg_abundance/absolute/Bacteria_core0.2_quantitative.csv", row.names = 1,check.names = F)
datasets <- list(
  list(name = "bac_ko_qa_log", otu_table = log(ko_qa+1))
)

for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  
  dataset <- microtable$new(
    sample_table = metadata_tree,
    otu_table = otu_table,
    tax_table = core_annotation
  )
  
  dataset$tidy_dataset()
  
  dataset$cal_abund(rel = FALSE)
  
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "metastat",
    group = "Order",
    taxa_level = "KOLabel"
  )
  
  write.csv(t1$res_diff, file = paste0("./stats/",data_name, "_tree_metastat_res_diff_KO_core0.2.csv"))
  write.csv(t1$abund_table, file = paste0("./stats/",data_name, "_tree_metastat_abund_table_KO_core0.2.csv"))
  write.csv(t1$res_abund, file = paste0("./stats/",data_name, "_tree_metastat_res_abund_KO_core0.2.csv"))
  
  message("已完成处理: ", data_name)
}


ko_qa <- read.csv(file = "../../kegg_abundance/absolute/Fungi_core0.2_quantitative.csv", row.names = 1,check.names = F)
datasets <- list(
  list(name = "fun_ko_qa_log", otu_table = log(ko_qa+1))
)

for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  
  dataset <- microtable$new(
    sample_table = metadata_tree,
    otu_table = otu_table,
    tax_table = core_annotation
  )
  
  dataset$tidy_dataset()
  
  dataset$cal_abund(rel = FALSE)
  
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "metastat",
    group = "Order",
    taxa_level = "KOLabel"
  )
  
  write.csv(t1$res_diff, file = paste0("./stats/",data_name, "_tree_metastat_res_diff_KO_core0.2.csv"))
  write.csv(t1$abund_table, file = paste0("./stats/",data_name, "_tree_metastat_abund_table_KO_core0.2.csv"))
  write.csv(t1$res_abund, file = paste0("./stats/",data_name, "_tree_metastat_res_abund_KO_core0.2.csv"))
  
  message("已完成处理: ", data_name)
}

set.seed(2024)
ko_qa <- read.csv(file = "../../kegg_abundance/absolute/Protists_core0.2_quantitative.csv", row.names = 1,check.names = F)
datasets <- list(
  list(name = "pro_ko_qa_log", otu_table = log(ko_qa+1))
)

for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  
  dataset <- microtable$new(
    sample_table = metadata_tree,
    otu_table = otu_table,
    tax_table = core_annotation
  )
  
  dataset$tidy_dataset()
  
  dataset$cal_abund(rel = FALSE)
  
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "metastat",
    group = "Order",
    taxa_level = "KOLabel"
  )
  
  write.csv(t1$res_diff, file = paste0("./stats/",data_name, "_tree_metastat_res_diff_KO_core0.2.csv"))
  write.csv(t1$abund_table, file = paste0("./stats/",data_name, "_tree_metastat_abund_table_KO_core0.2.csv"))
  write.csv(t1$res_abund, file = paste0("./stats/",data_name, "_tree_metastat_res_abund_KO_core0.2.csv"))
  
  message("已完成处理: ", data_name)
}
