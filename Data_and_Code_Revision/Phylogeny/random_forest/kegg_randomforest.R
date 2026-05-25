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


#### (1) KO

####数据准备####
plant_colors <- read.csv("../kegg_abundance/Tree_top_order_color.csv")

core_annotation <- read.csv('../kegg_abundance/K_gene_name.csv', row.names = 2)
names(core_annotation)
core_annotation <- core_annotation[,-1]
core_annotation$KOLabel <- paste(row.names(core_annotation), core_annotation$gene, sep = "|")
core_annotation <- core_annotation[,c(3,1,2)]

#####
metadata_tree <- read.csv(file = "../kegg_abundance/tree_metadata_merge_info.csv", row.names = 1)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")

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

###
datasets <- list(
  list(name = "pro_ko_core0.2_qa_log", otu_table = log(pro_ko_core0.2_qa+1)),
  list(name = "fun_ko_core0.2_qa_log", otu_table = log(fun_ko_core0.2_qa+1)),
  list(name = "bac_ko_core0.2_qa_log", otu_table = log(bac_ko_core0.2_qa+1)),
  # list(name = "pro_ko_core0.2_qa", otu_table = pro_ko_core0.2_qa),
  # list(name = "fun_ko_core0.2_qa", otu_table = fun_ko_core0.2_qa),
  # list(name = "bac_ko_core0.2_qa", otu_table = bac_ko_core0.2_qa),
  list(name = "pro_ko_core0.2_tpm", otu_table = pro_ko_core0.2_tpm),
  list(name = "fun_ko_core0.2_tpm", otu_table = fun_ko_core0.2_tpm),
  list(name = "bac_ko_core0.2_tpm", otu_table = bac_ko_core0.2_tpm)
  )

# 循环处理每个数据集
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
  t1 <- trans_diff$new(
    dataset = dataset,
    method = "rf",
    alpha = 0.2, #该参数输出只保留矫正P值小于0.2的KO
    group = "Order",
    taxa_level = "KOLabel"
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0("./metagenome/", data_name, "_tree_rf_res_diff_KO.csv"))
  write.csv(t1$abund_table, file = paste0("./metagenome/", data_name, "_tree_rf_abund_table_KO.csv"))
  write.csv(t1$res_abund, file = paste0("./metagenome/", data_name, "_tree_rf_res_abund_KO.csv"))
  
  # 打印进度信息
  message("已完成处理: ", data_name)
}


################(2) pathway

pathway_annotation <- read.csv('../kegg_abundance/kegg_pathway_for_1KPM.csv')
names(pathway_annotation)
pathway_annotation <- pathway_annotation[,c("L3","pathway")]
pathway_annotation <- unique.data.frame(pathway_annotation)
row.names(pathway_annotation) <- gsub("\\[PATH.*", "", pathway_annotation$pathway)

pathway_annotation <- pathway_annotation[,c(2,1)]

###relative abundance
bac_pathway_tpm <- read.csv(file = "../kegg_abundance/relative/Summary_bac_tpm_all_pathway.csv", row.names = 1)
names(bac_pathway_tpm) <- gsub("RS.", "", names(bac_pathway_tpm))
bac_pathway_tpm <- bac_pathway_tpm[-c(1:26), ]
row.names(bac_pathway_tpm) <- gsub("\\[PATH.*", "", row.names(bac_pathway_tpm))
row.names(bac_pathway_tpm) <- gsub(".*\\|", "", row.names(bac_pathway_tpm))

bac_pathway_qa <- read.csv(file = "../kegg_abundance/absolute/Summary_bac_qa_all_pathway.csv", row.names = 1)
names(bac_pathway_qa) <- gsub("RS.", "", names(bac_pathway_qa))
bac_pathway_qa <- bac_pathway_qa[-c(1:26), ]
row.names(bac_pathway_qa) <- gsub("\\[PATH.*", "", row.names(bac_pathway_qa))
row.names(bac_pathway_qa) <- gsub(".*\\|", "", row.names(bac_pathway_qa))

###
fun_pathway_tpm <- read.csv(file = "../kegg_abundance/relative/Summary_fun_tpm_all_pathway.csv", row.names = 1)
names(fun_pathway_tpm) <- gsub("RS.", "", names(fun_pathway_tpm))
fun_pathway_tpm <- fun_pathway_tpm[-c(1:26), ]
row.names(fun_pathway_tpm) <- gsub("\\[PATH.*", "", row.names(fun_pathway_tpm))
row.names(fun_pathway_tpm) <- gsub(".*\\|", "", row.names(fun_pathway_tpm))

fun_pathway_qa <- read.csv(file = "../kegg_abundance/absolute/Summary_fun_qa_all_pathway.csv", row.names = 1)
names(fun_pathway_qa) <- gsub("RS.", "", names(fun_pathway_qa))
fun_pathway_qa <- fun_pathway_qa[-c(1:26), ]
row.names(fun_pathway_qa) <- gsub("\\[PATH.*", "", row.names(fun_pathway_qa))
row.names(fun_pathway_qa) <- gsub(".*\\|", "", row.names(fun_pathway_qa))


###
pro_pathway_tpm <- read.csv(file = "../kegg_abundance/relative/Summary_pro_tpm_all_pathway.csv", row.names = 1)
names(pro_pathway_tpm) <- gsub("RS.", "", names(pro_pathway_tpm))
pro_pathway_tpm <- pro_pathway_tpm[-c(1:26), ]
row.names(pro_pathway_tpm) <- gsub("\\[PATH.*", "", row.names(pro_pathway_tpm))
row.names(pro_pathway_tpm) <- gsub(".*\\|", "", row.names(pro_pathway_tpm))

pro_pathway_qa <- read.csv(file = "../kegg_abundance/absolute/Summary_pro_qa_all_pathway.csv", row.names = 1)
names(pro_pathway_qa) <- gsub("RS.", "", names(pro_pathway_qa))
pro_pathway_qa <- pro_pathway_qa[-c(1:26), ]
row.names(pro_pathway_qa) <- gsub("\\[PATH.*", "", row.names(pro_pathway_qa))
row.names(pro_pathway_qa) <- gsub(".*\\|", "", row.names(pro_pathway_qa))

###
datasets <- list(
  list(name = "bac_pathway_tpm", otu_table = bac_pathway_tpm),
  list(name = "fun_pathway_tpm", otu_table = fun_pathway_tpm),
  list(name = "pro_pathway_tpm", otu_table = pro_pathway_tpm),
  
  list(name = "bac_pathway_qa", otu_table = bac_pathway_qa),
  list(name = "fun_pathway_qa", otu_table = fun_pathway_qa),
  list(name = "pro_pathway_qa", otu_table = pro_pathway_qa),
  
  list(name = "bac_pathway_qa_log", otu_table = log(bac_pathway_qa+1)),
  list(name = "fun_pathway_qa_log", otu_table = log(fun_pathway_qa+1)),
  list(name = "pro_pathway_qa_log", otu_table = log(pro_pathway_qa+1))
)


# 循环处理每个数据集
for (data_info in datasets) {
  data_name <- data_info$name
  otu_table <- data_info$otu_table
  
  # 创建microtable对象
  dataset <- microtable$new(
    sample_table = metadata_tree,
    otu_table = otu_table,
    tax_table = pathway_annotation
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
    taxa_level = "pathway",
    p_adjust_method = "none"
  )
  
  # 保存结果文件
  write.csv(t1$res_diff, file = paste0("./metagenome/", data_name, "_tree_rf_res_diff_Pathway.csv"))
  write.csv(t1$abund_table, file = paste0("./metagenome/", data_name, "_tree_rf_abund_table_Pathway.csv"))
  write.csv(t1$res_abund, file = paste0("./metagenome/", data_name, "_tree_rf_res_abund_Pathway.csv"))
  
  # 绘制MeanDecreaseGini条形图
  g1 <- t1$plot_diff_bar(
    use_number = 1:20,
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
    filename = paste0(data_name, "_tree_MeanDecreaseGini_abundance_barplot_Pathway.pdf"),
    plot = g12,
    width = 18,
    height = 9
  )
  
  # 打印进度信息
  message("已完成处理: ", data_name)
}



