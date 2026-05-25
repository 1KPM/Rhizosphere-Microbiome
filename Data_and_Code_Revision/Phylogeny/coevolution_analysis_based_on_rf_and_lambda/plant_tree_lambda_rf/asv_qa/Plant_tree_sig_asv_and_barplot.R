### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

######################################################################
library(dplyr)
library(tibble)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(reshape2)
library(ape)

### 定义绘图函数
plot_coevolved_asv_tree <- function(microbe_type, lambda_file, qa_file, ra_file, 
                                    tree_file, metadata_file, color_file, 
                                    output_prefix, p_value_col = "padj") {
  
  ### 读取颜色文件
  plant_colors <- read.csv(color_file)
  
  ### 读取lambda_rf结果
  lambda_rf_asv <- read.csv(file = lambda_file, row.names = 1)
  
  ### 读取丰度数据
  qa_mean <- read.csv(qa_file, header = T, row.names = 1, check.names = F)
  ra_mean <- read.csv(ra_file, header = T, row.names = 1, check.names = F)
  
  qa_mean <- as.data.frame(t(qa_mean)) %>% rownames_to_column("FeatureID")
  ra_mean <- as.data.frame(t(ra_mean)) %>% rownames_to_column("FeatureID")
  
  ### 使用相对丰度绘图
  qa_mean <- ra_mean
  ### 只是展示的情况用ra的值而已, coevolved ASV是基于logqa鉴定的
  ### 因为有的样本绝对丰度太高，导致所有asv丰度都很高，
  ### 不适合绘制绝对丰度，使用相对丰度绘制
  
  
  ### 筛选共进化的ASV
  order_coevoled_asv <- lambda_rf_asv[!is.na(lambda_rf_asv$Tree_Order) & lambda_rf_asv$padj < 0.05, ]
  
  ### 计算各Order的总丰度
  order_coevoled_asv_qa <- merge.data.frame(
    order_coevoled_asv[, c("FeatureID", "Tree_Order")], 
    qa_mean, 
    by = "FeatureID"
  )
  
  order_coevoled_asv_qa_total <- aggregate(
    order_coevoled_asv_qa[, -c(1:2)], 
    by = list(order_coevoled_asv_qa$Tree_Order), 
    FUN = sum
  )
  
  order_coevoled_asv_qa_total <- order_coevoled_asv_qa_total %>% 
    column_to_rownames("Group.1") %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
  
  # 不是每个目都有共进化ASV，对目颜色文件进行筛选
  plant_colors <- plant_colors[plant_colors$Order %in% names(order_coevoled_asv_qa_total), ]
  
  order_coevoled_asv_qa_total <- order_coevoled_asv_qa_total[,  plant_colors$Order]
  names(order_coevoled_asv_qa_total) <- paste(names(order_coevoled_asv_qa_total), "_coevolved_ASV", sep = "")
  
  ### 读取并处理系统发育树
  tree <- read.tree(tree_file)
  drop_list <- tree$tip.label[-match(row.names(order_coevoled_asv_qa_total), tree$tip.label)]
  fin_tree <- drop.tip(tree, drop_list)
  fin_table <- order_coevoled_asv_qa_total[row.names(order_coevoled_asv_qa_total) %in% fin_tree$tip.label, ]
  
  ### 读取并处理元数据
  metadata <- read.csv(metadata_file)
  metadata <- subset(metadata, metadata$TreeID %in% fin_tree$tip.label)
  metadata$Order <- ifelse(metadata$Order %in% plant_colors$Order, metadata$Order, "Others")
  metadata$Order <- factor(metadata$Order, levels = plant_colors$Order)
  
  fin_table <- fin_table[fin_tree$tip.label, ]
  fin_table <- fin_table %>% rownames_to_column("label")
  
  ### 分组信息
  plant_groupInfo <- split(metadata$TreeID, metadata$Order)
  plant_tree_grouped <- groupOTU(fin_tree, plant_groupInfo, group_name = "Lineage")
  
  ### 绘制系统发育树
  plant_tree <- ggtree(
    plant_tree_grouped, 
    aes(color = Lineage), 
    size = 1, 
    ladderize = T, 
    linewidth = 0.1,
    layout = 'fan', 
    branch.length = "none", 
    right = TRUE, 
    open.angle = 180
  ) +
    scale_color_manual(values = plant_colors$Color2, breaks = plant_colors$Order) +
    new_scale_color()
  
  ### 定义柱状图颜色
  coevoled_asv_columns <- names(fin_table)[-1]
  color_vector <- setNames(plant_colors$Color2, coevoled_asv_columns)
  
  ### 循环添加柱状图
  final_plot <- plant_tree
  for (col in coevoled_asv_columns) {
    col_color <- color_vector[col]
    final_plot <- final_plot +
      geom_fruit(
        data = fin_table,
        geom = geom_bar,
        stat = 'identity',
        width = 2,
        aes_string(y = "label", x = col),
        pwidth = 0.1,
        fill = col_color
      )
  }
  
  ### 保存图形
  output_file <- paste0("./plant_tree_barplot/", output_prefix, ".pdf")
  ggsave(
    final_plot,
    filename = output_file,
    width = 12,
    height = 12,
    dpi = 900,
    bg = "white"
  )
  
  cat("Plot saved to:", output_file, "\n")
}

######################################################################
### 生成细菌图
plot_coevolved_asv_tree(
  microbe_type = "bacterial",
  lambda_file = "./stats/bac_core_asv_qa_log_lambda_indicators.csv",
  qa_file = "../../../feature_table/16S/core_feature_table_absolute_mean.csv",
  ra_file = "../../../feature_table/16S/core_feature_table_relative_mean.csv",
  tree_file = "../../../metadata/tree_metadata_merge_info_final_align_tree.nwk",
  metadata_file = "../../../metadata/tree_metadata_merge_info.csv",
  color_file = "../../../feature_table/Tree_top_order_color.csv",
  output_prefix = "plant_order_level_coevolved_bacterial_asv_loop"
)

######################################################################
### 生成真菌图
plot_coevolved_asv_tree(
  microbe_type = "fungal",
  lambda_file = "./stats/fun_core_asv_qa_log_lambda_indicators.csv",
  qa_file = "../../../feature_table/ITS/core_feature_table_absolute_mean.csv",
  ra_file = "../../../feature_table/ITS/core_feature_table_relative_mean.csv",
  tree_file = "../../../metadata/tree_metadata_merge_info_final_align_tree.nwk",
  metadata_file = "../../../metadata/tree_metadata_merge_info.csv",
  color_file = "../../../feature_table/Tree_top_order_color.csv",
  output_prefix = "plant_order_level_coevolved_fungal_asv_loop"
)

######################################################################
### 生成原生生物图
plot_coevolved_asv_tree(
  microbe_type = "protistan",
  lambda_file = "./stats/pro_core_asv_qa_log_lambda_indicators.csv",
  qa_file = "../../../feature_table/Protist/core_feature_table_absolute_mean.csv",
  ra_file = "../../../feature_table/Protist/core_feature_table_relative_mean.csv",
  tree_file = "../../../metadata/tree_metadata_merge_info_final_align_tree.nwk",
  metadata_file = "../../../metadata/tree_metadata_merge_info.csv",
  color_file = "../../../feature_table/Tree_top_order_color.csv",
  output_prefix = "plant_order_level_coevolved_protistan_asv_loop"
)
