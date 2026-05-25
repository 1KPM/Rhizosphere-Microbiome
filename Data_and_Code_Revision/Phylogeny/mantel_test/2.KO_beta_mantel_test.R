### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
library(ggtree)
library(ggplot2)
library(reshape2)
library(agricolae)
library(picante)
library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(gridExtra)
library(ggplotify)
library(geosphere)

### Import data ----------------------------------------------------------------
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

### Keep sample with tree_fun_KO_reads
metadata_rs <- metadata_rs[metadata_rs$TreeID %in% tree_file$tip.label, ]
match(tree_file$tip.label, unique(metadata_rs$TreeID)) ###说明有的树木样本没有扩增子测序数据

tree_dist <- as.data.frame(cophenetic.phylo(tree_file))

####################################################################################### 
###(1) bacteria KO
### Import and filter data
###
bray_bac_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Bacteria.quantitative.bray.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(bray_bac_KO) <- gsub("RS.", "", names(bray_bac_KO))
row.names(bray_bac_KO) <- gsub("RS-", "", row.names(bray_bac_KO))
bray_bac_KO <- bray_bac_KO[row.names(bray_bac_KO) %in% metadata_tree$TreeID, names(bray_bac_KO) %in% metadata_tree$TreeID]

###
jaccard_bac_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Bacteria.quantitative.jaccard.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(jaccard_bac_KO) <- gsub("RS.", "", names(jaccard_bac_KO))
row.names(jaccard_bac_KO) <- gsub("RS-", "", row.names(jaccard_bac_KO))
jaccard_bac_KO <- jaccard_bac_KO[row.names(jaccard_bac_KO) %in% metadata_tree$TreeID, names(jaccard_bac_KO) %in% metadata_tree$TreeID]

#######################################################################################
###(2) fungi KO
### Import and filter data
###
bray_fun_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Fungi.quantitative.bray.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(bray_fun_KO) <- gsub("RS.", "", names(bray_fun_KO))
row.names(bray_fun_KO) <- gsub("RS-", "", row.names(bray_fun_KO))
bray_fun_KO <- bray_fun_KO[row.names(bray_fun_KO) %in% metadata_tree$TreeID, names(bray_fun_KO) %in% metadata_tree$TreeID]

###
jaccard_fun_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Fungi.quantitative.jaccard.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(jaccard_fun_KO) <- gsub("RS.", "", names(jaccard_fun_KO))
row.names(jaccard_fun_KO) <- gsub("RS-", "", row.names(jaccard_fun_KO))
jaccard_fun_KO <- jaccard_fun_KO[row.names(jaccard_fun_KO) %in% metadata_tree$TreeID, names(jaccard_fun_KO) %in% metadata_tree$TreeID]

#######################################################################################
###(3) protist KO
### Import and filter data
###
bray_pro_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Protist.quantitative.bray.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(bray_pro_KO) <- gsub("RS.", "", names(bray_pro_KO))
row.names(bray_pro_KO) <- gsub("RS-", "", row.names(bray_pro_KO))
bray_pro_KO <- bray_pro_KO[row.names(bray_pro_KO) %in% metadata_tree$TreeID, names(bray_pro_KO) %in% metadata_tree$TreeID]

###
jaccard_pro_KO <- read.csv(file = "../ko_qa_alpha_beta/diversity/Protist.quantitative.jaccard.csv", fileEncoding = "GBK", row.names = 1, header = T, sep = ",")
names(jaccard_pro_KO) <- gsub("RS.", "", names(jaccard_pro_KO))
row.names(jaccard_pro_KO) <- gsub("RS-", "", row.names(jaccard_pro_KO))
jaccard_pro_KO <- jaccard_pro_KO[row.names(jaccard_pro_KO) %in% metadata_tree$TreeID, names(jaccard_pro_KO) %in% metadata_tree$TreeID]

#######################################################################################
###（4）Mantel tests
### Mantel tests between community similarity and phylogenetic relatedness 
### need to rearrange community similarity matrix to match phylogeny 

### 从metadata_tree中提取经纬度信息
# 假设metadata_tree数据框包含TreeID、Longitude、Latitude列
# 首先确保每个TreeID只有一个经纬度（如果有重复，取平均值）
tree_coords <- aggregate(cbind(Longitude, Latitude) ~ TreeID, 
                         data = metadata_tree, 
                         FUN = mean)

# 设置行名
row.names(tree_coords) <- tree_coords$TreeID

# 创建地理距离矩阵
geo_dist_all <- distm(tree_coords[, c("Longitude", "Latitude")], 
                      fun = distGeo) / 1000  # 转换为公里

# 设置行名和列名
row.names(geo_dist_all) <- tree_coords$TreeID
colnames(geo_dist_all) <- tree_coords$TreeID

### 定义函数进行Mantel检验
perform_mantel_tests <- function(distance_matrices_list, tree_dist_mat, group_name) {
  ### 准备存储结果的列表
  mantel_results <- list()
  
  ### 循环进行Mantel检验
  for (i in seq_along(distance_matrices_list)) {
    matrix_name <- names(distance_matrices_list)[i]
    cat(group_name, ":", matrix_name, "\n")
    
    # 获取当前距离矩阵
    dist_mat <- distance_matrices_list[[i]]
    
    # 确保距离矩阵与树距离矩阵的顺序一致
    dist_mat <- dist_mat[match(names(tree_dist_mat), names(dist_mat)), 
                         match(names(tree_dist_mat), names(dist_mat))]
    
    # 执行Mantel检验
    mantel_result <- mantel(as.dist(tree_dist_mat), as.dist(dist_mat), 
                            method = "spearman", permutations = 999)
    
    # 保存结果
    mantel_results[[matrix_name]] <- mantel_result
    
    # 打印结果
    cat("Mantel's r:", mantel_result$statistic, "\n")
    cat("P-value:", mantel_result$signif, "\n\n")
  }
  
  return(mantel_results)
}

#######################################################################################
###（5）partial Mantel tests
### 定义函数进行partial Mantel检验（控制地理距离）
perform_mantel_partial_tests <- function(distance_matrices_list, tree_dist_mat, geo_dist_mat, group_name) {
  ### 准备存储结果的列表
  mantel_results <- list()
  
  ### 循环进行偏Mantel检验
  for (i in seq_along(distance_matrices_list)) {
    matrix_name <- names(distance_matrices_list)[i]
    cat(group_name, ":", matrix_name, "(控制地理距离)", "\n")
    
    # 获取当前距离矩阵
    dist_mat <- distance_matrices_list[[i]]
    
    # 确保所有距离矩阵的顺序一致
    common_ids <- Reduce(intersect, list(
      names(tree_dist_mat), 
      names(geo_dist_mat), 
      names(dist_mat)
    ))
    
    # 重新排序所有矩阵
    tree_mat <- tree_dist_mat[common_ids, common_ids]
    geo_mat <- geo_dist_mat[common_ids, common_ids]
    microbe_mat <- dist_mat[common_ids, common_ids]
    
    # 执行偏Mantel检验
    mantel_result <- mantel.partial(as.dist(tree_mat), as.dist(microbe_mat), as.dist(geo_mat),
                                    method = "spearman", permutations = 999)
    
    # 保存结果
    mantel_results[[matrix_name]] <- mantel_result
    
    # 打印结果
    cat("Partial Mantel's r:", mantel_result$statistic, "\n")
    cat("P-value:", mantel_result$signif, "\n\n")
  }
  
  return(mantel_results)
}

#######################################################################################
### (6) 对每个组执行Mantel检验和偏Mantel检验（控制地理距离）
# Bacteria KO
bacteria_distances <- list(
  bray = bray_bac_KO,
  jaccard = jaccard_bac_KO
)

# 准备地理距离矩阵子集
geo_dist_bacteria <- geo_dist_all[row.names(geo_dist_all) %in% names(bray_bac_KO), 
                                  colnames(geo_dist_all) %in% names(bray_bac_KO)]
geo_dist_bacteria <- as.data.frame(geo_dist_bacteria)

tree_dist_bacteria <- tree_dist[row.names(tree_dist) %in% names(bray_bac_KO), 
                                names(tree_dist) %in% names(bray_bac_KO)]

# 普通Mantel检验
bacteria_results <- perform_mantel_tests(bacteria_distances, tree_dist_bacteria, "Bacteria KO")

# 偏Mantel检验（控制地理距离）
bacteria_results_partial <- perform_mantel_partial_tests(bacteria_distances, tree_dist_bacteria, 
                                                         geo_dist_bacteria, "Bacteria KO")

# Fungi KO
fungi_distances <- list(
  bray = bray_fun_KO,
  jaccard = jaccard_fun_KO
)

# 准备地理距离矩阵子集
geo_dist_fungi <- geo_dist_all[row.names(geo_dist_all) %in% names(bray_fun_KO), 
                               colnames(geo_dist_all) %in% names(bray_fun_KO)]

geo_dist_fungi <- as.data.frame(geo_dist_fungi)

tree_dist_fungi <- tree_dist[row.names(tree_dist) %in% names(bray_fun_KO), 
                             names(tree_dist) %in% names(bray_fun_KO)]

# 普通Mantel检验
fungi_results <- perform_mantel_tests(fungi_distances, tree_dist_fungi, "Fungi KO")

# 偏Mantel检验（控制地理距离）
fungi_results_partial <- perform_mantel_partial_tests(fungi_distances, tree_dist_fungi, 
                                                      geo_dist_fungi, "Fungi KO")

# Protist KO
protist_distances <- list(
  bray = bray_pro_KO,
  jaccard = jaccard_pro_KO
)

# 准备地理距离矩阵子集
geo_dist_protist <- geo_dist_all[row.names(geo_dist_all) %in% names(bray_pro_KO), 
                                 colnames(geo_dist_all) %in% names(bray_pro_KO)]

geo_dist_protist <- as.data.frame(geo_dist_protist)

tree_dist_protist <- tree_dist[row.names(tree_dist) %in% names(bray_pro_KO), 
                               names(tree_dist) %in% names(bray_pro_KO)]

# 普通Mantel检验
protist_results <- perform_mantel_tests(protist_distances, tree_dist_protist, "Protist KO")

# 偏Mantel检验（控制地理距离）
protist_results_partial <- perform_mantel_partial_tests(protist_distances, tree_dist_protist, 
                                                        geo_dist_protist, "Protist KO")

### 整理所有普通Mantel检验结果
beta_all_mantel <- as.data.frame(rbind(
  c("Bacteria", "Bray-Curtis", bacteria_results$bray$statistic, bacteria_results$bray$signif),
  c("Bacteria", "Jaccard", bacteria_results$jaccard$statistic, bacteria_results$jaccard$signif),
  c("Fungi", "Bray-Curtis", fungi_results$bray$statistic, fungi_results$bray$signif),
  c("Fungi", "Jaccard", fungi_results$jaccard$statistic, fungi_results$jaccard$signif),
  c("Protist", "Bray-Curtis", protist_results$bray$statistic, protist_results$bray$signif),
  c("Protist", "Jaccard", protist_results$jaccard$statistic, protist_results$jaccard$signif)
))

names(beta_all_mantel) <- c("Group", "Distance_Metric", "Mantel's r", "P")

### 整理所有偏Mantel检验结果
beta_all_mantel_partial <- as.data.frame(rbind(
  c("Bacteria", "Bray-Curtis", bacteria_results_partial$bray$statistic, bacteria_results_partial$bray$signif),
  c("Bacteria", "Jaccard", bacteria_results_partial$jaccard$statistic, bacteria_results_partial$jaccard$signif),
  c("Fungi", "Bray-Curtis", fungi_results_partial$bray$statistic, fungi_results_partial$bray$signif),
  c("Fungi", "Jaccard", fungi_results_partial$jaccard$statistic, fungi_results_partial$jaccard$signif),
  c("Protist", "Bray-Curtis", protist_results_partial$bray$statistic, protist_results_partial$bray$signif),
  c("Protist", "Jaccard", protist_results_partial$jaccard$statistic, protist_results_partial$jaccard$signif)
))

names(beta_all_mantel_partial) <- c("Group", "Distance_Metric", "Partial_Mantel_r", "P")

### 合并两个结果表
beta_all_mantel_combined <- merge(beta_all_mantel, beta_all_mantel_partial, 
                                  by = c("Group", "Distance_Metric"), 
                                  all = TRUE)

# 重命名列
names(beta_all_mantel_combined) <- c("Group", "Distance_Metric", 
                                     "Mantel_r", "Mantel_P", 
                                     "Partial_Mantel_r", "Partial_Mantel_P")

# 显示结果
print("普通Mantel检验结果：")
print(beta_all_mantel_combined[, c("Group", "Distance_Metric", "Mantel_r", "Mantel_P")])

print("\n偏Mantel检验结果（控制地理距离）：")
print(beta_all_mantel_combined[, c("Group", "Distance_Metric", "Partial_Mantel_r", "Partial_Mantel_P")])

###结果表明控制地点并没有对结果产生多大的影响

### 保存结果到文件
write.csv(beta_all_mantel_combined, file = "./meta_stats/ko_all_groups_mantel_and_geo_controlled_partial_mantel_test.csv",
          fileEncoding = "GBK", row.names = FALSE)

#######################################################################################
###(7) 按照植物目（Order）进行分组，并对每个植物目分别进行细菌、真菌、原生生物的Mantel和偏Mantel检验

# 读取top10植物目信息
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]

# 从metadata_tree中获取每个植物目对应的TreeID
order_trees <- split(metadata_tree$TreeID, metadata_tree$Order)
order_trees <- order_trees[names(order_trees) %in% top10orders]

# 对每个植物目进行分析
all_order_results <- list()

for (order_name in names(order_trees)) {
  # 获取该植物目的TreeID
  order_tree_ids <- order_trees[[order_name]]
  
  # 提取该植物目的距离矩阵子集
  order_mats <- list()
  
  # 为每个距离矩阵提取子集
  all_matrices <- list(
    bray_bac_KO = bray_bac_KO,
    jaccard_bac_KO = jaccard_bac_KO,
    bray_fun_KO = bray_fun_KO,
    jaccard_fun_KO = jaccard_fun_KO,
    bray_pro_KO = bray_pro_KO,
    jaccard_pro_KO = jaccard_pro_KO
  )
  
  for (mat_name in names(all_matrices)) {
    mat <- all_matrices[[mat_name]]
    common_ids <- intersect(order_tree_ids, row.names(mat))
    order_mats[[mat_name]] <- mat[common_ids, common_ids]
  }
  
  # 进行Mantel检验
  cat("\n植物目:", order_name, "\n")
  
  # 存储该植物目的结果
  order_results <- list()
  
  # 微生物组定义
  microbe_groups <- list(
    Bacteria = c("bray_bac_KO", "jaccard_bac_KO"),
    Fungi = c("bray_fun_KO", "jaccard_fun_KO"),
    Protist = c("bray_pro_KO", "jaccard_pro_KO")
  )
  
  for (group_name in names(microbe_groups)) {
    cat("  ", group_name, ":\n")
    
    for (mat_name in microbe_groups[[group_name]]) {
      microbe_mat <- order_mats[[mat_name]]
      
      # 获取共有的TreeID
      common_ids <- Reduce(intersect, list(
        row.names(microbe_mat), 
        row.names(tree_dist), 
        row.names(geo_dist_all)
      ))
      
      tree_subset <- tree_dist[common_ids, common_ids]
      geo_subset <- geo_dist_all[common_ids, common_ids]
      microbe_subset <- microbe_mat[common_ids, common_ids]
      
      # 普通Mantel检验
      mantel_result <- mantel(as.dist(tree_subset), as.dist(microbe_subset), 
                              method = "spearman", permutations = 999)
      
      # 偏Mantel检验（控制地理距离）
      mantel_partial_result <- mantel.partial(as.dist(tree_subset), as.dist(microbe_subset), 
                                              as.dist(geo_subset), 
                                              method = "spearman", permutations = 999)
      
      # 提取距离度量名称
      dist_metric <- gsub("_\\d+S$|_fun_KO$", "", mat_name)
      dist_metric <- gsub("_", " ", dist_metric)
      dist_metric <- tools::toTitleCase(dist_metric)
      
      # 保存结果
      result_row <- data.frame(
        Order = order_name,
        Microbe_Group = group_name,
        Distance_Metric = dist_metric,
        Sample_Size = length(common_ids),
        Mantel_r = mantel_result$statistic,
        Mantel_P = mantel_result$signif,
        Partial_Mantel_r = mantel_partial_result$statistic,
        Partial_Mantel_P = mantel_partial_result$signif,
        stringsAsFactors = FALSE
      )
      
      order_results[[paste(group_name, mat_name)]] <- result_row
      
      cat("    ", dist_metric, ": ", 
          "普通: r =", round(mantel_result$statistic, 3), "P =", round(mantel_result$signif, 4), 
          " | 偏: r =", round(mantel_partial_result$statistic, 3), "P =", round(mantel_partial_result$signif, 4), 
          "\n")
    }
  }
  
  # 合并该植物目的结果
  order_df <- do.call(rbind, order_results)
  all_order_results[[order_name]] <- order_df
}

# 整理所有结果
combined_results <- do.call(rbind, all_order_results)
rownames(combined_results) <- NULL

# 保存结果
write.csv(combined_results, file = "./meta_stats/ko_mantel_results_by_plant_order_with_partial.csv",
          fileEncoding = "GBK", row.names = FALSE)

cat("\n分析完成！结果已保存。\n")


#######################################################################################
###(8) 按照植物纲（Class）进行分组，并对每个植物纲分别进行细菌、真菌、原生生物的Mantel和偏Mantel检验

# 指定要分析的植物纲
target_classes <- c("Magnoliopsida", "Liliopsida")

# 从metadata_tree中获取每个植物纲对应的TreeID
class_trees <- split(metadata_tree$TreeID, metadata_tree$Class)
class_trees <- class_trees[names(class_trees) %in% target_classes]

# 对每个植物纲进行分析
all_class_results <- list()

for (class_name in names(class_trees)) {
  # 获取该植物纲的TreeID
  class_tree_ids <- class_trees[[class_name]]
  
  # 提取该植物纲的距离矩阵子集
  class_mats <- list()
  
  # 为每个距离矩阵提取子集
  all_matrices <- list(
    bray_bac_KO = bray_bac_KO,
    jaccard_bac_KO = jaccard_bac_KO,
    bray_fun_KO = bray_fun_KO,
    jaccard_fun_KO = jaccard_fun_KO,
    bray_pro_KO = bray_pro_KO,
    jaccard_pro_KO = jaccard_pro_KO
  )
  
  for (mat_name in names(all_matrices)) {
    mat <- all_matrices[[mat_name]]
    common_ids <- intersect(class_tree_ids, row.names(mat))
    class_mats[[mat_name]] <- mat[common_ids, common_ids]
  }
  
  # 进行Mantel检验
  cat("\n植物纲:", class_name, "\n")
  
  # 存储该植物纲的结果
  class_results <- list()
  
  # 微生物组定义
  microbe_groups <- list(
    Bacteria = c("bray_bac_KO", "jaccard_bac_KO"),
    Fungi = c("bray_fun_KO", "jaccard_fun_KO"),
    Protist = c("bray_pro_KO", "jaccard_pro_KO")
  )
  
  for (group_name in names(microbe_groups)) {
    cat("  ", group_name, ":\n")
    
    for (mat_name in microbe_groups[[group_name]]) {
      microbe_mat <- class_mats[[mat_name]]
      
      # 获取共有的TreeID
      common_ids <- Reduce(intersect, list(
        row.names(microbe_mat), 
        row.names(tree_dist), 
        row.names(geo_dist_all)
      ))
      
      tree_subset <- tree_dist[common_ids, common_ids]
      geo_subset <- geo_dist_all[common_ids, common_ids]
      microbe_subset <- microbe_mat[common_ids, common_ids]
      
      # 普通Mantel检验
      mantel_result <- mantel(as.dist(tree_subset), as.dist(microbe_subset), 
                              method = "spearman", permutations = 999)
      
      # 偏Mantel检验（控制地理距离）
      mantel_partial_result <- mantel.partial(as.dist(tree_subset), as.dist(microbe_subset), 
                                              as.dist(geo_subset), 
                                              method = "spearman", permutations = 999)
      
      # 提取距离度量名称
      dist_metric <- gsub("_\\d+S$|_fun_KO$", "", mat_name) #不影响结果
      dist_metric <- gsub("_", " ", dist_metric)
      dist_metric <- tools::toTitleCase(dist_metric)
      
      # 保存结果
      result_row <- data.frame(
        Class = class_name,
        Microbe_Group = group_name,
        Distance_Metric = dist_metric,
        Sample_Size = length(common_ids),
        Mantel_r = mantel_result$statistic,
        Mantel_P = mantel_result$signif,
        Partial_Mantel_r = mantel_partial_result$statistic,
        Partial_Mantel_P = mantel_partial_result$signif,
        stringsAsFactors = FALSE
      )
      
      class_results[[paste(group_name, mat_name)]] <- result_row
      
      cat("    ", dist_metric, ": ", 
          "普通: r =", round(mantel_result$statistic, 3), "P =", round(mantel_result$signif, 4), 
          " | 偏: r =", round(mantel_partial_result$statistic, 3), "P =", round(mantel_partial_result$signif, 4), 
          "\n")
    }
  }
  
  # 合并该植物纲的结果
  class_df <- do.call(rbind, class_results)
  all_class_results[[class_name]] <- class_df
}

# 整理所有结果
combined_results <- do.call(rbind, all_class_results)
rownames(combined_results) <- NULL

# 保存结果
write.csv(combined_results, file = "./meta_stats/ko_mantel_results_by_plant_class_with_partial.csv", fileEncoding = "GBK", row.names = FALSE)

cat("\n分析完成！结果已保存。\n")


#######################################################################################
###(9) 可视化总体样本及按照植物目（Order）进行分组，并对每个植物目分别进行细菌、真菌、原生生物的Mantel和偏Mantel检验

###y轴表示不同目，x轴表示不同距离矩阵，不同界微生物分页
###颜色表示P值，气泡大小表示Mantel'R

###
beta_all_mantel <- read.csv(file = "./meta_stats/ko_all_groups_mantel_and_geo_controlled_partial_mantel_test.csv", fileEncoding = "GBK")
beta_class_mantel <- read.csv(file = "./meta_stats/ko_mantel_results_by_plant_class_with_partial.csv", fileEncoding = "GBK")
beta_orders_mantel <- read.csv(file = "./meta_stats/ko_mantel_results_by_plant_order_with_partial.csv", fileEncoding = "GBK")

beta_all_mantel$Plant_Group <- "All"
names(beta_all_mantel)[c(1)] <- c("Microbe_Group")
names(beta_class_mantel)[1] <- c("Plant_Group")
names(beta_orders_mantel)[c(1)] <- c("Plant_Group")

beta_all_mantel <- as.data.frame(rbind(beta_all_mantel[,c(7,1:6)], beta_class_mantel[,-4], beta_orders_mantel[,-4]))
beta_all_mantel$Distance_Metric <- gsub("^Bray.*", "Bray-Curtis", beta_all_mantel$Distance_Metric)
beta_all_mantel$Distance_Metric <- gsub("^Jaccard.*", "Jaccard", beta_all_mantel$Distance_Metric)
beta_all_mantel$Microbe_Group <- gsub("^Protist$", "Protists", beta_all_mantel$Microbe_Group)


beta_all_mantel$Distance_Metric <- factor(beta_all_mantel$Distance_Metric, levels = c("Jaccard", "Bray-Curtis"))
beta_all_mantel$Plant_Group <- factor(beta_all_mantel$Plant_Group, levels = c("All", target_classes, top10orders))

beta_all_mantel$Signficant <- ifelse(beta_all_mantel$Mantel_P < 0.01, "**", 
                                     ifelse(beta_all_mantel$Mantel_P < 0.05, "*", ""))

###
mantel_plot <- ggplot(beta_all_mantel, aes(x=Distance_Metric, y=Plant_Group)) +
  geom_point(aes(size= Mantel_r, color=-log2(Mantel_P))) +
  geom_text(aes(label = Signficant, vjust = 0.7)) +
  scale_colour_gradient2(name="-log2(Mantel's P)",low="darkgreen", mid = "white", high="red") +
  scale_size_continuous(name="Mantel's R") +
  labs(title="", x="", y="") +
  facet_grid(~ Microbe_Group) +
  theme(legend.position="right") + 
  theme(legend.title = element_text(size=8,color="black"),
        legend.text= element_text(size=8,color="black"),
        axis.text.x = element_text(size=8,color="black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=8,color="black"),
        axis.title= element_text(size=8,color="black")) 

mantel_plot

name <- "./meta_plots/Mantel_tests_between_microbial_functional_similarity_and_plant_phylogenetic_relatedness"
width <- 10
height <- 9.2
ggsave(paste0(name, ".png"), mantel_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), mantel_plot, width = width, height = height, units = "cm")

##################################################################
###partial mantel test
beta_all_mantel$Signficant_partial <- ifelse(beta_all_mantel$Partial_Mantel_P < 0.01, "**", 
                                             ifelse(beta_all_mantel$Partial_Mantel_P < 0.05, "*", ""))


###
partial_mantel_plot <- ggplot(beta_all_mantel, aes(x=Distance_Metric, y=Plant_Group)) +
  geom_point(aes(size= Partial_Mantel_r, color=-log2(Partial_Mantel_P))) +
  geom_text(aes(label = Signficant_partial, vjust = 0.7)) +
  scale_colour_gradient2(name="-log2(Partial Mantel's P)",low="darkgreen", mid = "white", high="red") +
  scale_size_continuous(name="Partial Mantel's R") +
  labs(title="", x="", y="") +
  facet_grid(~ Microbe_Group) +
  theme(legend.position="right") + 
  theme(legend.title = element_text(size=8,color="black"),
        legend.text= element_text(size=8,color="black"),
        axis.text.x = element_text(size=8,color="black", angle = 45, hjust = 1),
        axis.text.y = element_text(size=8,color="black"),
        axis.title= element_text(size=8,color="black")) 

partial_mantel_plot

name <- "./meta_plots/Partial_Mantel_tests_between_microbial_functional_similarity_and_plant_phylogenetic_relatedness"
width <- 10
height <- 9.2
ggsave(paste0(name, ".png"), partial_mantel_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), partial_mantel_plot, width = width, height = height, units = "cm")

write.csv(beta_all_mantel, file = "./meta_stats/ko_mantel_and_partial_mantel_results_all_and_by_plant_class_order.csv", fileEncoding = "GBK", row.names = FALSE)