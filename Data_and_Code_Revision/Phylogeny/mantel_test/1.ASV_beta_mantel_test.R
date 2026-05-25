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

### Import data ----------------------------------------------------------------
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

### Keep sample with tree_ITS_reads
metadata_rs <- metadata_rs[metadata_rs$TreeID %in% tree_file$tip.label, ]
match(tree_file$tip.label, unique(metadata_rs$TreeID)) ###说明有的树木样本没有扩增子测序数据

tree_dist <- as.data.frame(cophenetic.phylo(tree_file))

####################################################################################### 
###(1) bacteria 16S
### Import and filter data
### average beta diversity index based on tree ID
bray_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__bray_curtis.csv", fileEncoding = "GBK", row.names = 1)
bray_16S <- bray_16S[row.names(bray_16S) %in% metadata_rs$FileID, names(bray_16S) %in% metadata_rs$FileID]

jaccard_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__jaccard.csv", fileEncoding = "GBK", row.names = 1)
jaccard_16S <- jaccard_16S[row.names(jaccard_16S) %in% metadata_rs$FileID, names(jaccard_16S) %in% metadata_rs$FileID]

weighted_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__weighted.csv", fileEncoding = "GBK", row.names = 1)
weighted_16S <- weighted_16S[row.names(weighted_16S) %in% metadata_rs$FileID, names(weighted_16S) %in% metadata_rs$FileID]

unweighted_16S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/16S_beta_diversity_absolute__unweighted.csv", fileEncoding = "GBK", row.names = 1)
unweighted_16S <- unweighted_16S[row.names(unweighted_16S) %in% metadata_rs$FileID, names(unweighted_16S) %in% metadata_rs$FileID]

# 定义矩阵求平均处理函数
calculate_mean_distance <- function(dist_matrix, metadata) {
  dist_matrix$ID <- row.names(dist_matrix)
  dist_matrix_double <- melt(dist_matrix)
  dist_matrix_double <- merge(dist_matrix_double, unique(metadata[,1:2]), 
                              by.x = "ID", by.y = "FileID")
  dist_matrix_double <- merge(dist_matrix_double, unique(metadata[,1:2]), 
                              by.x = "variable", by.y = "FileID")
  dist_matrix_mean <- aggregate(value ~ TreeID.x + TreeID.y, 
                                data = dist_matrix_double, FUN = mean)
  dist_matrix_mean <- dcast(dist_matrix_mean, TreeID.x ~ TreeID.y)
  row.names(dist_matrix_mean) <- dist_matrix_mean$TreeID.x
  dist_matrix_mean <- dist_matrix_mean[,-1]
  return(dist_matrix_mean)
}

# 对bray_16S取平均值
bray_16S_mean <- calculate_mean_distance(bray_16S, metadata_rs)
# write.csv(bray_16S_mean, file = "../asv_qa_alpha_beta/rhizo_16S_bray_curtis_tree_mean.csv", fileEncoding = "GBK")

# 对jaccard_16S取平均值
jaccard_16S_mean <- calculate_mean_distance(jaccard_16S, metadata_rs)
# write.csv(jaccard_16S_mean, file = "../asv_qa_alpha_beta/rhizo_16S_jaccard_tree_mean.csv", fileEncoding = "GBK")

# 对weighted_16S取平均值
weighted_16S_mean <- calculate_mean_distance(weighted_16S, metadata_rs)
# write.csv(weighted_16S_mean, file = "../asv_qa_alpha_beta/rhizo_16S_weighted_unifrac_tree_mean.csv", fileEncoding = "GBK")

# 对unweighted_16S取平均值
unweighted_16S_mean <- calculate_mean_distance(unweighted_16S, metadata_rs)
# write.csv(unweighted_16S_mean, file = "../asv_qa_alpha_beta/rhizo_16S_unweighted_unifrac_tree_mean.csv", fileEncoding = "GBK")

#######################################################################################
###(2) fungi ITS
### Import and filter data
### average beta diversity index based on tree ID
bray_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__bray_curtis.csv", fileEncoding = "GBK", row.names = 1)
bray_ITS <- bray_ITS[row.names(bray_ITS) %in% metadata_rs$FileID, names(bray_ITS) %in% metadata_rs$FileID]

jaccard_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__jaccard.csv", fileEncoding = "GBK", row.names = 1)
jaccard_ITS <- jaccard_ITS[row.names(jaccard_ITS) %in% metadata_rs$FileID, names(jaccard_ITS) %in% metadata_rs$FileID]

weighted_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__weighted.csv", fileEncoding = "GBK", row.names = 1)
weighted_ITS <- weighted_ITS[row.names(weighted_ITS) %in% metadata_rs$FileID, names(weighted_ITS) %in% metadata_rs$FileID]

unweighted_ITS <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/ITS_beta_diversity_absolute__unweighted.csv", fileEncoding = "GBK", row.names = 1)
unweighted_ITS <- unweighted_ITS[row.names(unweighted_ITS) %in% metadata_rs$FileID, names(unweighted_ITS) %in% metadata_rs$FileID]

# 对bray_ITS取平均值
bray_ITS_mean <- calculate_mean_distance(bray_ITS, metadata_rs)
# write.csv(bray_ITS_mean, file = "../asv_qa_alpha_beta/rhizo_ITS_bray_curtis_tree_mean.csv", fileEncoding = "GBK")

# 对jaccard_ITS取平均值
jaccard_ITS_mean <- calculate_mean_distance(jaccard_ITS, metadata_rs)
# write.csv(jaccard_ITS_mean, file = "../asv_qa_alpha_beta/rhizo_ITS_jaccard_tree_mean.csv", fileEncoding = "GBK")

# 对weighted_ITS取平均值
weighted_ITS_mean <- calculate_mean_distance(weighted_ITS, metadata_rs)
# write.csv(weighted_ITS_mean, file = "../asv_qa_alpha_beta/rhizo_ITS_weighted_unifrac_tree_mean.csv", fileEncoding = "GBK")

# 对unweighted_ITS取平均值
unweighted_ITS_mean <- calculate_mean_distance(unweighted_ITS, metadata_rs)
# write.csv(unweighted_ITS_mean, file = "../asv_qa_alpha_beta/rhizo_ITS_unweighted_unifrac_tree_mean.csv", fileEncoding = "GBK")


#######################################################################################
###(3) protist 18S
### Import and filter data
### average beta diversity index based on tree ID
bray_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__bray_curtis.csv", fileEncoding = "GBK", row.names = 1)
bray_18S <- bray_18S[row.names(bray_18S) %in% metadata_rs$FileID, names(bray_18S) %in% metadata_rs$FileID]

jaccard_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__jaccard.csv", fileEncoding = "GBK", row.names = 1)
jaccard_18S <- jaccard_18S[row.names(jaccard_18S) %in% metadata_rs$FileID, names(jaccard_18S) %in% metadata_rs$FileID]

weighted_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__weighted.csv", fileEncoding = "GBK", row.names = 1)
weighted_18S <- weighted_18S[row.names(weighted_18S) %in% metadata_rs$FileID, names(weighted_18S) %in% metadata_rs$FileID]

unweighted_18S <- read.csv(file = "../asv_qa_alpha_beta/02-get_absolute_beta_diversity/Protist_beta_diversity_absolute__unweighted.csv", fileEncoding = "GBK", row.names = 1)
unweighted_18S <- unweighted_18S[row.names(unweighted_18S) %in% metadata_rs$FileID, names(unweighted_18S) %in% metadata_rs$FileID]

# 对bray_18S取平均值
bray_18S_mean <- calculate_mean_distance(bray_18S, metadata_rs)
# write.csv(bray_18S_mean, file = "../asv_qa_alpha_beta/rhizo_18S_bray_curtis_tree_mean.csv", fileEncoding = "GBK")

# 对jaccard_18S取平均值
jaccard_18S_mean <- calculate_mean_distance(jaccard_18S, metadata_rs)
# write.csv(jaccard_18S_mean, file = "../asv_qa_alpha_beta/rhizo_18S_jaccard_tree_mean.csv", fileEncoding = "GBK")

# 对weighted_18S取平均值
weighted_18S_mean <- calculate_mean_distance(weighted_18S, metadata_rs)
# write.csv(weighted_18S_mean, file = "../asv_qa_alpha_beta/rhizo_18S_weighted_unifrac_tree_mean.csv", fileEncoding = "GBK")

# 对unweighted_18S取平均值
unweighted_18S_mean <- calculate_mean_distance(unweighted_18S, metadata_rs)
# write.csv(unweighted_18S_mean, file = "../asv_qa_alpha_beta/rhizo_18S_unweighted_unifrac_tree_mean.csv", fileEncoding = "GBK")



#######################################################################################
###（4）Mantel tests
### Mantel tests between community similarity and phylogenetic relatedness 
### need to rearrange community similarity matrix to match phylogeny 

# ###Bacteria 16S
# ###subset
# tree_dist_16S <- tree_dist[row.names(tree_dist) %in% names(bray_16S_mean), names(tree_dist) %in% names(bray_16S_mean)]
# ###rearrange
# bray_16S_mean <- bray_16S_mean[match(names(tree_dist_16S), names(bray_16S_mean)), match(names(tree_dist_16S), names(bray_16S_mean))]
# 
# bray_16S_mantel <- mantel(as.dist(tree_dist_16S), as.dist(bray_16S_mean), method="spearman", permutations=999) 
# bray_16S_mantel$statistic
# bray_16S_mantel$signif

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

### 从metadata_tree中提取经纬度信息
# 假设metadata_tree数据框包含TreeID、Longitude、Latitude列
# 首先确保每个TreeID只有一个经纬度（如果有重复，取平均值）
tree_coords <- aggregate(cbind(Longitude, Latitude) ~ TreeID, 
                         data = metadata_tree, 
                         FUN = mean)

# 设置行名
row.names(tree_coords) <- tree_coords$TreeID

# 使用geosphere包计算大圆距离
library(geosphere)

# 创建地理距离矩阵
geo_dist_all <- distm(tree_coords[, c("Longitude", "Latitude")], 
                      fun = distGeo) / 1000  # 转换为公里

# 设置行名和列名
row.names(geo_dist_all) <- tree_coords$TreeID
colnames(geo_dist_all) <- tree_coords$TreeID


#######################################################################################
### (6) 对每个组执行Mantel检验和偏Mantel检验（控制地理距离）
# Bacteria 16S
bacteria_distances <- list(
  bray = bray_16S_mean,
  jaccard = jaccard_16S_mean,
  weighted = weighted_16S_mean,
  unweighted = unweighted_16S_mean
)

# 准备地理距离矩阵子集
geo_dist_bacteria <- geo_dist_all[row.names(geo_dist_all) %in% names(bray_16S_mean), 
                                  colnames(geo_dist_all) %in% names(bray_16S_mean)]
geo_dist_bacteria <- as.data.frame(geo_dist_bacteria)

tree_dist_bacteria <- tree_dist[row.names(tree_dist) %in% names(bray_16S_mean), 
                                names(tree_dist) %in% names(bray_16S_mean)]

# 普通Mantel检验
bacteria_results <- perform_mantel_tests(bacteria_distances, tree_dist_bacteria, "Bacteria 16S")

# 偏Mantel检验（控制地理距离）
bacteria_results_partial <- perform_mantel_partial_tests(bacteria_distances, tree_dist_bacteria, 
                                                         geo_dist_bacteria, "Bacteria 16S")

# Fungi ITS
fungi_distances <- list(
  bray = bray_ITS_mean,
  jaccard = jaccard_ITS_mean,
  weighted = weighted_ITS_mean,
  unweighted = unweighted_ITS_mean
)

# 准备地理距离矩阵子集
geo_dist_fungi <- geo_dist_all[row.names(geo_dist_all) %in% names(bray_ITS_mean), 
                               colnames(geo_dist_all) %in% names(bray_ITS_mean)]

geo_dist_fungi <- as.data.frame(geo_dist_fungi)

tree_dist_fungi <- tree_dist[row.names(tree_dist) %in% names(bray_ITS_mean), 
                             names(tree_dist) %in% names(bray_ITS_mean)]

# 普通Mantel检验
fungi_results <- perform_mantel_tests(fungi_distances, tree_dist_fungi, "Fungi ITS")

# 偏Mantel检验（控制地理距离）
fungi_results_partial <- perform_mantel_partial_tests(fungi_distances, tree_dist_fungi, 
                                                      geo_dist_fungi, "Fungi ITS")

# Protist 18S
protist_distances <- list(
  bray = bray_18S_mean,
  jaccard = jaccard_18S_mean,
  weighted = weighted_18S_mean,
  unweighted = unweighted_18S_mean
)

# 准备地理距离矩阵子集
geo_dist_protist <- geo_dist_all[row.names(geo_dist_all) %in% names(bray_18S_mean), 
                                 colnames(geo_dist_all) %in% names(bray_18S_mean)]

geo_dist_protist <- as.data.frame(geo_dist_protist)

tree_dist_protist <- tree_dist[row.names(tree_dist) %in% names(bray_18S_mean), 
                               names(tree_dist) %in% names(bray_18S_mean)]

# 普通Mantel检验
protist_results <- perform_mantel_tests(protist_distances, tree_dist_protist, "Protist 18S")

# 偏Mantel检验（控制地理距离）
protist_results_partial <- perform_mantel_partial_tests(protist_distances, tree_dist_protist, 
                                                        geo_dist_protist, "Protist 18S")

### 整理所有普通Mantel检验结果
beta_all_mantel <- as.data.frame(rbind(
  c("Bacteria", "Bray-Curtis", bacteria_results$bray$statistic, bacteria_results$bray$signif),
  c("Bacteria", "Jaccard", bacteria_results$jaccard$statistic, bacteria_results$jaccard$signif),
  c("Bacteria", "Weighted UniFrac", bacteria_results$weighted$statistic, bacteria_results$weighted$signif),
  c("Bacteria", "Unweighted UniFrac", bacteria_results$unweighted$statistic, bacteria_results$unweighted$signif),
  c("Fungi", "Bray-Curtis", fungi_results$bray$statistic, fungi_results$bray$signif),
  c("Fungi", "Jaccard", fungi_results$jaccard$statistic, fungi_results$jaccard$signif),
  c("Fungi", "Weighted UniFrac", fungi_results$weighted$statistic, fungi_results$weighted$signif),
  c("Fungi", "Unweighted UniFrac", fungi_results$unweighted$statistic, fungi_results$unweighted$signif),
  c("Protist", "Bray-Curtis", protist_results$bray$statistic, protist_results$bray$signif),
  c("Protist", "Jaccard", protist_results$jaccard$statistic, protist_results$jaccard$signif),
  c("Protist", "Weighted UniFrac", protist_results$weighted$statistic, protist_results$weighted$signif),
  c("Protist", "Unweighted UniFrac", protist_results$unweighted$statistic, protist_results$unweighted$signif)
))

names(beta_all_mantel) <- c("Group", "Distance_Metric", "Mantel's r", "P")

### 整理所有偏Mantel检验结果
beta_all_mantel_partial <- as.data.frame(rbind(
  c("Bacteria", "Bray-Curtis", bacteria_results_partial$bray$statistic, bacteria_results_partial$bray$signif),
  c("Bacteria", "Jaccard", bacteria_results_partial$jaccard$statistic, bacteria_results_partial$jaccard$signif),
  c("Bacteria", "Weighted UniFrac", bacteria_results_partial$weighted$statistic, bacteria_results_partial$weighted$signif),
  c("Bacteria", "Unweighted UniFrac", bacteria_results_partial$unweighted$statistic, bacteria_results_partial$unweighted$signif),
  c("Fungi", "Bray-Curtis", fungi_results_partial$bray$statistic, fungi_results_partial$bray$signif),
  c("Fungi", "Jaccard", fungi_results_partial$jaccard$statistic, fungi_results_partial$jaccard$signif),
  c("Fungi", "Weighted UniFrac", fungi_results_partial$weighted$statistic, fungi_results_partial$weighted$signif),
  c("Fungi", "Unweighted UniFrac", fungi_results_partial$unweighted$statistic, fungi_results_partial$unweighted$signif),
  c("Protist", "Bray-Curtis", protist_results_partial$bray$statistic, protist_results_partial$bray$signif),
  c("Protist", "Jaccard", protist_results_partial$jaccard$statistic, protist_results_partial$jaccard$signif),
  c("Protist", "Weighted UniFrac", protist_results_partial$weighted$statistic, protist_results_partial$weighted$signif),
  c("Protist", "Unweighted UniFrac", protist_results_partial$unweighted$statistic, protist_results_partial$unweighted$signif)
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
write.csv(beta_all_mantel_combined, file = "./stats/ASV_all_groups_mantel_and_geo_controlled_partial_mantel_test.csv",
          fileEncoding = "GBK", row.names = FALSE)

# 保存地理距离矩阵
write.csv(geo_dist_all, file = "../asv_qa_alpha_beta/geographic_distance_matrix.csv", fileEncoding = "GBK")
saveRDS(geo_dist_all, file = "../asv_qa_alpha_beta/geographic_distance_matrix.rds")

#######################################################################################
###(7) 按照植物目（Order）进行分组，并对每个植物目分别进行细菌、真菌、原生生物的Mantel和偏Mantel检验

# 读取top10植物目信息
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]

# 从metadata_tree中获取每个植物目对应的TreeID
order_trees <- split(metadata_tree$TreeID, metadata_tree$Order)
order_trees <- order_trees[names(order_trees) %in% top10orders]

# 创建地理距离矩阵
library(geosphere)
geo_dist_all <- distm(metadata_tree[, c("Longitude", "Latitude")], 
                      fun = distGeo) / 1000
row.names(geo_dist_all) <- metadata_tree$TreeID
colnames(geo_dist_all) <- metadata_tree$TreeID

# 对每个植物目进行分析
all_order_results <- list()

for (order_name in names(order_trees)) {
  # 获取该植物目的TreeID
  order_tree_ids <- order_trees[[order_name]]
  
  # 提取该植物目的距离矩阵子集
  order_mats <- list()
  
  # 为每个距离矩阵提取子集
  all_matrices <- list(
    bray_16S = bray_16S_mean,
    jaccard_16S = jaccard_16S_mean,
    weighted_16S = weighted_16S_mean,
    unweighted_16S = unweighted_16S_mean,
    bray_ITS = bray_ITS_mean,
    jaccard_ITS = jaccard_ITS_mean,
    weighted_ITS = weighted_ITS_mean,
    unweighted_ITS = unweighted_ITS_mean,
    bray_18S = bray_18S_mean,
    jaccard_18S = jaccard_18S_mean,
    weighted_18S = weighted_18S_mean,
    unweighted_18S = unweighted_18S_mean
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
    Bacteria = c("bray_16S", "jaccard_16S", "weighted_16S", "unweighted_16S"),
    Fungi = c("bray_ITS", "jaccard_ITS", "weighted_ITS", "unweighted_ITS"),
    Protist = c("bray_18S", "jaccard_18S", "weighted_18S", "unweighted_18S")
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
      dist_metric <- gsub("_\\d+S$|_ITS$", "", mat_name)
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
write.csv(combined_results, file = "./stats/ASV_mantel_results_by_plant_order_with_partial.csv",
          fileEncoding = "GBK", row.names = FALSE)

cat("\n分析完成！结果已保存。\n")


#######################################################################################
###(8) 按照植物纲（Class）进行分组，并对每个植物纲分别进行细菌、真菌、原生生物的Mantel和偏Mantel检验

# 指定要分析的植物纲
target_classes <- c("Magnoliopsida", "Liliopsida")

# 从metadata_tree中获取每个植物纲对应的TreeID
class_trees <- split(metadata_tree$TreeID, metadata_tree$Class)
class_trees <- class_trees[names(class_trees) %in% target_classes]

# 创建地理距离矩阵
library(geosphere)
geo_dist_all <- distm(metadata_tree[, c("Longitude", "Latitude")], 
                      fun = distGeo) / 1000
row.names(geo_dist_all) <- metadata_tree$TreeID
colnames(geo_dist_all) <- metadata_tree$TreeID

# 对每个植物纲进行分析
all_class_results <- list()

for (class_name in names(class_trees)) {
  # 获取该植物纲的TreeID
  class_tree_ids <- class_trees[[class_name]]
  
  # 提取该植物纲的距离矩阵子集
  class_mats <- list()
  
  # 为每个距离矩阵提取子集
  all_matrices <- list(
    bray_16S = bray_16S_mean,
    jaccard_16S = jaccard_16S_mean,
    weighted_16S = weighted_16S_mean,
    unweighted_16S = unweighted_16S_mean,
    bray_ITS = bray_ITS_mean,
    jaccard_ITS = jaccard_ITS_mean,
    weighted_ITS = weighted_ITS_mean,
    unweighted_ITS = unweighted_ITS_mean,
    bray_18S = bray_18S_mean,
    jaccard_18S = jaccard_18S_mean,
    weighted_18S = weighted_18S_mean,
    unweighted_18S = unweighted_18S_mean
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
    Bacteria = c("bray_16S", "jaccard_16S", "weighted_16S", "unweighted_16S"),
    Fungi = c("bray_ITS", "jaccard_ITS", "weighted_ITS", "unweighted_ITS"),
    Protist = c("bray_18S", "jaccard_18S", "weighted_18S", "unweighted_18S")
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
      dist_metric <- gsub("_\\d+S$|_ITS$", "", mat_name)
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
write.csv(combined_results, file = "./stats/ASV_mantel_results_by_plant_class_with_partial.csv",
          fileEncoding = "GBK", row.names = FALSE)

cat("\n分析完成！结果已保存。\n")


#######################################################################################
###(9) 可视化总体样本及按照植物目（Order）进行分组，并对每个植物目分别进行细菌、真菌、原生生物的Mantel和偏Mantel检验

###y轴表示不同目，x轴表示不同距离矩阵，不同界微生物分页
###颜色表示P值，气泡大小表示Mantel'R

###
beta_all_mantel <- read.csv(file = "./stats/ASV_all_groups_mantel_and_geo_controlled_partial_mantel_test.csv", fileEncoding = "GBK")
beta_class_mantel <- read.csv(file = "./stats/ASV_mantel_results_by_plant_class_with_partial.csv", fileEncoding = "GBK")
beta_orders_mantel <- read.csv(file = "./stats/ASV_mantel_results_by_plant_order_with_partial.csv", fileEncoding = "GBK")

beta_all_mantel$Plant_Group <- "All"
names(beta_all_mantel)[c(1)] <- c("Microbe_Group")
names(beta_class_mantel)[1] <- c("Plant_Group")
names(beta_orders_mantel)[c(1)] <- c("Plant_Group")

names(beta_class_mantel)
names(beta_orders_mantel)

beta_all_mantel <- as.data.frame(rbind(beta_all_mantel[,c(7,1:6)], beta_class_mantel[,-4], beta_orders_mantel[,-4]))
beta_all_mantel$Distance_Metric <- gsub("^Bray$", "Bray-Curtis", beta_all_mantel$Distance_Metric)
beta_all_mantel$Distance_Metric <- gsub("^Weighted$", "Weighted UniFrac", beta_all_mantel$Distance_Metric)
beta_all_mantel$Distance_Metric <- gsub("^Unweighted$", "Unweighted UniFrac", beta_all_mantel$Distance_Metric)
beta_all_mantel$Microbe_Group <- gsub("^Protist$", "Protists", beta_all_mantel$Microbe_Group)


beta_all_mantel$Distance_Metric <- factor(beta_all_mantel$Distance_Metric, levels = c("Jaccard", "Bray-Curtis", "Unweighted UniFrac", "Weighted UniFrac"))
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

name <- "./plots/Mantel_tests_between_microbial_community_similarity_and_plant_phylogenetic_relatedness"
width <- 12
height <- 10
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

name <- "./plots/Partial_Mantel_tests_between_microbial_community_similarity_and_plant_phylogenetic_relatedness"
width <- 12
height <- 10
ggsave(paste0(name, ".png"), partial_mantel_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), partial_mantel_plot, width = width, height = height, units = "cm")

write.csv(beta_all_mantel, file = "./stats/ASV_mantel_and_partial_mantel_results_all_and_by_plant_class_order.csv", fileEncoding = "GBK", row.names = FALSE)
