### Title: Get Network Robustness
### Author: Mingxing Wang
### Date: 2024-07-01
### Reference: https://doi.org/10.1038/s41558-021-00989-9


# cor_matrix: 数据框，相关性matrix
# feature_tab: 数据框，物种特征表
# r_threshold: 0到1的浮点型数值，NULL表示输入的cor_matrix已经按相关性阈值过滤
# weighted: 布尔型，是否计算加权后的robustness鲁棒性结果
# table_type: 字符型，当weighted=TRUE时该参数有效，可选类型为absolute、relative、rarefied
# sample_depth: 整型数值，当table_type='rarefied'时该参数有效
# n_permutation: 整型数值，置换次数，即重复计算次数，默认为100次
# rm_percent: 0到1的浮点型数值，随机删除网络节点的占比
# rm_feature: 字符向量型，需要删除的目标网络节点名称


data_processing <- function(cor_matrix, feature_tab, r_threshold, table_type, sample_depth) {
    if (!is.null(r_threshold)) {
        cor_matrix <- cor_matrix * (abs(cor_matrix) >= r_threshold)  # only keep links above the r_threshold
    }
    diag(cor_matrix) <- 0  # remove links for self to self
    cor_matrix <- cor_matrix[colSums(abs(cor_matrix)) > 0, colSums(abs(cor_matrix)) > 0]
    
    if (table_type == 'rarefied') {
        mean_feature_abund <- rowMeans(feature_tab) / sample_depth
    } else {
        mean_feature_abund <- rowMeans(feature_tab)
    }
    mean_feature_abund <- mean_feature_abund[colnames(cor_matrix)]
    
    return(list(cor_matrix, mean_feature_abund))
}


remove_feature <- function(cor_matrix, rm_list, weighted, mean_feature_abund) {
    cor_matrix[rm_list,] <- 0
    cor_matrix[,rm_list] <- 0
    pruned_matrix <- cor_matrix
    if (weighted) {cor_matrix <- cor_matrix * mean_feature_abund}
    
    # remove the features that have negative interaction or no interaction with others
    rm_node <- which(colMeans(cor_matrix) <= 0)
    remain_percent <- (nrow(cor_matrix) - length(rm_node)) / nrow(cor_matrix)
    return(remain_percent)
    
    # For simplicity, we only consider the immediate effects of removing the 'rm_list' feature.
    # We don't consider the sequential effects of extinction of the 'rm_node' feature.
    # Write the pruned network:
    # pruned_matrix[rm_node, ] <- 0
    # pruned_matrix[, rm_node] <- 0
    # write.csv(pruned_matrix, 'pruned_matrix.csv')
}


random_remove <- function(cor_matrix, rm_percent, weighted, mean_feature_abund) {
    rm_list <- sample(1:nrow(cor_matrix), round(nrow(cor_matrix) * rm_percent))
    remain_percent <- remove_feature(cor_matrix, rm_list, weighted, mean_feature_abund)
    return(remain_percent)
}


target_remove <- function(cor_matrix, rm_feature, weighted, mean_feature_abund) {
    rm_list <- rm_feature
    remain_percent <- remove_feature(cor_matrix, rm_list, weighted, mean_feature_abund)
    return(remain_percent)
}


random_remove_simulation <- function(cor_matrix, n_permutation, weighted, mean_feature_abund) {
    t(sapply(seq(0.05, 1.00, 0.05), function(percent) {
        remain_percent <- sapply(1:n_permutation, function(i) {
            random_remove(cor_matrix, percent, weighted, mean_feature_abund)
        })
        
        remove_percent <- percent
        remain_mean <- mean(remain_percent)
        remain_sd <- sd(remain_percent)
        remain_se <- sd(remain_percent) / (n_permutation ^ 0.5)
        
        result <- c(remove_percent, remain_mean, remain_sd, remain_se)
        names(result) <- c('remove_percent', 'remain_mean', 'remain_sd', 'remain_se')
        return(result)
    }))
}


get_random_remove_robustness <- function(
        cor_matrix, feature_tab, r_threshold = NULL, weighted = FALSE, table_type = NULL, sample_depth = NULL, 
        n_permutation = 100, rm_percent = 0.5) {
    data_list <- data_processing(cor_matrix, feature_tab, r_threshold, table_type, sample_depth)
    
    remain_percent <- sapply(1:n_permutation, function(i) {
        random_remove(data_list[[1]], rm_percent, weighted, data_list[[2]])
    })
    
    remove_simulation <- random_remove_simulation(data_list[[1]], n_permutation, weighted, data_list[[2]])
    
    return(list(remain_percent = remain_percent, remove_simulation = remove_simulation))
}


get_target_remove_robustness <- function(
        cor_matrix, feature_tab, r_threshold = NULL, weighted = FALSE, table_type = NULL, sample_depth = NULL, 
        rm_feature = NULL) {
    data_list <- data_processing(cor_matrix, feature_tab, r_threshold, table_type, sample_depth)
    
    remain_percent <- target_remove(data_list[[1]], rm_feature, weighted, data_list[[2]])

    return(remain_percent)
}
