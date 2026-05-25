### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

## 加载和整理数据
plant_colors <- read.csv("./data/Tree_top10_order_color.csv")
metadata <- read.csv('./data/tree_metadata_merge_info.csv')
row.names(metadata) <- metadata$TreeID
metadata$Order <- ifelse(metadata$Order %in% plant_colors$Order, metadata$Order, "Others")
metadata$TreeID <- paste("RS-", metadata$TreeID, sep = "")
row.names(metadata) <- metadata$TreeID

### 加载ko相对丰度数据（TPM）
tpm_ko_bacteria <- read.csv('./data/Bacteria_core0.2_tpm.csv', row.names = 1, check.names = FALSE)
tpm_ko_fungi <- read.csv('./data/Fungi_core0.2_tpm.csv', row.names = 1, check.names = FALSE)
tpm_ko_protists <- read.csv('./data/Protists_core0.2_tpm.csv', row.names = 1, check.names = FALSE)

### 加载ko绝对丰度数据
qa_ko_bacteria <- read.csv('./data/Bacteria_core0.2_quantitative.csv', row.names = 1, check.names = FALSE)
qa_ko_fungi <- read.csv('./data/Fungi_core0.2_quantitative.csv', row.names = 1, check.names = FALSE)
qa_ko_protists <- read.csv('./data/Protists_core0.2_quantitative.csv', row.names = 1, check.names = FALSE)


### 加载ko绝对丰度数据，log转化的
logqa_ko_bacteria <- log(qa_ko_bacteria+1)
logqa_ko_fungi <- log(qa_ko_fungi+1)
logqa_ko_protists <- log(qa_ko_protists+1)

###################################################################
### 批量生成lefse文件

# 创建数据框列表
data_list <- list(
  tpm_bacteria = tpm_ko_bacteria,
  tpm_fungi = tpm_ko_fungi,
  tpm_protists = tpm_ko_protists,
  qa_bacteria = qa_ko_bacteria,
  qa_fungi = qa_ko_fungi,
  qa_protists = qa_ko_protists,
  logqa_bacteria = logqa_ko_bacteria,
  logqa_fungi = logqa_ko_fungi,
  logqa_protists = logqa_ko_protists
)

# 对应的文件名
file_names <- c(
  "tpm_ko_bacteria_lefse.txt",
  "tpm_ko_fungi_lefse.txt", 
  "tpm_ko_protists_lefse.txt",
  "qa_ko_bacteria_lefse.txt",
  "qa_ko_fungi_lefse.txt",
  "qa_ko_protists_lefse.txt",
  "logqa_ko_bacteria_lefse.txt",
  "logqa_ko_fungi_lefse.txt",
  "logqa_ko_protists_lefse.txt"
)

# 循环生成lefse文件
for (i in 1:length(data_list)) {
  data_name <- names(data_list)[i]
  data <- data_list[[i]]
  file_name <- file_names[i]
  
  # 创建lefse格式的数据
  lefse_data <- rbind(
    Order = metadata[names(data),]$Order,
    TreeID = metadata[names(data),]$TreeID,
    data
  )
  
  # 写入文件
  write.table(lefse_data, 
              file = paste0("./data/", file_name), 
              sep = "\t", 
              col.names = FALSE, 
              quote = FALSE)
  
  cat("已生成:", file_name, "\n")
}

cat("所有lefse文件已生成完成！\n")

#############################################################
####换成两两比较

# 循环生成
# 获取所有Order（除了Fabales）
other_orders <- setdiff(plant_colors$Order, "Fabales")

# 循环处理每个数据集
for (data_name in names(data_list)) {
  data <- data_list[[data_name]]
  
  # 创建lefse格式的数据
  lefse_data <- rbind(
    Order = metadata[colnames(data), ]$Order,
    TreeID = metadata[colnames(data), ]$TreeID,
    data
  )
  
  # 为每个比较组创建文件
  for (i in other_orders) {
    # 选择样本
    samples_to_keep <- lefse_data["Order", ] %in% c("Fabales", i)
    
    # 创建子集
    subset_data <- lefse_data[, samples_to_keep, drop = FALSE]
    
    # 修改Order标签，确保Fabales在前
    subset_data["Order", ] <- ifelse(
      subset_data["Order", ] == "Fabales",
      "1Fabales",
      subset_data["Order", ]
    )
    
    # 创建文件名
    file_name <- paste0(
      "./data/Fabales_", 
      i, 
      "_", 
      gsub("_ko", "", data_name),  # 移除_ko
      "_lefse.txt"
    )
    
    # 写入文件
    write.table(
      subset_data,
      file = file_name,
      sep = "\t",
      col.names = FALSE,
      quote = FALSE
    )
    
    cat("已生成:", basename(file_name), "\n")
  }
}

cat("所有两两比较文件已生成完成！\n")
