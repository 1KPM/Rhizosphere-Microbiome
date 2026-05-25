### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

## 加载和整理数据
plant_colors <- read.csv("./data/Tree_top10_order_color.csv")
metadata <- read.csv('./data/rhizosphere_metadata_merge_info.csv')
row.names(metadata) <- metadata$FileID
metadata$Order <- ifelse(metadata$Order %in% plant_colors$Order, metadata$Order, "Others")

### 加载拷贝数据
fin_copies_all <- read.csv('./data/All_sample_copies.csv', row.names = 1)

### 加载分层taxa数据
all_taxa_16s <- read.csv('./data/Bacterial_all_taxa_16S_reads.csv', row.names = 1, check.names = FALSE)
all_taxa_its <- read.csv('./data/Fungal_all_taxa_ITS_reads.csv', row.names = 1, check.names = FALSE)
all_taxa_pro <- read.csv('./data/Protistan_all_taxa_18S_reads.csv', row.names = 1, check.names = FALSE)

### 计算分层taxa的绝对丰度
### 16S
all_taxa_qa_16s <- all_taxa_16s[row.names(fin_copies_all),]

for (i in row.names(all_taxa_qa_16s)) {
  all_taxa_qa_16s[i,-1] <- round(all_taxa_qa_16s[i,-1] / fin_copies_all[i, 'Bacteria_16S'] * fin_copies_all[i, 'Copies_16S'])
  all_taxa_qa_16s[i,-1] <- log(all_taxa_qa_16s[i,-1] + 1)
}
write.csv(all_taxa_qa_16s, './data/Bacterial_all_taxa_16S_qa_log.csv', quote = F)

### ITS
all_taxa_qa_its <- all_taxa_its[row.names(fin_copies_all),]
for (i in row.names(all_taxa_qa_its)) {
  all_taxa_qa_its[i,-1] <- round(all_taxa_qa_its[i,-1] / fin_copies_all[i, 'Fungi_ITS'] * fin_copies_all[i, 'Copies_ITS'])
  all_taxa_qa_its[i,-1] <- log(all_taxa_qa_its[i,-1] + 1)
}
write.csv(all_taxa_qa_its, './data/Fungal_all_taxa_ITS_qa_log.csv', quote = F)

### Protist
all_taxa_qa_pro <- all_taxa_pro[row.names(fin_copies_all),]
for (i in row.names(all_taxa_qa_pro)) {
  all_taxa_qa_pro[i,-1] <- round(all_taxa_qa_pro[i,-1] / fin_copies_all[i, 'Protist_18S'] * fin_copies_all[i, 'Copies_Protist'])
  all_taxa_qa_pro[i,-1] <- log(all_taxa_qa_pro[i,-1] + 1)
}
write.csv(all_taxa_qa_pro, './data/Protistan_all_taxa_18S_qa_log.csv', quote = F)

###################################################################
###准备lefse文件
###bacteria
all_taxa_read_16s <- merge(metadata[row.names(fin_copies_all),c("TreeID","Order")],
                           all_taxa_16s[row.names(fin_copies_all),], 
                           by = "row.names")

all_taxa_read_16s_lefse <- as.data.frame(t(all_taxa_read_16s[,-c(1,2,4)]))
row.names(all_taxa_read_16s_lefse) <- gsub(";[pcofg]__", "|", row.names(all_taxa_read_16s_lefse))
row.names(all_taxa_read_16s_lefse) <- gsub("d__", "", row.names(all_taxa_read_16s_lefse))

write.table(all_taxa_read_16s_lefse, file = "./data/Bacterial_all_taxa_16S_reads_lefse.txt", 
            sep = "\t", col.names = F, quote = F)

###
all_taxa_qa_16s <- merge(metadata[row.names(fin_copies_all),c("TreeID","Order")],
                           all_taxa_qa_16s[row.names(fin_copies_all),], 
                           by = "row.names")

all_taxa_qa_16s_lefse <- as.data.frame(t(all_taxa_qa_16s[,-c(1,2,4)]))
row.names(all_taxa_qa_16s_lefse) <- gsub(";[pcofg]__", "|", row.names(all_taxa_qa_16s_lefse))
row.names(all_taxa_qa_16s_lefse) <- gsub("d__", "", row.names(all_taxa_qa_16s_lefse))


write.table(all_taxa_qa_16s_lefse, file = "./data/Bacterial_all_taxa_16S_qa_log_lefse.txt", 
            sep = "\t", col.names = F, quote = F)


###fungi
all_taxa_read_its <- merge(metadata[row.names(fin_copies_all),c("TreeID","Order")],
                           all_taxa_its[row.names(fin_copies_all),], 
                           by = "row.names")

all_taxa_read_its_lefse <- as.data.frame(t(all_taxa_read_its[,-c(1,2,4)]))
row.names(all_taxa_read_its_lefse) <- gsub(";[pcofg]__", "|", row.names(all_taxa_read_its_lefse))
row.names(all_taxa_read_its_lefse) <- gsub("k__", "", row.names(all_taxa_read_its_lefse))


write.table(all_taxa_read_its_lefse, file = "./data/Fungal_all_taxa_ITS_reads_lefse.txt", 
            sep = "\t", col.names = F, quote = F)

###
all_taxa_qa_its <- merge(metadata[row.names(fin_copies_all),c("TreeID","Order")],
                         all_taxa_qa_its[row.names(fin_copies_all),], 
                         by = "row.names")

all_taxa_qa_its_lefse <- as.data.frame(t(all_taxa_qa_its[,-c(1,2,4)]))
row.names(all_taxa_qa_its_lefse) <- gsub(";[pcofg]__", "|", row.names(all_taxa_qa_its_lefse))
row.names(all_taxa_qa_its_lefse) <- gsub("k__", "", row.names(all_taxa_qa_its_lefse))


write.table(all_taxa_qa_its_lefse, file = "./data/Fungal_all_taxa_ITS_qa_log_lefse.txt", 
            sep = "\t", col.names = F, quote = F)


###protist
all_taxa_read_pro <- merge(metadata[row.names(fin_copies_all),c("TreeID","Order")],
                           all_taxa_pro[row.names(fin_copies_all),], 
                           by = "row.names")

all_taxa_read_pro_lefse <- as.data.frame(t(all_taxa_read_pro[,-c(1,2,4)]))
row.names(all_taxa_read_pro_lefse) <- gsub(";", "|", row.names(all_taxa_read_pro_lefse))
row.names(all_taxa_read_pro_lefse) <- gsub("Eukaryota", "Protists", row.names(all_taxa_read_pro_lefse))


write.table(all_taxa_read_pro_lefse, file = "./data/Protistan_all_taxa_18S_reads_lefse.txt",
            sep = "\t", col.names = F, quote = F)

###
all_taxa_qa_pro <- merge(metadata[row.names(fin_copies_all),c("TreeID","Order")],
                         all_taxa_qa_pro[row.names(fin_copies_all),], 
                         by = "row.names")

all_taxa_qa_pro_lefse <- as.data.frame(t(all_taxa_qa_pro[,-c(1,2,4)]))
row.names(all_taxa_qa_pro_lefse) <- gsub(";", "|", row.names(all_taxa_qa_pro_lefse))
row.names(all_taxa_qa_pro_lefse) <- gsub("Eukaryota", "Protists", row.names(all_taxa_qa_pro_lefse))


write.table(all_taxa_qa_pro_lefse, file = "./data/Protistan_all_taxa_18S_qa_log_lefse.txt", 
            sep = "\t", col.names = F, quote = F)


####换成两两比较
# 16S
for (i in plant_colors$Order[-1]) {
  # 创建数据子集
  subset_data <- all_taxa_qa_16s_lefse[, all_taxa_qa_16s_lefse["Order",] %in% c("Fabales", i)]
  subset_data["Order",] <- gsub("Fabales", "1Fabales", subset_data["Order",]) # 确保lefse顺序都是Fabales在前
  # 写入文件
  write.table(subset_data, 
              file = paste0("./data/Bacterial_Fabales_", i, "_taxa_16S_qa_log_lefse.txt"), 
              sep = "\t", 
              col.names = FALSE, 
              quote = FALSE)
}

#its
for (i in plant_colors$Order[-1]) {
  # 创建数据子集
  subset_data <- all_taxa_qa_its_lefse[, all_taxa_qa_its_lefse["Order",] %in% c("Fabales", i)]
  subset_data["Order",] <- gsub("Fabales", "1Fabales", subset_data["Order",]) # 确保lefse顺序都是Fabales在前
  # 写入文件
  write.table(subset_data, 
              file = paste0("./data/Fungal_Fabales_", i, "_taxa_ITS_qa_log_lefse.txt"), 
              sep = "\t", 
              col.names = FALSE, 
              quote = FALSE)
}

#18S
for (i in plant_colors$Order[-1]) {
  # 创建数据子集
  subset_data <- all_taxa_qa_pro_lefse[, all_taxa_qa_pro_lefse["Order",] %in% c("Fabales", i)]
  subset_data["Order",] <- gsub("Fabales", "1Fabales", subset_data["Order",]) # 确保lefse顺序都是Fabales在前
  # 写入文件
  write.table(subset_data, 
              file = paste0("./data/Protistan_Fabales_", i, "_taxa_18S_qa_log_lefse.txt"), 
              sep = "\t", 
              col.names = FALSE, 
              quote = FALSE)
}




