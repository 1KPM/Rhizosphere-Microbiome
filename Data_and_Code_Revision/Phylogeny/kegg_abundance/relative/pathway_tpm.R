# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
library(tidyverse)

#Import Data
df_bac <- read.csv("./Bacteria_RS_tpm_abundance.csv",check.names = F,row.names = 1,header = T)
df_fun <- read.csv("./Fungi_RS_tpm_abundance.csv",check.names = F,row.names = 1,header = T)
df_pro <- read.csv("./Protist_RS_tpm_abundance.csv",check.names = F,row.names = 1,header = T)
ko_pathway <- read.csv("../kegg_pathway_for_1KPM.csv",check.names = F,row.names = 1,header = T)

#####################################################################
#分类汇总到Pathway和super pathway的TPM (相对丰度)
df_bac_tpm <- rownames_to_column(df_bac, var = "KO")
df_bac_tpm_pathway <- merge.data.frame(ko_pathway[c("KO","level1", "level2", "pathway")],
                                     df_bac_tpm, by = "KO", all.y = T)

df_fun_tpm <- rownames_to_column(df_fun, var = "KO")
df_fun_tpm_pathway <- merge.data.frame(ko_pathway[c("KO","level1", "level2", "pathway")],
                                     df_fun_tpm, by = "KO", all.y = T)

df_pro_tpm <- rownames_to_column(df_pro, var = "KO")
df_pro_tpm_pathway <- merge.data.frame(ko_pathway[c("KO","level1", "level2", "pathway")],
                                     df_pro_tpm, by = "KO", all.y = T)

# 按照 level1, level2, pathway 分组汇总
# 对每个样本列进行求和
# 使用循环处理
data_frames <- list(
  bac = df_bac_tpm_pathway,
  fun = df_fun_tpm_pathway,
  pro = df_pro_tpm_pathway
)

for (df_name in names(data_frames)) {
  df <- data_frames[[df_name]]
  
  # 创建文件名基础
  base_name <- paste0("Summary_", df_name, "_tpm")
  
  # 处理三个level
  summaries <- list()
  
  # Level 1
  level1 <- df %>%
    group_by(level1) %>%
    summarise(across(starts_with("RS-"), ~ sum(., na.rm = TRUE)), .groups = "drop") %>%
    slice_head(n = -1)  # 删除最后一行
  names(level1)[1] <- "Pathway"
  summaries[[1]] <- level1
  
  # Level 2
  level2 <- df %>%
    group_by(level1, level2) %>%
    summarise(across(starts_with("RS-"), ~ sum(., na.rm = TRUE)), .groups = "drop") %>%
    slice_head(n = -1)
  level2$Pathway <- paste(level2$level1, level2$level2, sep = "|")
  level2 <- level2 %>% select(Pathway, starts_with("RS-"))
  summaries[[2]] <- level2
  
  # Level 3
  level3 <- df %>%
    group_by(level1, level2, pathway) %>%
    summarise(across(starts_with("RS-"), ~ sum(., na.rm = TRUE)), .groups = "drop") %>%
    slice_head(n = -1)
  level3$Pathway <- paste(level3$level1, level3$level2, level3$pathway, sep = "|")
  level3 <- level3 %>% select(Pathway, starts_with("RS-"))
  summaries[[3]] <- level3
  
  # 合并所有level
  all_pathway <- bind_rows(summaries)
  
  # 保存结果
  write.csv(all_pathway, paste0(base_name, "_all_pathway.csv"), row.names = FALSE)
  
  cat("Processed", df_name, "\n")
}


#############################
###Supplementary Information Table 4a-c
#Import Data
df_bac <- read.csv("./Bacteria_core0.2_tpm.csv",check.names = F, header = T)
df_fun <- read.csv("./Fungi_core0.2_tpm.csv",check.names = F, header = T)
df_pro <- read.csv("./Protists_core0.2_tpm.csv",check.names = F, header = T)
ko_des <- read.csv("../K_gene_name.csv",check.names = F,row.names = 1,header = T)

names(df_bac)[1] <- "KO"
names(df_fun)[1] <- "KO"
names(df_pro)[1] <- "KO"

df_bac$KO <- gsub("b", "", df_bac$KO)
df_fun$KO <- gsub("f", "", df_fun$KO)
df_pro$KO <- gsub("p", "", df_pro$KO)

df_bac_des <- merge(ko_des, df_bac, by = "KO", all.y = T)
df_fun_des <- merge(ko_des, df_fun, by = "KO", all.y = T)
df_pro_des <- merge(ko_des, df_pro, by = "KO", all.y = T)

# 保存结果
write.csv(df_bac_des, "Bacteria_core0.2_tpm_des.csv", row.names = FALSE)
write.csv(df_fun_des, "Fungi_core0.2_tpm_des.csv", row.names = FALSE)
write.csv(df_pro_des, "Protists_core0.2_tpm_des.csv", row.names = FALSE)

