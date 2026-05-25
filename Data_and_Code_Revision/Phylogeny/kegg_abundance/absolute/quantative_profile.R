# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
library(tidyverse)

#Import Data
df_bac <- read.delim("./reads/Bacteria_KO_number.tsv",check.names = F,row.names = 1,header = T)
df_fun <- read.delim("./reads/Fungi_KO_number.tsv",check.names = F,row.names = 1,header = T)
df_pro <- read.delim("./reads/Protist_KO_number.tsv",check.names = F,row.names = 1,header = T)
just_df <- read.csv("./Rhizobiales_amplicon_meta.csv")
ko_pathway <- read.csv("../kegg_pathway_for_1KPM.csv",check.names = F,row.names = 1,header = T)

#Quantitative Profile
df_bac_q <- df_bac %>% 
  select(any_of(just_df$SampleName)) %>%
  mutate(across(everything(), 
                ~ . * just_df$ratio[just_df$SampleName == cur_column()],
                .names = "{.col}"))
write.csv(df_bac_q,"./Bacteria_RS_quantitative_abundance.csv")

df_fun_q <- df_fun %>% 
  select(any_of(just_df$SampleName)) %>%
  mutate(across(everything(), 
                ~ . * just_df$ratio[just_df$SampleName == cur_column()],
                .names = "{.col}"))
write.csv(df_fun_q,"./Fungi_RS_quantitative_abundance.csv")

df_pro_q <- df_pro %>% 
  select(any_of(just_df$SampleName)) %>%
  mutate(across(everything(), 
                ~ . * just_df$ratio[just_df$SampleName == cur_column()],
                .names = "{.col}"))
write.csv(df_pro_q,"./Protists_RS_quantitative_abundance.csv")

#####################################################################
#分类汇总到Pathway和super pathway的丰度
df_bac_qa <- rownames_to_column(df_bac_q, var = "KO")
df_bac_qa_pathway <- merge.data.frame(ko_pathway[c("KO","level1", "level2", "pathway")],
                                     df_bac_qa, by = "KO", all.y = T)

df_fun_qa <- rownames_to_column(df_fun_q, var = "KO")
df_fun_qa_pathway <- merge.data.frame(ko_pathway[c("KO","level1", "level2", "pathway")],
                                     df_fun_qa, by = "KO", all.y = T)

df_pro_qa <- rownames_to_column(df_pro_q, var = "KO")
df_pro_qa_pathway <- merge.data.frame(ko_pathway[c("KO","level1", "level2", "pathway")],
                                     df_pro_qa, by = "KO", all.y = T)

# 按照 level1, level2, pathway 分组汇总
# 对每个样本列进行求和
# 使用循环处理
data_frames <- list(
  bac = df_bac_qa_pathway,
  fun = df_fun_qa_pathway,
  pro = df_pro_qa_pathway
)

for (df_name in names(data_frames)) {
  df <- data_frames[[df_name]]
  
  # 创建文件名基础
  base_name <- paste0("Summary_", df_name, "_q")
  
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
  
  # 创建log转换版本
  log_pathway <- all_pathway
  log_pathway <- log_pathway %>%
    mutate(across(starts_with("RS-"), ~ log(. + 1)))
  
  # 保存结果
  write.csv(all_pathway, paste0(base_name, "_all_pathway.csv"), row.names = FALSE)
  write.csv(log_pathway, paste0(base_name, "_log_all_pathway.csv"), row.names = FALSE)
  
  cat("Processed", df_name, "\n")
}


#####################################################################
##Core Profile
#Core 0.2
df_bac_r <- read.csv("../relative/Bacteria_core0.2_tpm.csv",row.names = 1,check.names = F,header = T)
df_bac_q_core <- df_bac_q %>% filter(rownames(.) %in% rownames(df_bac_r))
write.csv(df_bac_q_core,"./Bacteria_core0.2_quantitative.csv")

df_fun_r <- read.csv("../relative/Fungi_core0.2_tpm.csv",row.names = 1,check.names = F,header = T)
df_fun_q_core <- df_fun_q %>% filter(rownames(.) %in% rownames(df_fun_r))
write.csv(df_fun_q_core,"./Fungi_core0.2_quantitative.csv")

df_pro_r <- read.csv("../relative/Protist_core0.2_tpm.csv",row.names = 1,check.names = F,header = T)
df_pro_q_core <- df_pro_q %>% filter(rownames(.) %in% rownames(df_pro_r))
write.csv(df_pro_q_core,"./Protists_core0.2_quantitative.csv")

#Core 0.5
df_bac_r <- read.csv("../relative/Bacteria_core0.5_tpm.csv",row.names = 1,check.names = F,header = T)
df_bac_q_core <- df_bac_q %>% filter(rownames(.) %in% rownames(df_bac_r))
write.csv(df_bac_q_core,"./Bacteria_core0.5_quantitative.csv")

df_fun_r <- read.csv("../relative/Fungi_core0.5_tpm.csv",row.names = 1,check.names = F,header = T)
df_fun_q_core <- df_fun_q %>% filter(rownames(.) %in% rownames(df_fun_r))
write.csv(df_fun_q_core,"./Fungi_core0.5_quantitative.csv")

df_pro_r <- read.csv("../relative/Protist_core0.5_tpm.csv",row.names = 1,check.names = F,header = T)
df_pro_q_core <- df_pro_q %>% filter(rownames(.) %in% rownames(df_pro_r))
write.csv(df_pro_q_core,"./Protists_core0.5_quantitative.csv")