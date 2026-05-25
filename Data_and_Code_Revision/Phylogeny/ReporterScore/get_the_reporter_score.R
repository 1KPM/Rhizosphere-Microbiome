pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(dplyr)
library(tibble)
library(ReporterScore)

###Bacteria
set.seed(2024)
KO_name <- read.csv("./data/K_gene_name.csv",header =  T,row.names = 1)
order_list <- c("Asparagales","Arecales","Fabales","Rosales","Lamiales",
                "Malpighiales","Sapindales","Gentianales","Malvales","Myrtales")
KO_abun <- read.csv("../kegg_abundance/absolute/Bacteria_core0.2_quantitative.csv",row.names = 1,check.names = F)
KO_abun <- log(KO_abun + 1)
metadata<-read.csv(file="./data/metadata.csv",header=T)
metadata_filtered <- metadata %>%
  filter(TreeID %in% colnames(KO_abun)) %>%
  mutate(Group = ifelse(Order %in% order_list, as.character(Order), "Others")) %>%
  select(TreeID,Group) %>%
  column_to_rownames("TreeID")

pathway_name <- read.csv(paste0("./data/kegg_pathway_for_1KPM.csv"),header =  T)
pathway_name$pathway <- sub("\\s*\\[.*", "", pathway_name$pathway)
other_order <- c("Asparagales","Arecales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Malvales","Myrtales","Others")
all_pathway_enrich <- data.frame()
all_enrich_KO <- data.frame()

for (o in other_order){
  metadata_sub <- metadata_filtered %>% filter(Group == "Fabales" | Group == o)
  KO_abun_sub <- KO_abun[, rownames(metadata_sub)]
  
  reporter_res <- reporter_score(KO_abun_sub, "Group", metadata_sub, 
                                 mode = "directed", method = "wilcox.test",
                                 perm = 999, p.adjust.method1 = "fdr", p.adjust.method2 = "fdr")
  
  pathway_enrich <- reporter_res$reporter_s %>% 
    filter(Description %in% pathway_name$pathway) %>% 
    mutate(GroupComparison = paste0("Fabales vs ", o))
  
  all_pathway_enrich <- bind_rows(all_pathway_enrich, pathway_enrich)
  
  enrich_KO <- reporter_res$ko_stat %>% 
    merge(KO_name, by.x = "KO_id", by.y = "KO", all.x = TRUE) %>% 
    mutate(GroupComparison = paste0("Fabales vs ", o)) %>%
    rename(
      average_Order = !!sym(paste0("average_", o)),
      sd_Order = !!sym(paste0("sd_", o))
    ) %>% 
    rename(
      average_Fabales = average_Fabales,
      sd_Fabales = sd_Fabales
    )
  
  all_enrich_KO <- bind_rows(all_enrich_KO, enrich_KO)
}

write.csv(all_pathway_enrich, file = "./data/Bacteria_Fabales_vs_all_orders_pathway_list.csv", row.names = FALSE)
write.csv(all_enrich_KO, file = "./data/Bacteria_Fabales_vs_all_orders_KO_stat.csv", row.names = FALSE)

all_enrich_KO <- all_enrich_KO %>% 
  mutate(
    across(c(diff_mean, sign, Z_score), ~ 
             if_else(
               !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
               . * -1,
               .
             )
    ),
    type = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      recode(type, "Depleted" = "Enriched", "Enriched" = "Depleted"),
      type
    ),
    Significantly = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      recode(Significantly, "Depleted" = "Enriched", "Enriched" = "Depleted"),
      Significantly
    )
  )
write.csv(all_enrich_KO, file = "./data/Bacteria_Fabales_vs_all_orders_KO_stat1.csv", row.names = FALSE)

all_pathway_enrich <- read.csv("./data/Bacteria_Fabales_vs_all_orders_pathway_list.csv")
all_pathway_enrich <- all_pathway_enrich %>% 
  mutate(
    across(c(Z_score, BG_Mean, ReporterScore), ~ 
             if_else(
               !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
               . * -1,
               .
             )
    ),
    temp = Significant_up_num,
    Significant_up_num = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      Significant_down_num,
      Significant_up_num
    ),
    Significant_down_num = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      temp,
      Significant_down_num
    ),
    temp = NULL
  )
write.csv(all_pathway_enrich, file = "./data/Bacteria_Fabales_vs_all_orders_pathway_list1.csv", row.names = FALSE)





###Fungi
set.seed(2024)
KO_name <- read.csv("./data/K_gene_name.csv",header =  T,row.names = 1)
order_list <- c("Asparagales","Arecales","Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Malvales","Myrtales")
KO_abun <- read.csv("../kegg_abundance/absolute/Fungi_core0.2_quantitative.csv",row.names = 1,check.names = F)
KO_abun <- log(KO_abun + 1)
metadata<-read.csv(file="./data/metadata.csv",header=T)
metadata_filtered <- metadata %>%
  filter(TreeID %in% colnames(KO_abun)) %>%
  mutate(Group = ifelse(Order %in% order_list, as.character(Order), "Others")) %>%
  select(TreeID,Group) %>%
  column_to_rownames("TreeID")

pathway_name <- read.csv(paste0("./data/kegg_pathway_for_1KPM.csv"),header =  T)
pathway_name$pathway <- sub("\\s*\\[.*", "", pathway_name$pathway)
other_order <- c("Asparagales","Arecales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Malvales","Myrtales","Others")
all_pathway_enrich <- data.frame()
all_enrich_KO <- data.frame()

for (o in other_order){
  metadata_sub <- metadata_filtered %>% filter(Group == "Fabales" | Group == o)
  KO_abun_sub <- KO_abun[, rownames(metadata_sub)]
  
  reporter_res <- reporter_score(KO_abun_sub, "Group", metadata_sub, 
                                 mode = "directed", method = "wilcox.test",
                                 perm = 999, p.adjust.method1 = "fdr", p.adjust.method2 = "fdr")
  
  pathway_enrich <- reporter_res$reporter_s %>% 
    filter(Description %in% pathway_name$pathway) %>% 
    mutate(GroupComparison = paste0("Fabales vs ", o))
  
  all_pathway_enrich <- bind_rows(all_pathway_enrich, pathway_enrich)
  
  enrich_KO <- reporter_res$ko_stat %>% 
    merge(KO_name, by.x = "KO_id", by.y = "KO", all.x = TRUE) %>% 
    mutate(GroupComparison = paste0("Fabales vs ", o)) %>%
    rename(
      average_Order = !!sym(paste0("average_", o)),
      sd_Order = !!sym(paste0("sd_", o))
    ) %>% 
    rename(
      average_Fabales = average_Fabales,
      sd_Fabales = sd_Fabales
    )
  
  all_enrich_KO <- bind_rows(all_enrich_KO, enrich_KO)
}

write.csv(all_pathway_enrich, file = "./data/Fungi_Fabales_vs_all_orders_pathway_list.csv", row.names = FALSE)
write.csv(all_enrich_KO, file = "./data/Fungi_Fabales_vs_all_orders_KO_stat.csv", row.names = FALSE)

all_enrich_KO <- all_enrich_KO %>% 
  mutate(
    across(c(diff_mean, sign, Z_score), ~ 
             if_else(
               !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
               . * -1,
               .
             )
    ),
    type = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      recode(type, "Depleted" = "Enriched", "Enriched" = "Depleted"),
      type
    ),
    Significantly = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      recode(Significantly, "Depleted" = "Enriched", "Enriched" = "Depleted"),
      Significantly
    )
  )
write.csv(all_enrich_KO, file = "./data/Fungi_Fabales_vs_all_orders_KO_stat1.csv", row.names = FALSE)

all_pathway_enrich <- read.csv("./data/Fungi_Fabales_vs_all_orders_pathway_list.csv")
all_pathway_enrich <- all_pathway_enrich %>% 
  mutate(
    across(c(Z_score, BG_Mean, ReporterScore), ~ 
             if_else(
               !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
               . * -1,
               .
             )
    ),
    temp = Significant_up_num,
    Significant_up_num = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      Significant_down_num,
      Significant_up_num
    ),
    Significant_down_num = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      temp,
      Significant_down_num
    ),
    temp = NULL
  )
write.csv(all_pathway_enrich, file = "./data/Fungi_Fabales_vs_all_orders_pathway_list1.csv", row.names = FALSE)



###Protist
set.seed(2024)
KO_name <- read.csv("./data/K_gene_name.csv",header =  T,row.names = 1)
order_list <- c("Asparagales","Arecales","Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Malvales","Myrtales")
KO_abun <- read.csv("../kegg_abundance/absolute/Protist_core0.2_quantitative.csv",row.names = 1,check.names = F)
KO_abun <- log(KO_abun + 1)
metadata<-read.csv(file="./data/metadata.csv",header=T)
metadata_filtered <- metadata %>%
  filter(TreeID %in% colnames(KO_abun)) %>%
  mutate(Group = ifelse(Order %in% order_list, as.character(Order), "Others")) %>%
  select(TreeID,Group) %>%
  column_to_rownames("TreeID")

pathway_name <- read.csv(paste0("./data/kegg_pathway_for_1KPM.csv"),header =  T)
pathway_name$pathway <- sub("\\s*\\[.*", "", pathway_name$pathway)
other_order <- c("Asparagales","Arecales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Malvales","Myrtales","Others")
all_pathway_enrich <- data.frame()
all_enrich_KO <- data.frame()

for (o in other_order){
  metadata_sub <- metadata_filtered %>% filter(Group == "Fabales" | Group == o)
  KO_abun_sub <- KO_abun[, rownames(metadata_sub)]
  
  reporter_res <- reporter_score(KO_abun_sub, "Group", metadata_sub, 
                                 mode = "directed", method = "wilcox.test",
                                 perm = 999, p.adjust.method1 = "fdr", p.adjust.method2 = "fdr")
  
  pathway_enrich <- reporter_res$reporter_s %>% 
    filter(Description %in% pathway_name$pathway) %>% 
    mutate(GroupComparison = paste0("Fabales vs ", o))
  
  all_pathway_enrich <- bind_rows(all_pathway_enrich, pathway_enrich)
  
  enrich_KO <- reporter_res$ko_stat %>% 
    merge(KO_name, by.x = "KO_id", by.y = "KO", all.x = TRUE) %>% 
    mutate(GroupComparison = paste0("Fabales vs ", o)) %>%
    rename(
      average_Order = !!sym(paste0("average_", o)),
      sd_Order = !!sym(paste0("sd_", o))
    ) %>% 
    rename(
      average_Fabales = average_Fabales,
      sd_Fabales = sd_Fabales
    )
  
  all_enrich_KO <- bind_rows(all_enrich_KO, enrich_KO)
}

write.csv(all_pathway_enrich, file = "./data/Protist_Fabales_vs_all_orders_pathway_list.csv", row.names = FALSE)
write.csv(all_enrich_KO, file = "./data/Protist_Fabales_vs_all_orders_KO_stat.csv", row.names = FALSE)

all_enrich_KO <- all_enrich_KO %>% 
  mutate(
    across(c(diff_mean, sign, Z_score), ~ 
             if_else(
               !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
               . * -1,
               .
             )
    ),
    type = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      recode(type, "Depleted" = "Enriched", "Enriched" = "Depleted"),
      type
    ),
    Significantly = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      recode(Significantly, "Depleted" = "Enriched", "Enriched" = "Depleted"),
      Significantly
    )
  )
write.csv(all_enrich_KO, file = "./data/Protist_Fabales_vs_all_orders_KO_stat1.csv", row.names = FALSE)

all_pathway_enrich <- read.csv("./data/Protist_Fabales_vs_all_orders_pathway_list.csv")
all_pathway_enrich <- all_pathway_enrich %>% 
  mutate(
    across(c(Z_score, BG_Mean, ReporterScore), ~ 
             if_else(
               !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
               . * -1,
               .
             )
    ),
    temp = Significant_up_num,
    Significant_up_num = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      Significant_down_num,
      Significant_up_num
    ),
    Significant_down_num = if_else(
      !GroupComparison %in% c("Fabales vs Asparagales", "Fabales vs Arecales"),
      temp,
      Significant_down_num
    ),
    temp = NULL
  )
write.csv(all_pathway_enrich, file = "./data/Protist_Fabales_vs_all_orders_pathway_list1.csv", row.names = FALSE)
