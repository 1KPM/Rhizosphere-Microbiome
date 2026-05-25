pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(dplyr)
library(tibble)

##Import Data
#-----------------------------------------------------------------
metadata <- read.csv("./data/metadata.csv")
kegg_info <- read.csv("./data/kegg_pathway_for_1KPM.csv")
kegg_info$L2 <- gsub(pattern = "l1_", replacement = "l2_", x = kegg_info$L2)
kegg_info$L3 <- gsub(pattern = "ko", replacement = "map", x = kegg_info$L3)
K_pathway <- read.table("./data/ko_pathway_mapping.txt")
pathway_super <- read.table("./data/level3_level2_mapping.txt")
#-----------------------------------------------------------------

#Profile Function
#-----------------------------------------------------------------
process_kegg_pathway <- function(kingdom, type, input_df) {
  unassign_df <- input_df[nrow(input_df), , drop = FALSE]
  
  #Pathway Profile
  raw_pathway_df <- input_df %>%
    rownames_to_column(var = "KO") %>%
    merge(K_pathway, by.x = "KO", by.y = "V1", all.x = FALSE) %>%
    group_by(V2) %>%
    summarise(across(matches("^RS|^S"), sum), .groups = "drop") %>% 
    column_to_rownames(var = "V2")
  
  fin_pathway_df <- raw_pathway_df %>%
    filter(rownames(.) %in% (kegg_info %>% select(L3, pathway) %>% distinct(L3, .keep_all = TRUE))$L3) %>%
    rownames_to_column(var = "L3") %>%
    left_join(kegg_info %>% select(L3, pathway) %>% distinct(L3, .keep_all = TRUE), by = "L3") %>%
    filter(!is.na(pathway) & pathway != "") %>%
    group_by(pathway) %>%
    summarise(across(starts_with(c("RS", "S")), sum), .groups = "drop") %>%
    column_to_rownames(var = "pathway") %>%
    `rownames<-`(gsub(pattern = "\\[.*", replacement = "", x = rownames(.))) %>%
    filter(rownames(.) != "" & !is.na(rownames(.)))
  fin_pathway_df <- rbind(fin_pathway_df, unassign_df)
  
  # Superpathway Profile
  fin_superpathway_df <- raw_pathway_df %>%
    rownames_to_column(var = "pathway") %>%
    merge(pathway_super, by.x = "pathway", by.y = "V1", all.x = FALSE) %>%
    group_by(V2) %>%
    summarise(across(matches("^RS|^S"), sum), .groups = "drop") %>%
    column_to_rownames(var = "V2") %>%
    filter(rownames(.) %in% unique(kegg_info$L2)) %>%
    rownames_to_column(var = "L2") %>%
    left_join(kegg_info %>% select(L2, level2) %>% distinct(L2, .keep_all = TRUE), by = "L2") %>%
    group_by(level2) %>%
    summarise(across(starts_with(c("RS", "S")), sum), .groups = "drop") %>%
    column_to_rownames(var = "level2")
  fin_superpathway_df <- rbind(fin_superpathway_df, unassign_df)
  
  # Export Result
  write.csv(fin_pathway_df, paste0("./data/", kingdom, "_", type, "_tpm_pathway_abundance.csv"))
  write.csv(fin_superpathway_df, paste0("./data/", kingdom, "_", type, "_tpm_superpathway_abundance.csv"))
}
#-----------------------------------------------------------------

#Work Process
#-----------------------------------------------------------------
kingdoms <- c("Bacteria", "Fungi", "Protist")
types <- c("RS", "bulk")

for (kin in kingdoms) {
  for (type in types) {
    input_file <- paste0("../kegg_abundance/relative/", kin, "_", type, "_tpm_abundance.csv")
    input_df <- read.csv(input_file, check.names = F, row.names = 1)
    process_kegg_pathway(kingdom = kin, type = type, input_df = input_df)
  }
}
#-----------------------------------------------------------------