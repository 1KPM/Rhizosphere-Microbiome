pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


library(reshape2)
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(purrr)
library(stringr)
library(tibble)
library(dplyr)
library(circlize)
library(tidyr)

##Enrich Data Prepare
source("enrichment_functions.R")
categories <- c("Bacteria" = "bac", "Fungi" = "fun", "Protists" = "pro")

for (category_name in names(categories)) {
  category_abbr <- categories[category_name]
  
  data_diff <- read.csv(paste0("./stats/", category_abbr, "_ko_qa_log_tree_metastat_res_diff_KO_core0.2.csv"), row.names = 1)  %>% 
    mutate(Significant = (pvalue < 0.05 & qvalue < 0.2)) %>% 
    filter(Significant)
  
  data_diff <- data_diff %>% mutate(KO = str_extract(Taxa, "^[^|]+"))
  
  taxonomy <- read.csv("./data/kegg_pathway_for_1KPM.csv")
  ko_pathway <- taxonomy %>% 
    select(KO, pathway) %>% 
    distinct(KO, .keep_all = TRUE)
  ko_pathway$pathway <- sub("\\[.*", "", ko_pathway$pathway)
  
  data_Fabales <- data_diff %>% 
    left_join(ko_pathway, by = "KO") %>% 
    mutate(pathway = replace_na(pathway, "Others")) %>%
    rename(FeatureID = KO, Order = pathway) %>%
    filter(grepl("Fabales", Comparison)) %>%
    filter(Order != "" & Order != "Others")
  
  abundance <- read.csv(
    file = paste0("../../kegg_abundance/absolute/", category_name, "_RS_quantitative_abundance.csv"), 
    row.names = 1, 
    check.names = FALSE
  )
  
  taxonomy_filtered <- taxonomy %>%
    filter(KO %in% rownames(abundance)) %>%
    select(KO, pathway)
  taxonomy_filtered$pathway <- sub("\\[.*", "", taxonomy_filtered$pathway)
  colnames(taxonomy_filtered) <- c("FeatureID", "Order")
  
  comparison_groups <- unique(data_Fabales$Comparison)
  
  enrich_results <- list()
  
  for (group in comparison_groups) {
    enrich_data <- data_Fabales %>% 
      filter(Comparison == group & Group == "Fabales")
    enrich_res <- perform_taxa_enrichment(
      tax_level = "Order",
      diff_data = enrich_data,
      tax_data = taxonomy_filtered,
      min_counts = 1
    )
    if (nrow(enrich_res) > 0) {
      enrich_res$GroupComparison <- group
      enrich_res$Direction <- "Enriched"
      enrich_results[[paste0(group, "_Enriched")]] <- enrich_res
    }
    
    deplete_data <- data_Fabales %>% 
      filter(Comparison == group & Group != "Fabales")
    deplete_res <- perform_taxa_enrichment(
      tax_level = "Order",
      diff_data = deplete_data,
      tax_data = taxonomy_filtered,
      min_counts = 1
    )
    if (nrow(deplete_res) > 0) {
      deplete_res$GroupComparison <- group
      deplete_res$Direction <- "Depleted"
      enrich_results[[paste0(group, "_Deplete")]] <- deplete_res
    }
  }
    final_result <- do.call(rbind, enrich_results)
    write.csv(final_result,paste0("./stats/", category_abbr, "_ko_qa_log_tree_metastat_res_diff_KO_enrich(core0.2p0.05q0.2).csv"))
}
  


##Enrich Plot
microbes <- c("bac", "fun", "pro")

target_columns <-c("Rosales - Fabales","Lamiales - Fabales","Malpighiales - Fabales","Sapindales - Fabales",
                   "Fabales - Gentianales","Fabales - Asparagales","Malvales - Fabales","Myrtales - Fabales",
                   "Arecales - Fabales","Others - Fabales")

new_columns <- c("Fabales vs Rosales", "Fabales vs Lamiales", "Fabales vs Malpighiales", "Fabales vs Sapindales", 
                 "Fabales vs Gentianales", "Fabales vs Asparagales", "Fabales vs Malvales", "Fabales vs Myrtales", 
                 "Fabales vs Arecales", "Fabales vs Others")

microbe_colors <- setNames(
  c("#8DD3C7","#BEBADA","#FFFFB3"),
  microbes
)

combined_data <- map_dfr(microbes, function(microbe) {
  read.csv(paste0("./stats/", microbe, "_ko_qa_log_tree_metastat_res_diff_KO_enrich(core0.2p0.05q0.2).csv")) %>%
    filter(FDR < 0.05) %>%
    filter(str_detect(GroupComparison, "Fabales")) %>%
    distinct(Taxon, GroupComparison, .keep_all = TRUE) %>%
    mutate(EnrichmentFactor = if_else(Direction == "Depleted", 
                                      -EnrichmentFactor, 
                                      EnrichmentFactor)) %>%
    mutate(Kingdom = microbe)
}) 

combined_data <- combined_data %>% filter(!Taxon %in% c("Ribosome biogenesis in eukaryotes"))

combined_data <- combined_data %>%
  mutate(
    Taxon = ifelse(Kingdom == "bac", trimws(Taxon, "right"), Taxon)
  )

full_matrix <- combined_data %>%
  pivot_wider(
    id_cols = c(Taxon, Kingdom),
    names_from = GroupComparison,
    values_from = EnrichmentFactor,
    values_fill = 0
  ) %>%
  arrange(factor(Kingdom, levels = microbes)) %>% 
  column_to_rownames("Taxon") %>%
  select(-Kingdom) %>%
  as.matrix()

new_matrix <- matrix(0, 
                     nrow = nrow(full_matrix),
                     ncol = length(target_columns),
                     dimnames = list(rownames(full_matrix), target_columns))

common_cols <- intersect(colnames(full_matrix), target_columns)
for (col in common_cols) {
  new_matrix[, col] <- full_matrix[, col]
}

full_matrix <- as.matrix(new_matrix)
colnames(full_matrix) <- new_columns

row_annotation <- combined_data %>%
  distinct(Taxon, Kingdom) %>%
  arrange(factor(Kingdom, levels = microbes)) %>%
  column_to_rownames("Taxon") %>%
  select(Kingdom)

row_ha <- rowAnnotation(
  df = row_annotation["Kingdom"],
  col = list(Kingdom = microbe_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm"),
    ncol = 1
  )
)

final_plot <- Heatmap(
  matrix = full_matrix,
  col = colorRamp2(
    breaks = c(-15, 0, 15),
    c("#5BA2D8", "white", "#EC499A")
  ),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  column_title_side = "top",  
  column_title_gp = gpar(fontsize = 6),
  rect_gp = gpar(col = "black", lwd = 0.2),
  
  left_annotation = row_ha,
  row_title = NULL,
  row_split = factor(rep(seq_along(rle(row_annotation$Kingdom)$lengths), rle(row_annotation$Kingdom)$lengths)),
  row_gap = unit(2, "mm"),
  width = ncol(full_matrix) * unit(6, "mm"),
  #  height = nrow(full_matrix) * unit(4, "mm"),
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    title = "EnrichmentFactor",
    at = seq(-15, 15, by = 5),
    labels = seq(-15, 15, by = 5),
    title_gp = gpar(fontsize = 6),
    labels_gp = gpar(fontsize = 6)
  )
)


final_plot_grob <- grid.grabExpr(draw(final_plot))

ggsave("./plots/diff_ko_enrich_pathway.pdf", final_plot_grob, width = 12, height = 24, dpi = 300,units = "cm")
ggsave("./plots/diff_ko_enrich_pathway.png", final_plot_grob, width = 12, height = 24, dpi = 300,units = "cm")
