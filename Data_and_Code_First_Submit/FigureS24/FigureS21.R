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

kingdom <- c("bac","fun","pro")
subtitle_map <- c(
  "bac" = "Bacteria",
  "fun" = "Fungi",
  "pro" = "Protist"
)
plot_list <- list()
for (kin in kingdom){
  data <- read.csv(paste0("./data/",kin,"_qa_log_metastat_res_diff_ASV_matrix(p0.05q0.2).csv"), row.names = 1)
  data_matrix <- as.matrix(data)
  n <- nrow(data_matrix)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i > j) {
        data_matrix[i, j] <- -data_matrix[i, j]
      }
    }
  }
  data_long <- melt(data_matrix)
  
  data_long$row <- as.numeric(factor(data_long$Var1, levels = rownames(data)))
  data_long$col <- as.numeric(factor(data_long$Var2, levels = colnames(data)))
  
  data_long$transformed_value <- with(data_long, ifelse(
    row < col, 
    log1p(value),
    ifelse(
      row > col, 
      -log1p(-value), 
      NA  
    )
  ))
  custom_subtitle <- subtitle_map[[kin]]
  p <- ggplot(data_long, aes(x = Var2, y = Var1, fill = transformed_value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(
      low = "#5BA2D8", 
      mid = "white", 
      high = "#EC499A", 
      midpoint = 0,
      na.value = "gray",
      name = "log(number)",
      breaks = c(min(data_long$transformed_value,na.rm = T),max(data_long$transformed_value,na.rm = T)),
      labels = c("Depleted", "Enriched"),
      limits = c(min(data_long$transformed_value,na.rm = T),max(data_long$transformed_value,na.rm = T))
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      plot.margin = margin(t = 3, r = 0, b = 0, l = 0, unit = "pt"),
      axis.text.x = element_text(size = 7,angle = 45, hjust = 1, vjust = 1,color = "black"),
      axis.text.y = element_text(size = 7,hjust = 1, vjust = 0.5,color = "black"),
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid = element_blank(),
      legend.key.width = unit(4, "mm"), 
      legend.key.height = unit(4, "mm"), 
    ) +
    labs(title = "", x = "", y = custom_subtitle) +
    geom_text(
      aes(label = abs(value)),
      color = "black",
      size = 2,
      check_overlap = TRUE
    ) +
    scale_x_discrete(limits = colnames(data)) +
    scale_y_discrete(limits = rev(rownames(data))) +
    geom_tile(
      data = subset(data_long, Var1 == Var2),
      fill = "gray",
      color = "black",
      show.legend = FALSE
    )
  plot_list[[kin]] <- p
}
final_plot1 <- grid.arrange(
  grobs = plot_list,
  heights = rep(1, 3),
  ncol = 1,
  padding = unit(0, "line")
)

microbes <- c("bac", "fun", "pro")

target_columns <-c("Rosales - Fabales","Lamiales - Fabales","Malpighiales - Fabales","Sapindales - Fabales","Fabales - Gentianales","Fabales - Asparagales","Malvales - Fabales","Myrtales - Fabales","Arecales - Fabales","Others - Fabales")

new_columns <- c("Fabales vs Rosales", "Fabales vs Lamiales", "Fabales vs Malpighiales", "Fabales vs Sapindales", 
                 "Fabales vs Gentianales", "Fabales vs Asparagales", "Fabales vs Malvales", "Fabales vs Myrtales", 
                 "Fabales vs Arecales", "Fabales vs Others")

microbe_colors <- setNames(
  c("#8DD3C7","#BEBADA","#FFFFB3"),
  microbes
)

combined_data <- map_dfr(microbes, function(microbe) {
  read.csv(paste0("./data/", microbe, "_qa_log_tree_metastat_res_diff_ASV_enrich(p0.05q0.2).csv")) %>%
    filter(FDR < 0.05) %>%
    filter(str_detect(GroupComparison, "Fabales")) %>%
    distinct(Taxon, GroupComparison, .keep_all = TRUE) %>%
    mutate(EnrichmentFactor = if_else(Direction == "Depleted", 
                                      -EnrichmentFactor, 
                                      EnrichmentFactor)) %>%
    mutate(Kingdom = microbe)
}) 

combined_data <- combined_data %>% filter(!Taxon %in% c("AD3", "0319-6G20", "JG36-TzT-191", "Unassign","Subgroup_2","MB-A2-108"))

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
  show_legend = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm"),
    ncol = 1
  )
)

final_plot2 <- Heatmap(
  matrix = full_matrix,
  col = colorRamp2(c(min(full_matrix), 0, max(full_matrix)), c("#5BA2D8", "white", "#EC499A")),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  #  column_title = " ",
  column_title_side = "top",  
  column_title_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 0.2),
  
  left_annotation = row_ha,
  row_title = NULL,
  row_split = factor(rep(seq_along(rle(row_annotation$Kingdom)$lengths), rle(row_annotation$Kingdom)$lengths)),
  row_gap = unit(2, "mm"),
  width = ncol(full_matrix) * unit(5.5, "mm"),
  #  height = nrow(full_matrix) * unit(4, "mm"),
  
  heatmap_legend_param = list(
    title = "EnrichmentFactor",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    gap = unit(5, "cm") 
  )
)

final_plot2_grob <- grid.grabExpr(draw(final_plot2))

combined_plot <- grid.arrange(
  grobs = list(final_plot1, final_plot2_grob),
  ncol = 2,
  heights = unit(1, "null"),
  padding = unit(0, "line")
)

ggsave("FigureS21.pdf", combined_plot, width = 21, height = 21,unit = "cm",dpi = 900,bg="white")
ggsave("FigureS21.jpg", combined_plot, width = 21, height = 21,unit = "cm",dpi = 900,bg="white")
ggsave("FigureS21.tiff", combined_plot, width = 21, height = 21,unit = "cm",dpi = 900,bg="white")
