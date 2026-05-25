# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(dplyr)
library(igraph)
library(ggraph)
library(tibble)
library(tidyr)
library(ggplot2)


ko_table <- read.csv("./data/kegg_pathway_for_1KPM.csv")
ko_table$pathway <- sub("\\s*\\[.*", "", ko_table$pathway)
ko_table <- ko_table %>% filter(L3 !="map00942")
top_orders <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales")
pathway_abundance <- read.csv("./data/All_pathway_order_abundance.csv",row.names = 1,check.names = F)
rownames(pathway_abundance) <- trimws(rownames(pathway_abundance), which = "right")
level2_abundance <- read.csv("./data/All_level2_order_abundance.csv",row.names = 1,check.names = F)
level1_abundance <- read.csv("./data/All_level1_order_abundance.csv",row.names = 1,check.names = F)
ko_table <- ko_table %>% filter(pathway %in% rownames(pathway_abundance))

edges <- bind_rows(
  data.frame(from = "root", to = unique(ko_table$level1)),
  ko_table %>% distinct(level1, level2) %>% transmute(from = level1, to = level2),
  ko_table %>% distinct(level2, pathway) %>% transmute(from = level2, to = pathway)
) %>% 
  slice(-101)

nodes <- data.frame(
  name = unique(c(edges$from, edges$to))
) %>% 
  mutate(
    level = case_when(
      name == "root" ~ "root",
      name %in% ko_table$level1 ~ "level1",
      name %in% ko_table$level2 ~ "level2",
      TRUE ~ "pathway"
    )
  )
mygraph <- graph_from_data_frame(edges, vertices = nodes)

calculate_total_abundance <- function(nodes, abundance_table, level_type) {
  nodes %>%
    filter(level == level_type) %>%
    mutate(
      total_abundance = sapply(name, function(x) {
        if (x %in% rownames(abundance_table)) {
          rowSums(abundance_table[x, , drop = FALSE], na.rm = TRUE)
        } else {
          0
        }
      })
    )
}
tree_layout <- create_layout(mygraph, 'dendrogram', circular = TRUE)

nodes_abundance <- bind_rows(
  pathway_abundance %>% 
    mutate(total = rowSums(.), level = "pathway") %>% 
    tibble::rownames_to_column("name"),
  
  level2_abundance %>% 
    mutate(total = rowSums(.), level = "level2") %>% 
    tibble::rownames_to_column("name"),
  
  level1_abundance %>% 
    mutate(total = rowSums(.), level = "level1") %>% 
    tibble::rownames_to_column("name")
) %>%
  group_by(name) %>% 
  summarise(total_abundance = max(total, na.rm = TRUE))

tree_data <- tree_layout %>% 
  left_join(nodes_abundance, by = "name") %>% 
  mutate(total_abundance = ifelse(is.na(total_abundance), 0, total_abundance))
order_levels <- c(top_orders, "Others")

column_sums <- colSums(pathway_abundance[, order_levels])

result_df <- data.frame(
  Order =  names(pathway_abundance[, order_levels]),
  Sum = column_sums/sum(column_sums)
)

pathway_percent <- pathway_abundance %>%
  rownames_to_column("pathway") %>%
  rowwise() %>%
  mutate(across(-pathway, ~ {
    if (sum(c_across(-pathway)) == 0) . 
    else . / sum(c_across(-pathway))
  })) %>%
  column_to_rownames("pathway")
pathway_percent <- pathway_percent[, order_levels]

level2_percent <- level2_abundance %>%
  rownames_to_column("level2") %>%
  rowwise() %>%
  mutate(across(-level2, ~ {
    if (sum(c_across(-level2)) == 0) . 
    else . / sum(c_across(-level2))
  })) %>%
  column_to_rownames("level2")
level2_percent <- level2_percent[, order_levels]




bar_data_pathway <- pathway_percent %>%
  rownames_to_column(var = "pathway_name") %>%
  inner_join(
    tree_layout %>% filter(level == "pathway") %>% select(name, x, y),
    by = c("pathway_name" = "name")
  ) %>%
  pivot_longer(
    cols = -c(pathway_name, x, y),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(category = factor(category, levels = order_levels)) %>%
  group_by(pathway_name) %>%
  mutate(
    angle = atan2(y, x),
    cumulative = cumsum(value),
    start = lag(cumulative, default = 0) * 0.3,
    end = cumulative * 0.3
  ) %>%
  arrange(pathway_name, category) %>%
  mutate(
    start_x = x + start * cos(angle),
    start_y = y + start * sin(angle),
    end_x = x + end * cos(angle),
    end_y = y + end * sin(angle)
  )

bar_data_level2 <- level2_percent %>%
  rownames_to_column(var = "level2_name") %>%
  inner_join(
    tree_layout %>% filter(level == "level2") %>% select(name, x, y),
    by = c("level2_name" = "name")
  ) %>%
  pivot_longer(
    cols = -c(level2_name, x, y),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(category = factor(category, levels = order_levels)) %>%
  group_by(level2_name) %>%
  mutate(
    angle = atan2(y, x),
    cumulative = cumsum(value),
    start = lag(cumulative, default = 0) * 0.3,
    end = cumulative * 0.3
  ) %>%
  arrange(level2_name, category) %>%
  mutate(
    start_x = x + start * cos(angle),
    start_y = y + start * sin(angle),
    end_x = x + end * cos(angle),
    end_y = y + end * sin(angle)
  )

bar_data_combined <- bind_rows(
  bar_data_pathway %>% mutate(level_type = "pathway"),
  bar_data_level2 %>% mutate(level_type = "level2")
)
order_colors <- c("#E41A1C","#596A98","#449B75","#6B886D","#AC5782","#FF7F00","#FFE528","#C9992C","#C66764","#E485B7","#999999")
#order_colors <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#D9D9D9")
level_colors <- c("level1" = "#66C2A5", "level2" = "#FC8D62", "pathway" = "#8DA0CB")

p <- ggraph(tree_data, layout = 'dendrogram', circular = TRUE) +
  geom_edge_diagonal(
    aes(color = node2.level),
    alpha = 0.6,
    width = 0.6,
    show.legend = FALSE 
  ) +
  
  geom_segment(
    data = bar_data_combined %>% filter(level_type == "level2"),
    aes(x = start_x, y = start_y,
        xend = end_x, yend = end_y,
        color = category),
    linewidth = 4,
    alpha = 1,
    lineend = "butt"
  ) +
  
  geom_segment(
    data = bar_data_combined %>% filter(level_type == "pathway"),
    aes(x = start_x, y = start_y,
        xend = end_x, yend = end_y,
        color = category),
    linewidth = 1.5,
    lineend = "butt"
  ) +
  
  geom_node_point(
    aes(size = total_abundance, fill = level),
    shape = 21,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    name = "Plant Order",
    values = order_colors,
    guide = guide_legend(
      ncol = 1,
      override.aes = list(size = 4),
      title.position = "top"
    )
  ) +
  
  scale_fill_manual(values = level_colors, guide = "none") +
  scale_edge_color_manual(values = level_colors, guide = "none") +
  scale_size(guide = "none") +
  
  geom_node_text(
    aes(label = name, filter = level %in% c("level1", "level2")),
    size = 1,
    hjust = 0.5,
    vjust = -1,
    check_overlap = TRUE
  ) +
  
  labs(title = "KEGG Annotation Summary") +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14,
      margin = margin(b = 10)
    ),
    legend.margin = margin(t = 10)
  )
width <- 30
height <- 30
name <- paste0(dir_name,"/Figure1B-kegg_tree_summary")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")