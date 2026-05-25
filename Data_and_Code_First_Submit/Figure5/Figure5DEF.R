pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(igraph)
library(tidyverse)
library(ggraph)
library(dplyr)
library(patchwork)  # 加载patchwork包
set.seed(123)

metadata <- read.csv("./data/metadata0425.csv")
kegg_name <- read.csv("./data/K_gene_name.csv", row.names = 1)
top10_orders <- c("Asparagales", "Arecales", "Fabales", "Rosales", "Lamiales",
                  "Malpighiales", "Sapindales", "Gentianales", "Malvales", "Myrtales")

lineage_colors <- c("1" = "#D53E4F", "2" = "#FFFFBF", "3" = "#3288BD")

microbe_types <- c("Bacteria", "Fungi", "Protist")

# 创建空列表存储三种类型的图
p1_list <- list()
p2_list <- list()
p3_list <- list()

for (microbe_type in microbe_types) {
  
  ko_file <- paste0("./data/", microbe_type, "_KO_core0.5.txt")
  signal_file <- paste0("./data/", microbe_type, "_KO_phylogenetic_signal_qa(log).csv")
  
  KO_table <- read.delim(ko_file, check.names = F, row.names = 1)
  signal <- read.csv(signal_file)
  
  signal_KO <- signal %>% 
    filter(P < 0.05) %>% 
    select(KO)
  
  KO_table_long <- KO_table %>%
    rownames_to_column(var = "KO") %>%
    pivot_longer(
      cols = -KO,
      names_to = "TreeID",
      values_to = "abundance",
      values_drop_na = TRUE
    ) %>%
    filter(abundance > 0) %>% 
    merge(., metadata %>% select(TreeID, Order), by = 'TreeID')
  
  KO_table_long_top_order <- KO_table_long %>%
    filter(Order %in% top10_orders, abundance >= 1)
  
  KO_table_long_top_order_sum <- KO_table_long_top_order %>%
    group_by(KO) %>%
    summarise(
      n_order = n_distinct(Order),
      abundance_sum = sum(abundance),
      .groups = "drop"
    ) %>%
    mutate(log_abundance_sum = log10(abundance_sum + 1))
  
  filtered_ko_ids <- KO_table_long_top_order_sum %>%
    filter(n_order <= 3) %>%
    pull(KO)
  
  KO_table_long_top_order_filtered <- KO_table_long_top_order %>%
    filter(KO %in% filtered_ko_ids)
  
  KO_table_long_top_order_sum_filtered <- KO_table_long_top_order_sum %>%
    filter(KO %in% filtered_ko_ids)
  
  edges <- KO_table_long_top_order_filtered %>%
    select(KO, Order) %>%
    distinct() %>%
    rename(from = KO, to = Order)
  
  nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
    left_join(KO_table_long_top_order_sum_filtered, by = c("name" = "KO")) %>%
    mutate(
      node_type = ifelse(name %in% edges$to, "Plant_Order", "KO"),
      Relative_abundance = log_abundance_sum,
      number_of_host_lineages = n_order,
      has_signal = ifelse(name %in% signal_KO$KO & node_type == "KO", TRUE, FALSE)
    ) %>%
    mutate(
      number_of_host_lineages = as.character(number_of_host_lineages),
      color_manual = ifelse(
        node_type == "KO",
        lineage_colors[number_of_host_lineages],
        "black"
      )
    ) %>%
    select(name, node_type, Relative_abundance, number_of_host_lineages, color_manual, has_signal)
  nodes <- nodes %>% left_join(kegg_name, by = c("name" = "KO"))
  write.csv(edges, paste0("./result/", microbe_type, "_edge.csv"), row.names = F)
  write.csv(nodes, paste0("./result/", microbe_type, "_node.csv"), row.names = F)
  
  net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # p1: 只显示植物目名称
  p1 <- ggraph(net, layout = "fr") +
    geom_edge_link(alpha = 1, color = "grey") +
    geom_node_point(aes(
      size = ifelse(node_type == "KO", Relative_abundance, 4),
      shape = ifelse(node_type == "KO",
                     ifelse(has_signal, "has_signal", "no_signal"),
                     "plant")
    ),
    color = nodes$color_manual,
    stroke = 0.8,
    show.legend = TRUE) +
    scale_size_continuous(
      name = "Relative abundance(log)",
      range = c(1, 6)
    ) +
    scale_shape_manual(
      name = "Node type",
      values = c("has_signal" = 16, "no_signal" = 1, "plant" = 1),
      labels = c("has_signal" = "KO (have signal)", 
                 "no_signal" = "KO (no signal)", 
                 "plant" = "Plant Order")
    ) +
    theme_void() +
    labs(title = microbe_type) +
    theme(legend.position = "none") +
    geom_node_text(
      aes(label = ifelse(node_type == "Plant_Order", name, NA)),
      repel = TRUE, size = 2
    )
  
  # p2: 显示植物目和有信号的KO名称
  p2 <- ggraph(net, layout = "fr") +
    geom_edge_link(alpha = 1, color = "grey") +
    geom_node_point(aes(
      size = ifelse(node_type == "KO", Relative_abundance, 4),
      shape = ifelse(node_type == "KO",
                     ifelse(has_signal, "has_signal", "no_signal"),
                     "plant")
    ),
    color = nodes$color_manual,
    stroke = 0.8,
    show.legend = TRUE) +
    scale_size_continuous(
      name = "Relative abundance(log)",
      range = c(1, 6)
    ) +
    scale_shape_manual(
      name = "Node type",
      values = c("has_signal" = 16, "no_signal" = 1, "plant" = 1),
      labels = c("has_signal" = "KO (have signal)", 
                 "no_signal" = "KO (no signal)", 
                 "plant" = "Plant Order")
    ) +
    theme_void() +
    labs(title = microbe_type) +
    theme(legend.position = "none") +
    geom_node_text(
      aes(label = ifelse(number_of_host_lineages == "1" & has_signal, name, 
                         ifelse(node_type == "Plant_Order", name, NA))),
      repel = TRUE, size = 2
    )
  
  # p3: 显示植物目和有信号KO的基因名称
  p3 <- ggraph(net, layout = "fr") +
    geom_edge_link(alpha = 1, color = "grey") +
    geom_node_point(aes(
      size = ifelse(node_type == "KO", Relative_abundance, 4),
      shape = ifelse(node_type == "KO",
                     ifelse(has_signal, "has_signal", "no_signal"),
                     "plant")
    ),
    color = nodes$color_manual,
    stroke = 0.8,
    show.legend = TRUE) +
    scale_size_continuous(
      name = "Relative abundance(log)",
      range = c(1, 6)
    ) +
    scale_shape_manual(
      name = "Node type",
      values = c("has_signal" = 16, "no_signal" = 1, "plant" = 1),
      labels = c("has_signal" = "KO (have signal)", 
                 "no_signal" = "KO (no signal)", 
                 "plant" = "Plant Order")
    ) +
    theme_void() +
    labs(title = microbe_type) +
    theme(legend.position = "none") +
    geom_node_text(
      aes(label = ifelse(number_of_host_lineages == "1" & has_signal, gene, 
                         ifelse(node_type == "Plant_Order", name, NA))),
      repel = TRUE, size = 2
    )
  
  # 将当前微生物的图添加到对应列表中
  p1_list[[microbe_type]] <- p1
  p2_list[[microbe_type]] <- p2
  p3_list[[microbe_type]] <- p3
}

# 使用patchwork纵向拼接三种微生物的p1图
combined_p1 <- wrap_plots(p1_list, ncol = 1) + 
  plot_annotation(tag_levels = 'A')

# 使用patchwork纵向拼接三种微生物的p2图
combined_p2 <- wrap_plots(p2_list, ncol = 1) + 
  plot_annotation(tag_levels = 'A')

# 使用patchwork纵向拼接三种微生物的p3图
combined_p3 <- wrap_plots(p3_list, ncol = 1) + 
  plot_annotation(tag_levels = 'A')


ggsave("Figure5BDE.jpg", combined_p1, width = 10, height = 30, dpi = 900, units = "cm")
ggsave("Figure5BDE.pdf", combined_p1, width = 10, height = 30, dpi = 900,units = "cm")
ggsave("Figure5BDE.tiff", combined_p1, width = 10, height = 30, dpi = 900,units = "cm")

ggsave("Figure5BDE_KO.jpg", combined_p2, width = 10, height = 30, dpi = 900, units = "cm")
ggsave("Figure5BDE_KO.pdf", combined_p2, width = 10, height = 30, dpi = 900,units = "cm")
ggsave("Figure5BDE_KO.tiff", combined_p2, width = 10, height = 30, dpi = 900,units = "cm")

ggsave("Figure5BDE_gene.jpg", combined_p3, width = 10, height = 30, dpi = 900, units = "cm")
ggsave("Figure5BDE_gene.pdf", combined_p3, width = 10, height = 30, dpi = 900,units = "cm")
ggsave("Figure5BDE_gene.tiff", combined_p3, width = 10, height = 30, dpi = 900,units = "cm")