# 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# 加载必要的包
library(igraph)
library(tidyverse)
library(ggraph)
library(dplyr)

set.seed(123)

metadata <- read.csv("./data/metadata0425.csv")
kegg_name <- read.csv("./data/K_gene_name.csv",row.names = 1)
top10_orders <- c("Asparagales", "Arecales", "Fabales", "Rosales", "Lamiales",
                  "Malpighiales", "Sapindales", "Gentianales", "Malvales", "Myrtales")

lineage_colors <- c("1" = "#D53E4F", "2" = "#FFFFBF", "3" = "#3288BD")

microbe_types <- c("Bacteria", "Fungi", "Protists")

for (microbe_type in microbe_types) {
  
  # test
  # microbe_type <- "Fungi"
  
  # core 0.5
  # ko_file <- paste0("./data/", microbe_type, "_KO_core0.5.txt")
  # signal_file <- paste0("./data/", microbe_type, "_KO_phylogenetic_signal_qa(log).csv")
  
  #core 0.2
  # 二分网络是基于有无的，使用TPM数据不影响；系统发育信号检验使用的是拷贝数矫正的绝对丰度数据
  ko_file <- paste0("./data/", microbe_type, "_core0.2_tpm.csv")
  signal_file <- paste0("./data/", microbe_type, "_core0.2_qa_log_phylogenetic_signal.csv")
  
  KO_table <- read.csv(ko_file, check.names = F, row.names = 1)
  signal <- read.csv(signal_file)
  names(signal)[1] <- "KO"
  
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
    filter(abundance > 0.01) %>%  #### 过滤一下丰度
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
    dplyr::rename(from = KO, to = Order)
  
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
  nodes_edges <- nodes %>% left_join(edges, by = c("name" = "from"))
  
  # 只保留有显著性的KO节点
  nodes_signal_only <- nodes %>%
    filter(node_type == "Plant_Order" | (node_type == "KO" & has_signal))
  
  edges_signal_only <- edges %>%
    filter(from %in% nodes_signal_only$name & to %in% nodes_signal_only$name)
  
  # 输出有注释的edge文件
  edges_signal_only_note <- merge(edges_signal_only, nodes_signal_only, by.x = "from", by.y = "name", all.x = T)
  
  write.csv(edges_signal_only_note, paste0("./result/", microbe_type, "_edge_signal_only.csv"), row.names = F)
  write.csv(nodes_signal_only, paste0("./result/", microbe_type, "_node_signal_only.csv"), row.names = F)
  
  net_signal_only <- graph_from_data_frame(d = edges_signal_only, vertices = nodes_signal_only, directed = FALSE)
  
  p1 <- ggraph(net_signal_only, layout = 'fr') +
    geom_edge_link(alpha = 1, color = "grey") +
    geom_node_point(aes(
      size = ifelse(node_type == "KO", Relative_abundance, 3),
      shape = node_type
    ),
    color = nodes_signal_only$color_manual,
    stroke = 0.8,
    show.legend = FALSE) +
    scale_size_continuous(
      name = "logTPM",
      range = c(1, 3)
    ) +
    scale_shape_manual(
      name = "Node type",
      values = c("KO" = 16, "Plant_Order" = 1),
      labels = c("KO" = "KO (have signal)", 
                 "Plant_Order" = "Plant Order")
    ) +
    theme_void() +
    # labs(title = paste0(microbe_type, " - Signal only")) +
    # theme(legend.position = "right") +
    # 为植物节点添加标签，使用更大的字体
    geom_node_text(
      aes(label = ifelse(node_type == "Plant_Order", name, NA)),
      repel = TRUE, 
      size = 2,  # 增大植物标签字体大小
      # fontface = "bold",  # 加粗
      max.overlaps = Inf
    ) +
    # 为KO节点添加标签
    geom_node_text(
      aes(label = ifelse(number_of_host_lineages == "1" & node_type == "KO", gene, NA)),
      repel = TRUE, 
      size = 1,  # KO标签使用较小字体
      max.overlaps = Inf
    )
  
  width <- 9
  height <- 7
  ggsave(paste0("./bipartite_network/", microbe_type, '-KO-Order-network_less3_signal_only.png'), p1, width = width, height = height, dpi = 900, type = 'cairo', units = 'cm')
  ggsave(paste0("./bipartite_network/", microbe_type, '-KO-Order-network_less3_signal_only.pdf'), p1, width = width, height = height, units = 'cm')
  
  # p2 <- ggraph(net_signal_only, layout = 'fr') +
  #   geom_edge_link(alpha = 1, color = "grey") +
  #   geom_node_point(aes(
  #     size = ifelse(node_type == "KO", Relative_abundance, 4),
  #     shape = node_type
  #   ),
  #   color = nodes_signal_only$color_manual,
  #   stroke = 0.8,
  #   show.legend = TRUE) +
  #   scale_size_continuous(
  #     name = "logTPM",
  #     range = c(1, 6)
  #   ) +
  #   scale_shape_manual(
  #     name = "Node type",
  #     values = c("KO" = 16, "Plant_Order" = 1),
  #     labels = c("KO" = "KO (have signal)", 
  #                "Plant_Order" = "Plant Order")
  #   ) +
  #   theme_void() +
  #   labs(title = paste0(microbe_type, " - Signal only")) +
  #   theme(legend.position = "right") +
  #   geom_node_text(
  #     aes(label = ifelse(node_type == "Plant_Order", name, 
  #                        ifelse(number_of_host_lineages == "1", name, NA))),
  #     repel = TRUE, size = 2
  #   )
  # 
  # width <- 20.5
  # height <- 15.5
  # ggsave(paste0("./bipartite_network/", microbe_type, '-KO-Order-network_less3_withKO_signal_only.png'), p2, width = width, height = height, dpi = 900, type = 'cairo', units = 'cm')
  # ggsave(paste0("./bipartite_network/", microbe_type, '-KO-Order-network_less3_withKO_signal_only.pdf'), p2, width = width, height = height, units = 'cm')
  # 
  # 
  # p3 <- ggraph(net_signal_only, layout = "fr") +
  #   geom_edge_link(alpha = 1, color = "grey") +
  #   geom_node_point(aes(
  #     size = ifelse(node_type == "KO", Relative_abundance, 4),
  #     shape = node_type
  #   ),
  #   color = nodes_signal_only$color_manual,
  #   stroke = 0.8,
  #   show.legend = TRUE) +
  #   scale_size_continuous(
  #     name = "logTPM",
  #     range = c(1, 6)
  #   ) +
  #   scale_shape_manual(
  #     name = "Node type",
  #     values = c("KO" = 16, "Plant_Order" = 1),
  #     labels = c("KO" = "KO (have signal)", 
  #                "Plant_Order" = "Plant Order")
  #   ) +
  #   theme_void() +
  #   labs(title = paste0(microbe_type, " - Signal only")) +
  #   theme(legend.position = "right") +
  #   geom_node_text(
  #     aes(label = ifelse(node_type == "Plant_Order", name, 
  #                        ifelse(number_of_host_lineages == "1", gene, NA))),
  #     repel = TRUE, size = 2
  #   )
  # 
  # width <- 20.5
  # height <- 15.5
  # ggsave(paste0("./bipartite_network/", microbe_type, '-KO-Order-network_less3_withgene_signal_only.png'), p3, width = width, height = height, dpi = 900, type = 'cairo', units = 'cm')
  # ggsave(paste0("./bipartite_network/", microbe_type, '-KO-Order-network_less3_withgene_signal_only.pdf'), p3, width = width, height = height, units = 'cm')
}