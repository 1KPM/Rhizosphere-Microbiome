pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(ggvenn)
library(tidyverse)
library(RColorBrewer)
library(patchwork)

##Figure4A
order_list <- c("Fabales", "Rosales", "Malpighiales", "Gentianales")
hub_lists <- list()
for (o in order_list) {
  hub_info <- read.csv(paste0('./data/Figure4-', o, '_Cross_kingdom_network_hub_info.csv'), row.names = 1)
  hub_list <- hub_info[hub_info$roles != 'Peripherals',]$label
  hub_lists[[o]] <- hub_list
}
common_hubs <- Reduce(intersect, hub_lists)

p1 <- ggvenn::ggvenn(
    hub_lists,
    columns = NULL,
    show_elements = FALSE,
    label_sep = "\n",
    show_percentage = FALSE,
    digits = 1,
    fill_color = c("#E41A1C", "#596A98", "#6B886D", "#FF7F00"),
    fill_alpha = 0.5,
    stroke_color = "white",
    stroke_alpha = 0.5,
    stroke_size = 0.5,
    stroke_linetype = "solid",
    set_name_color = "black",
    text_color = "black",
    set_name_size = 6 / 2.835,
    text_size = 5 / 2.835
)
p1

##Figure4A2
all_hub_info <- read.csv(paste0('./data/Figure4-Cross_kingdom_network_hub_info.csv'), row.names = 1)
all_hub_list <- all_hub_info[all_hub_info$roles != 'Peripherals',]$label

n1 <- length(all_hub_list)
n2 <- length(common_hubs)
n12 <- length(intersect(all_hub_list, common_hubs))
n1_not2 <- n1 - n12
n2_not1 <- n2 -n12
n_not1_not2 <- length(union(row.names(all_hub_info), row.names(common_hubs)))
tmp_df <- data.frame(
  'SubHub' = c(n12, n2_not1), 'Non_subHub' = c(n1_not2, n_not1_not2),
  row.names = c('AllHub', 'Non_AllHub')
)
sig_res <- chisq.test(tmp_df)

tmp_list <- list()
tmp_list[['All hub']] <- all_hub_list
tmp_list[['Stable hub']] <- common_hubs

p2 <- ggvenn(tmp_list,
             columns = NULL,
             show_elements = FALSE,
             label_sep = "\n",
             show_percentage = TRUE,
             digits = 1,
             fill_color = c("#6EB9C3", "#e17e16"),
             fill_alpha = 0.8,
             stroke_color = "white",
             stroke_alpha = 0.5,
             stroke_size = 0.5,
             stroke_linetype = "solid",
             set_name_color = "black",
             text_color = "black",
             set_name_size = 6 / 2.835,
             text_size = 5 / 2.835
) +
    labs(caption = bquote(
        chi^2 ~ "test:" ~ .(ifelse(sig_res$p.value < 0.001, 'p < 0.001',
                                   ifelse(sig_res$p.value < 0.01, 'p < 0.01',
                                          ifelse(sig_res$p.value < 0.05, 'p < 0.05', 'p > 0.05'))))
    )) + 
    theme(plot.title = element_blank(),
          plot.caption = element_text(
              size = 6, 
              color = 'black',
              hjust = 0.5
          ),
          panel.spacing = unit(0, "cm"),
          plot.margin = margin(2, 0, 0, 0, "cm")
    )
p2

##Figure4B
all_hub <- read.csv("./data/Figure4-Cross_kingdom_network_hub_info.csv")
inter_hub <- intersect(common_hubs,all_hub$label)
KO_list <- substr(inter_hub, 2, nchar(inter_hub))

result_all<-read.csv("./data/Figure4-34_ko_enrichment_results.csv")
kegg<-read.csv("./data/K_gene_name.csv")
result_all <- result_all %>%
  arrange(desc(Count)) %>%
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    hub = factor(hub)
  )
result_all_processed <- result_all %>%
  mutate(KO_num = str_sub(hub, start = 2)) %>% 
  left_join(kegg %>% select(KO, kegg_desc = gene),
            by = c("KO_num" = "KO")) %>%
  mutate(
    hub_label = paste0(hub, "|", coalesce(kegg_desc, "Unannotated")),
    hub = factor(hub, levels = unique(hub)),
    hub_label_clean = str_split(hub_label, fixed("["), n = 2, simplify = TRUE)[, 1]
  ) 

result_all_processed <- result_all_processed %>%
  add_row(
    Description = "None",
    hub_label_clean = "fK01064|WNT9",
    Count = 0,
    qvalue = 1
  )

mat_count <- result_all_processed %>%
  select(Description, hub_label_clean, Count) %>%
  group_by(Description, hub_label_clean) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  complete(Description, hub_label_clean, fill = list(Count = 0)) %>%
  pivot_wider(
    names_from = hub_label_clean,
    values_from = Count
  ) %>%
  column_to_rownames(var = "Description") %>%
  as.matrix() %>%
  {storage.mode(.) <- "numeric"; .}

mat_count <- mat_count[rownames(mat_count) != "None", , drop = FALSE]

row_clust <- hclust(dist(mat_count), method = "average")
row_order <- rownames(mat_count)[row_clust$order]

col_clust <- hclust(dist(t(mat_count)), method = "single")
col_order <- colnames(mat_count)[col_clust$order]
col_order[c(9, 10, 33, 34)] <- col_order[c(33, 34, 9, 10)]

result_all_processed <- result_all_processed %>%
  mutate(
    Description = factor(Description) %>% 
      fct_relevel(row_order) %>% 
      fct_drop(),
    hub_label_clean = factor(hub_label_clean) %>% 
      fct_relevel(col_order) %>% 
      fct_drop()
  )


p3 <- ggplot(result_all_processed,aes(x = Description, y = hub_label_clean)) +
    geom_point(
        aes(size = Count),
        position = position_nudge(x = 0.08, y = -0.08),
        color = "gray20",
        alpha = 0.4,
        shape = 19,
        show.legend = FALSE
    ) +
    geom_point(
        aes(size = Count, fill = -log10(qvalue)),
        color = "black",
        shape = 21,
        stroke = 0.3,
        alpha = 1 
    )+
    scale_x_discrete(
        limits = row_order,
        expand = expansion(add = 1)
    ) +
    scale_y_discrete(
        limits = col_order,
        expand = expansion(add = 1),
        drop = FALSE
    ) +
    scale_fill_gradientn(
        name = expression(-log[10](q[value])),
        colors = c("#FFFFE0", "#E41A1C"),
        values = scales::rescale(c(0, 5, 10, 20, 50)),
        limits = c(0, 50)
    ) +
    scale_size_continuous(
        name = expression(KO_count),
        range = c(0.3, 3),
        limits = c(0, 50),
        breaks = scales::extended_breaks(n = 5),
        guide = guide_legend(order = 1)
    ) +
    scale_color_gradientn(
        name = expression(-log[10](q[value])),
        colors = viridis::viridis(10),
        values = scales::rescale(c(0, 5, 10, 20, 50)),
        limits = c(0, 50),
        guide = guide_colorbar(
            barwidth = unit(0.6, "cm"),
            barheight = unit(4, "cm"),
            frame.colour = "black",
            ticks.colour = "black"
        )
    )+
    labs(
        x = NULL,
        y = NULL
    ) +
theme_bw() + theme(
    text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 7, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
)+
    #  coord_equal(ratio = 0.8) +
    ggnewscale::new_scale("size") 
p3



##Figure4C
final_result <- read.csv(file="./data/Figure4-KO_interaction_distribution.csv")
type_counts <- final_result %>%
  group_by(correlation) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

color_mapping <- type_counts %>%
  mutate(
    fill_color = case_when(
      correlation %in% c("UPUU", "UNUU") ~ "#596A98",
      correlation %in% c("UUPU", "UUNU") ~ "#6B886D",
      correlation %in% c("PPPU", "NNNU") ~ "#6e8072",
      correlation %in% c("UPPU", "UNNU") ~ "#6f8779",
      correlation %in% c("PPUU", "NNUU") ~ "#756488",
      correlation %in% c("PUPU", "NUNU") ~ "#837b64",
      correlation %in% c("PPPP", "NNNN") ~ "#e17e16",
      correlation %in% c("UPPP", "UNNN") ~ "#e18017",
      correlation %in% c("PUUU", "NUUU") ~ "#E41A1C",
      correlation %in% c("UPUP", "UNUN") ~ "#e48022",
      correlation %in% c("PUPP", "NUNN") ~ "#e77d14",
      correlation %in% c("UUPP", "UUNN") ~ "#e7851a",
      correlation %in% c("PUUP", "NUUN") ~ "#fb730f",
      correlation %in% c("UUUP", "UUUN") ~ "#FF7F00",
      correlation == "PPUP" ~ "#e4791b",
    )
  ) %>%
  pull(fill_color, name = correlation)

type_counts <- type_counts %>%
  mutate(
    border_color = if_else(correlation %in% c("PPPP", "NNNN"), "focus", "no_focus")
  )

p4<-ggplot(type_counts, 
           aes(x = reorder(correlation, count), 
               y = count, 
               fill = correlation,
               color = border_color)) +
    geom_bar(stat = "identity", width = 0.8,size = 1) +
    coord_flip() +
    geom_text(aes(label = count), hjust = -0.1,size = 6 / 2.835,color = "black") +
    #  geom_text(aes(label = star, y = count * 4),vjust = 0.78,color = "red",size = 10,fontface = "bold") +
    scale_fill_manual(values = color_mapping) +
    scale_color_manual(
        values = c("focus" = "red", "no_focus" = "transparent"),
        guide = "none"
    )+
    labs(
        title = "KO interaction with cross-kingdom network", 
        x = NULL, 
        y = NULL) +
theme_classic() + theme(
    text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 6, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5, color = "black"),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(color = "black", size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
)+
    scale_y_log10(
        expand = expansion(mult = c(0, 0.1)),
        breaks = 10^(0:5),
        labels = function(x) {
            ifelse(x >= 1e5, 
                   format(x, scientific = FALSE, big.mark = ","),
                   format(x, scientific = FALSE))
        }
    )
p4

##Figure4D
data <- read.csv("./data/Figure4-Group-enrichment.csv")
cor_values <- unique(data$Correlation)
combined_clean <- data %>%
  mutate(
    Correlation = factor(Correlation, levels = unique(cor_values)),
    Clade = factor(Clade, levels = c("Bacteria", "Fungi", "Protist"))
  )

heatmap_data <- combined_clean %>%
  tidyr::complete(
    nesting(Clade, Description),
    Correlation = unique(combined_clean$Correlation),
    fill = list(
      qvalue = 1,
      Count = 0,
      mean_qvalue = 1,
      max_count = 0
    )
  ) %>%
  mutate(
    log_q = -log10(ifelse(is.na(qvalue), 1, qvalue)),
    Significance = cut(
      qvalue,
      breaks = c(0, 0.001, 0.01, 0.05, 0.2, 1),
      labels = c("***", "**", "*", ".", "ns"),
      include.lowest = TRUE
    )
  ) %>%
  group_by(Clade, Correlation, Description) %>%
  summarise(
    log_q = max(log_q, na.rm = TRUE),
    max_count = max(Count, na.rm = TRUE), 
    Significance = last(na.omit(Significance)),
    .groups = "drop"
  ) 

nuun_data <- expand.grid(
  Clade = unique(heatmap_data$Clade),
  Correlation = "NUUN",
  Description = unique(heatmap_data$Description)
) %>%
  mutate(
    log_q = 0,
    max_count = 0,
    Significance = "ns",
  ) %>%
  distinct()

heatmap_data <- heatmap_data %>% bind_rows(nuun_data)
heatmap_data <- heatmap_data[!grepl("yeast", heatmap_data$Description), ]

facet_data <- heatmap_data %>%
  group_by(Clade) %>%
  group_modify(~{
    desc_nnnn <- .x %>% filter(Correlation == "NNNN" & Significance != "ns") %>% pull(Description)
    desc_pppp <- .x %>% filter(Correlation == "PPPP" & Significance != "ns") %>% pull(Description)
    desc_union <- union(desc_nnnn, desc_pppp)
    filter(.x, Description %in% desc_union)
  }) %>%
  ungroup()

corr_order <- levels(with(type_counts, reorder(correlation, count)))
facet_data$Correlation <- factor(facet_data$Correlation, levels = corr_order)

p5 <- ggplot(facet_data, aes(x = Description, y = Correlation)) +
    geom_point(
        data = ~ subset(.x, max_count > 0),
        aes(size = max_count),
        position = position_nudge(x = 0.08, y = -0.08),
        color = "gray20",
        alpha = 0.4,
        shape = 19,
        show.legend = FALSE
    ) +
    geom_point(
        data = ~ subset(.x, max_count > 0),
        aes(size = max_count, fill = log_q),
        color = "black",
        shape = 21,
        stroke = 0.3,
        alpha = 1 
    ) +
    scale_fill_gradientn(
        name = expression(-log[10](q[value])),
        colors = c("#FFFFE0", "#E41A1C"),
        values = scales::rescale(c(0, 5, 10, 20, 50)),
        limits = c(0, 50)
    ) +
    scale_size_continuous(
        name = expression(KO_count),
        range = c(0.3, 3),
        limits = c(0, 50),
        breaks = scales::extended_breaks(n = 5),
        guide = guide_legend(order = 1)
    ) +
    scale_y_discrete(drop = FALSE) +
    facet_grid(. ~ Clade, scales = "free_x", space = "free") +
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 6, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "none"
    )+
    labs(
        x = NULL,
        y = NULL
    ) 
p5



p_left <- cowplot::plot_grid(
    p1 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    NULL,
    p2 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    NULL,
    ncol = 1, nrow = 4, rel_heights = c(1, 0.1, 0.8, 0.1) 
)
p_top <- cowplot::plot_grid(
    p_left + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p3 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    ncol = 2, nrow = 1, rel_widths = c(1, 2.5) 
)
p_top

p_bottom <- cowplot::plot_grid(
    p4 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    NULL,
    p5 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 3, nrow = 1, rel_widths = c(1.15, 0.1, 2.2) 
)
p_bottom


p <- cowplot::plot_grid(
    p_top + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p_bottom + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    ncol = 1, nrow = 2, rel_widths = c(3, 2) 
)
p


width <- 17.5
height <- 20
ggsave(paste0("Figure4", ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0("Figure4", ".pdf"), p, width = width, height = height, units = "cm")
