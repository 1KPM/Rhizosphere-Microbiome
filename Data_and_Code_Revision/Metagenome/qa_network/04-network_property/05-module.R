# ******************************************************************************
# @File: 07-module.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-31 10:56:35
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "05-module"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(patchwork)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
kingdom <- c('Bacteria', 'Fungi', 'Protist')
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. Module Correlation
all_data_df <- all_main_module <- all_main_edge <- NULL

for (typ in type) {
    first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
    second_prefix <- ifelse(typ == "all", "02", "01")
    network <- ifelse(typ == "all", "Whole", ifelse(typ == "inter", "Interkingdom", "Intrakingdom"))
    
    dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
    
    node <- read.csv(paste0(dir_path, '/', typ, '_node_mt2.csv'))
    edge <- read.csv(paste0(dir_path, '/', typ, '_edge_mt2.csv'))
    module <- read.csv(paste0(dir_path, '/', typ, '_module_mt2.csv'))
    
    module_freq <- as.data.frame(table(module$group)) %>%
        arrange(-Freq) %>%
        filter(Var1 != "mother_no" & Freq >= nrow(module) * 0.1)
        
    main_module <- module %>%
        select(-degree) %>%
        filter(group %in% module_freq$Var1) %>%
        mutate(
            Network = paste0(network, " network"),
            group = gsub("model_", "M", group)
            ) 
    
    main_edge <- main_module %>%
        rename(OTU_1 = ID, Module_OTU_1 = group) %>%
        inner_join(edge, by = "OTU_1") %>%
        inner_join(
            main_module %>%
                rename(OTU_2 = ID, Module_OTU_2 = group), by = "OTU_2"
            ) %>%
        select(X1, Y1, OTU_1, Module_OTU_1, X2, Y2, OTU_2, Module_OTU_2, weight, Correlation)
    
    data_df <- main_edge %>%
        mutate(
            Module_1 = pmin(Module_OTU_1, Module_OTU_2),
            Module_2 = pmax(Module_OTU_1, Module_OTU_2)
        ) %>%
        count(Module_1, Module_2, Correlation, name = "Count") %>%
        rename(
            Module_OTU_1 = Module_1, 
            Module_OTU_2 = Module_2
        ) %>%
        mutate(
            Network = paste0(network, " network"),
            Correlation = case_when(Correlation == "-" ~ "Negative", TRUE ~ "Positive")
        )
    
    all_data_df <- bind_rows(all_data_df, data_df)
    all_main_module <- bind_rows(all_main_module, main_module)
    all_main_edge <- bind_rows(all_main_edge, main_edge %>% mutate(Network = paste0(network, " network")))
}

data_df$Network <- factor(data_df$Network, levels = c("Whole network", "Interkingdom network", "Intrakingdom network"))

p_summary <- all_data_df %>%
    mutate(Interaction_Pair = paste(Module_OTU_1, Module_OTU_2, sep = "-")) %>%
    
    ggplot(aes(x = Interaction_Pair, y = Count, fill = Correlation)) +
    geom_text(
        mapping = aes(x = Interaction_Pair, y = Count + 200000, label = Count),
        position = position_dodge(0.9),
        size = 6 / 2.835
    ) +
    facet_grid(cols = vars(Network), scales = "free", space = "free_x") +
    geom_col(position = "dodge", color = "black", alpha = 0.8) +
    labs(
        title = NULL,
        x = NULL,
        y = "Count",
        fill = "Correlation"
    )+
    scale_fill_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F8766D")) +
    
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right"
    )

name <- paste0(dir_name, "/module_correlation")
width <- 17
height <- 6
ggsave(paste0(name, ".png"), p_summary, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p_summary, width = width, height = height, units = "cm")

write.csv(all_data_df, paste0(name, "_data.csv"), quote = F, row.names = F)
write.csv(all_main_module, paste0(name, "_main_module.csv"), quote = F, row.names = F)
write.csv(all_main_edge, paste0(name, "_main_edge.csv"), quote = F, row.names = F)


# 2. Module Composition
pathway2ko <- read.csv("../05-top10order_network/00-data/pathway2ko.csv")
pathway2superpathway <- read.csv("../05-top10order_network/00-data/pathway2superpathway.csv")
prefix_map <- tibble(
  prefix = c("b", "f", "p"),
  Clade = c("Bacteria", "Fungi", "Protist")
)
KO_info <- merge(pathway2ko, pathway2superpathway, by = "Pathway", all = TRUE) %>% 
  select(KO,Superpathway) %>%
  distinct() %>%
  crossing(prefix_map) %>%
  mutate(KO = str_c(prefix, KO)) %>%
  select(-prefix)
type <- c("all")
for (typ in type) {
  network <- ifelse(typ == "all", "Whole", ifelse(typ == "inter", "Interkingdom", "Intrakingdom"))
  main_module <- all_main_module %>%
    filter(Network == paste0(network, " network"))
  
  tmp_df <- merge(main_module, KO_info, by.x = 'ID', by.y = "KO")
  
  p_all_list <- NULL
  for (kin in kingdom) {
    raw_df <- tmp_df %>%
      filter(Clade == kin)
    
    all_df <- NULL
    for (mod in unique(raw_df$group)) {
      res_df <- raw_df %>%
        filter(group == mod)
      fin_df <- data.frame(t(table(res_df$Superpathway))) %>%
        select(-Var1) %>%
        rename(Superpathway = Var2, Count = Freq) %>%
        mutate(Module = paste0(mod, " (", sum(Count), ")"))
      
      all_df <- bind_rows(all_df, fin_df)
    }
    
    
    superpathway_color <- read.csv("../05-top10order_network/00-data/superpathway_color.csv")
    color_df <- data.frame(
      Superpathway = superpathway_color$superpathway,
      Color = superpathway_color$fill_color2
    )
    all_df <- all_df %>%
      mutate(Superpathway = factor(Superpathway, levels = color_df$Superpathway)) %>%
      tidyr::complete(Module, Superpathway, fill = list(Count = 0))
    
    p <- ggplot(all_df, aes(x = Module, y = Count, fill = Superpathway)) +
      geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
      labs(
        title = kin,
        subtitle = NULL,
        x = NULL,
        y = paste0("Proportion of KOs (", network, " network)")
      ) + 
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = color_df$Color, name = "Superpathway", drop = FALSE) +
      theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 8, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right"
      ) +
      guides(fill = guide_legend(ncol = 1))
    p_all_list[[kin]] <- p
    
    name <- paste0(dir_name, "/", typ, "_", kin, "_module_feature_distribution")
    write.csv(all_df, paste0(name, "_data.csv"), quote = F, row.names = F)
  }

  p1 <-  p_all_list[["Bacteria"]] + theme(legend.position = "none")
  p2 <-  p_all_list[["Fungi"]] + theme(legend.position = "none")
  p3 <-  p_all_list[["Protist"]] + theme(legend.position = "none")
  p_all <- p1 | p2 | p3 +
    plot_layout(guides = "collect") +
    theme(
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(10, 15, 10, 10)
    )
  name <- paste0(dir_name, "/", typ, "_module_feature_distribution")
  width <- 17
  height <- 7
  ggsave(paste0(name, ".png"), p_all, width = width, height = height, dpi = 600, units = "cm")
  ggsave(paste0(name, ".pdf"), p_all, width = width, height = height, units = "cm")
}

type <- c("inter")
for (typ in type) {
  network <- ifelse(typ == "all", "Whole", ifelse(typ == "inter", "Interkingdom", "Intrakingdom"))
  main_module <- all_main_module %>%
    filter(Network == paste0(network, " network"))
  
  tmp_df <- merge(main_module, KO_info, by.x = 'ID', by.y = "KO")
  
  p_inter_list <- NULL
  for (kin in kingdom) {
    raw_df <- tmp_df %>%
      filter(Clade == kin)
    
    all_df <- NULL
    for (mod in unique(raw_df$group)) {
      res_df <- raw_df %>%
        filter(group == mod)
      fin_df <- data.frame(t(table(res_df$Superpathway))) %>%
        select(-Var1) %>%
        rename(Superpathway = Var2, Count = Freq) %>%
        mutate(Module = paste0(mod, " (", sum(Count), ")"))
      
      all_df <- bind_rows(all_df, fin_df)
    }
    
    
    superpathway_color <- read.csv("../05-top10order_network/00-data/superpathway_color.csv")
    color_df <- data.frame(
      Superpathway = superpathway_color$superpathway,
      Color = superpathway_color$fill_color2
    )
    all_df <- all_df %>%
      mutate(Superpathway = factor(Superpathway, levels = color_df$Superpathway)) %>%
      tidyr::complete(Module, Superpathway, fill = list(Count = 0))
    
    p <- ggplot(all_df, aes(x = Module, y = Count, fill = Superpathway)) +
      geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
      labs(
        title = kin,
        subtitle = NULL,
        x = NULL,
        y = paste0("Proportion of KOs (", network, " network)")
      ) + 
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = color_df$Color, name = "Superpathway", drop = FALSE) +
      theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 8, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right"
      ) +
      guides(fill = guide_legend(ncol = 1))
    p_inter_list[[kin]] <- p
    
    name <- paste0(dir_name, "/", typ, "_", kin, "_module_feature_distribution")
    write.csv(all_df, paste0(name, "_data.csv"), quote = F, row.names = F)
  }
  
  rel_widths <- switch(
    typ,
    "all"   = c(1, 0.95, 1.15),
    "inter" = c(1, 0.95, 1.15),
    c(1, 0.8, 0.95)
  )
  
  p1 <-  p_inter_list[["Bacteria"]] + theme(legend.position = "none")
  p2 <-  p_inter_list[["Fungi"]] + theme(legend.position = "none")
  p3 <-  p_inter_list[["Protist"]] + theme(legend.position = "none")
  p_inter <- p1 | p2 | p3 +
    plot_layout(guides = "collect") +
    theme(
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(10, 15, 10, 10)
    )
  name <- paste0(dir_name, "/", typ, "_module_feature_distribution")
  width <- 17
  height <- 7
  ggsave(paste0(name, ".png"), p_inter, width = width, height = height, dpi = 600, units = "cm")
  ggsave(paste0(name, ".pdf"), p_inter, width = width, height = height, units = "cm")
}




type <- c("intra")
kingdom <- c('Bacteria', 'Fungi')
for (typ in type) {
  network <- ifelse(typ == "all", "Whole", ifelse(typ == "inter", "Interkingdom", "Intrakingdom"))
  main_module <- all_main_module %>%
    filter(Network == paste0(network, " network"))
  
  tmp_df <- merge(main_module, KO_info, by.x = 'ID', by.y = "KO")
  
  p_intra_list <- NULL
  for (kin in kingdom) {
    raw_df <- tmp_df %>%
      filter(Clade == kin)
    
    all_df <- NULL
    for (mod in unique(raw_df$group)) {
      res_df <- raw_df %>%
        filter(group == mod)
      fin_df <- data.frame(t(table(res_df$Superpathway))) %>%
        select(-Var1) %>%
        rename(Superpathway = Var2, Count = Freq) %>%
        mutate(Module = paste0(mod, " (", sum(Count), ")"))
      
      all_df <- bind_rows(all_df, fin_df)
    }
    
    
    superpathway_color <- read.csv("../05-top10order_network/00-data/superpathway_color.csv")
    color_df <- data.frame(
      Superpathway = superpathway_color$superpathway,
      Color = superpathway_color$fill_color2
    )
    all_df <- all_df %>%
      mutate(Superpathway = factor(Superpathway, levels = color_df$Superpathway)) %>%
      tidyr::complete(Module, Superpathway, fill = list(Count = 0))
    
    p <- ggplot(all_df, aes(x = Module, y = Count, fill = Superpathway)) +
      geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
      labs(
        title = kin,
        subtitle = NULL,
        x = NULL,
        y = paste0("Proportion of KOs (", network, " network)")
      ) + 
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = color_df$Color, name = "Superpathway", drop = FALSE) +
      theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 8, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8), 
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right"
      ) +
      guides(fill = guide_legend(ncol = 1))
    p_intra_list[[kin]] <- p
    
    name <- paste0(dir_name, "/", typ, "_", kin, "_module_feature_distribution")
    write.csv(all_df, paste0(name, "_data.csv"), quote = F, row.names = F)
  }

  p1 <-  p_intra_list[["Bacteria"]] + theme(legend.position = "none")
  p2 <-  p_intra_list[["Fungi"]] + theme(legend.position = "none")
  p_intra <- (p1 | p2 | plot_spacer()) + 
    plot_layout(guides = "collect", widths = c(10, 5, 10)) +
    theme(
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(10, 15, 10, 10)
    )
  name <- paste0(dir_name, "/", typ, "_module_feature_distribution")
  width <- 17
  height <- 7
  ggsave(paste0(name, ".png"), p_intra, width = width, height = height, dpi = 600, units = "cm")
  ggsave(paste0(name, ".pdf"), p_intra, width = width, height = height, units = "cm")
}
# ------------------------------------------------------------------------------

p_bac_all  <- p_all_list[["Bacteria"]]
p_fun_all  <- p_all_list[["Fungi"]]
p_pro_all  <- p_all_list[["Protist"]]
p_bac_inter <- p_inter_list[["Bacteria"]]
p_fun_inter <- p_inter_list[["Fungi"]]
p_pro_inter <- p_inter_list[["Protist"]]
p_bac_intra <- p_intra_list[["Bacteria"]]
p_fun_intra <- p_intra_list[["Fungi"]]


layout <- "
AABBCC
DDEEFF
GGH###
"
row1 <- wrap_plots(
  p_bac_all, p_fun_all, p_pro_all,
  p_bac_inter, p_fun_inter, p_pro_inter,
  p_bac_intra, p_fun_intra, plot_spacer()
) +
  plot_layout(
    design = layout,
    guides = "collect"
  ) +
  guides(fill = guide_legend(ncol = 1))

final_plot <- wrap_plots(row1, p_summary, ncol = 1, heights = c(4, 1))

name <- paste0(dir_name, "/All_module_distribution")
width <- 17
height <- 20
ggsave(paste0(name, ".png"), final_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), final_plot, width = width, height = height, units = "cm")




all <- all_main_edge %>%
    filter(Network == "Whole network") %>%
    group_by(Module_OTU_1) %>%
    summarise(across(c(X1, Y1), mean, na.rm = TRUE))
all
inter <- all_main_edge %>%
    filter(Network == "Interkingdom network") %>%
    group_by(Module_OTU_1) %>%
    summarise(across(c(X1, Y1), mean, na.rm = TRUE))
inter
intra <- all_main_edge %>%
    filter(Network == "Intrakingdom network") %>%
    group_by(Module_OTU_1) %>%
    summarise(across(c(X1, Y1), mean, na.rm = TRUE))
intra
