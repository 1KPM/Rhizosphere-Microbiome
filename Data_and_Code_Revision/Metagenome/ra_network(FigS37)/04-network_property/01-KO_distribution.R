# ******************************************************************************
# @File: 01-ko_distribution.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-02-28 14:08:49
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
dir_name <- "01-ko_distribution"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
target_clades <- c('Bacteria', 'Fungi', 'Protist')
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (typ in type) {
    first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
    second_prefix <- ifelse(typ == "all", "02", "01")
    
    dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
    
    hub_info <- read.csv(paste0(dir_path, '/', typ, '_network_hub_info.csv'),row.names = 1)
    
    hub_info <- hub_info %>%
      mutate(Clade = case_when(
        substr(.[[1]], 1, 1) == "b" ~ "Bacteria",
        substr(.[[1]], 1, 1) == "f" ~ "Fungi",
        substr(.[[1]], 1, 1) == "p" ~ "Protist",
        TRUE ~ NA_character_
      ))
    
    # All KO中各Kingdom的数量
    all_label <- paste0("All KOs (", nrow(hub_info), ")")
    all_ko_df <- hub_info %>%
        mutate(Clade = factor(Clade, levels = target_clades)) %>%
        filter(!is.na(Clade)) %>% 
        count(Clade, .drop = FALSE) %>%
        rename(Kingdom = Clade, Value = n) %>%
        mutate(Percent = Value / sum(Value), KO = all_label)
    
    
    # Hub KO中各Kingdom的数量
    hub_ko_list <- rownames(hub_info)[hub_info$roles != "Peripherals"]
    hub_ko_tax <- hub_info[hub_ko_list, ]
    hub_label <- paste0("Hub KOs (", nrow(hub_ko_tax), ")")
    
    hub_ko_df <- hub_ko_tax %>%
        mutate(Clade = factor(Clade, levels = target_clades)) %>%
        filter(!is.na(Clade)) %>% 
        count(Clade, .drop = FALSE) %>%
        rename(Kingdom = Clade, Value = n) %>%
        mutate(Percent = Value / sum(Value), KO = hub_label)
    
    ko_df <- bind_rows(all_ko_df, hub_ko_df)
    
    p <- ggplot(ko_df, aes(x = "", y = Percent, fill = Kingdom)) +
        geom_bar(stat = "identity") +
        labs(
            title = NULL,
            subtitle = NULL,
            x = NULL,
            y = NULL
        ) + 
        geom_text(
            aes(label = scales::percent(Percent, 0.1)),
            color = "black", position = position_stack(vjust = 0.5),
            size = 7 / 2.835
        ) +
        coord_polar(theta = "y") +
        facet_wrap(~ KO, strip.position = "bottom") +
        theme_void() +
        scale_fill_manual(values = color_manual) +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            panel.border = element_blank(), 
            panel.grid = element_blank(), 
            axis.ticks = element_blank(), 
            plot.title = element_text(size = 7, hjust = 0.5),
            strip.text = element_text(size = 7),
            panel.spacing = unit(-1, "lines"),
            legend.title = element_text(size = 7, color = 'black'), 
            legend.text = element_text(size = 6,  color = 'black'), 
            legend.key.size = unit(0.25, 'cm'),
            legend.position = "none",
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
        ) 
 
    width = 6; height = 3.8
    name <- paste0(dir_name, "/", typ, "_ko_distribution")
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}
# ------------------------------------------------------------------------------
