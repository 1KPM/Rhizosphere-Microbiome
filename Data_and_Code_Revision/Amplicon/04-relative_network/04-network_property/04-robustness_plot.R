# ******************************************************************************
# @File: 04-robustness_plot.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-09 18:49:53
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
dir_name <- "04-robustness_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(reshape2)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. Core ASV 累积去除后的Robustness
all_data_df <- NULL
for (typ in type) {
    file_path <- paste0("03-robustness/", typ, "_robustness_by_cumulative_asv.csv")
    
    tmp_df <- read.csv(file_path, header = T)
    names(tmp_df) <- c("Counts", "Target", "Random")
    
    data_df <- reshape2::melt(tmp_df, id.vars = 'Counts', variable.name = 'Type', value.name = 'Robustness')
    data_df$Network <- ifelse(typ == "all", "Whole network", ifelse(typ == "inter", "Interkingdom network", "Intrakingdom network"))
    
    all_data_df <- rbind(all_data_df, data_df)
}

name <- paste0(dir_name, "/robustness_by_cumulative_asv")
write.csv(all_data_df, paste0(name, ".csv"), quote = F, row.names = F)

all_data_df$Network <- factor(all_data_df$Network, levels = c("Whole network", "Interkingdom network", "Intrakingdom network"))
p <- ggplot(all_data_df, aes(x = Counts, y = Robustness, color = Type)) + 
    geom_line() +
    geom_point(size = 0.5) +
    labs(
        x = 'Removed ASV number',
        y = 'Roubustness'
    ) +
    theme_bw() +
    facet_wrap(~ Network, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c('#c01f25','#57657c'), name = "Remove type") +
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          
          legend.position = 'right')



width <- 17
height <- 5
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

# 2. Core ASV 按Kingdom去除后的Robustness
all_data_df <- NULL
for (typ in type) {
    file_path <- paste0("03-robustness/", typ, "_robustness_by_kingdom.csv")
    
    tmp_df <- read.csv(file_path, header = T)
    names(tmp_df) <- c("Kingdom", "Counts", "Target", "Random")
    tmp_df <- tmp_df[tmp_df$Counts != 0,]
    tmp_df <- tmp_df %>%
        mutate(Kingdom = case_when(Kingdom == "Protist" ~ "Protists", TRUE ~ Kingdom))
    
    
    tmp_df$Group <- paste0(tmp_df$Kingdom, " (", tmp_df$Counts, ")")
    
    data_df <- reshape2::melt(tmp_df[c("Group", "Random", "Target")], id.vars = "Group", variable.name = 'Type', value.name = 'Robustness')
    data_df$Network <- ifelse(typ == "all", "Whole network", ifelse(typ == "inter", "Interkingdom network", "Intrakingdom network"))
    
    all_data_df <- rbind(all_data_df, data_df)
}

name <- paste0(dir_name, "/robustness_by_kingdom")
write.csv(all_data_df, paste0(name, ".csv"), quote = F, row.names = F)

all_data_df$Network <- factor(all_data_df$Network, levels = c("Whole network", "Interkingdom network", "Intrakingdom network"))

p <- ggplot(
    data = all_data_df, 
    mapping = aes(x = Group, y = Robustness, fill = Type)
) +
    labs(
        x = "",
        y = "Roubustness"
    ) + 
    geom_bar(
        stat = "identity",
        position = "dodge",
        width = 0.68
    ) +
    geom_text(
        data = all_data_df,
        mapping = aes(x = Group, y = Robustness * 1.1, label = round(Robustness, 3)),
        position = position_dodge(0.9),
        size = 6 / 2.835
    ) +
    scale_fill_manual(values = c('#57657c', '#c01f25'), name = "Remove type") +
    facet_grid(cols = vars(Network), scales = "free", space = "free_x") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 7, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 6,  color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )
width <- 17
height <- 6.5
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------