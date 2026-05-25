### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "03-community_composition_Order_in_plant_top10order"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)

### Define variable -----------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
metadata_tree <- read.csv("../00-rawdata/metadata/tree_metadata_merge_info.csv")
tree_order_color <- read.csv("../01-sort_data/01-tree_color/Tree_top10_order_color.csv")
color_manual <- c(colorRampPalette(brewer.pal(9, "Set1"))(11), "#666666")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. 16S
amp <- "16S"

color_df <- read.csv(paste0("02-community_composition_Order/", amp, "_top10Order_color_meanValue.csv"))
raw_df <- read.csv(paste0("02-community_composition_Order/", amp, "_top10Order_barplot_data.csv"))

res_df <- merge(raw_df, metadata_tree, by.x = "Tree", by.y = "TreeID")
res_df$Order <- ifelse(res_df$Order %in% tree_order_color$Order, res_df$Order, "Others")

fin_df <- aggregate(res_df["Value"], by = list(Order = res_df$Order, Taxonomy = res_df$Taxonomy), FUN = mean)
fin_df$Order <- factor(fin_df$Order, levels = tree_order_color$Order)
fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)

p_16s <- ggplot(fin_df, aes(x = Order, y = Value, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = position_fill(reverse = F), width = 0.68) +
    labs(
        x = "",
        y = "Relative abundance (16S)"
    ) +
    
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_df$Color, name = "Bacteria") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 7, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 6,  color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "right"
    )

# 2. ITS
amp <- "ITS"

color_df <- read.csv(paste0("02-community_composition_Order/", amp, "_top10Order_color_meanValue.csv"))
raw_df <- read.csv(paste0("02-community_composition_Order/", amp, "_top10Order_barplot_data.csv"))

res_df <- merge(raw_df, metadata_tree, by.x = "Tree", by.y = "TreeID")
res_df$Order <- ifelse(res_df$Order %in% tree_order_color$Order, res_df$Order, "Others")

fin_df <- aggregate(res_df["Value"], by = list(Order = res_df$Order, Taxonomy = res_df$Taxonomy), FUN = mean)
fin_df$Order <- factor(fin_df$Order, levels = tree_order_color$Order)
fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)

p_its <- ggplot(fin_df, aes(x = Order, y = Value, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = position_fill(reverse = F), width = 0.68) +
    labs(
        x = "",
        y = "Relative abundance (ITS)"
    ) +
    
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_df$Color, name = "Fungi") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 7, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 6,  color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "right"
    )

# 3. Protist
amp <- "Protist"

color_df <- read.csv(paste0("02-community_composition_Order/", amp, "_top10Order_color_meanValue.csv"))
raw_df <- read.csv(paste0("02-community_composition_Order/", amp, "_top10Order_barplot_data.csv"))

res_df <- merge(raw_df, metadata_tree, by.x = "Tree", by.y = "TreeID")
res_df$Order <- ifelse(res_df$Order %in% tree_order_color$Order, res_df$Order, "Others")

fin_df <- aggregate(res_df["Value"], by = list(Order = res_df$Order, Taxonomy = res_df$Taxonomy), FUN = mean)
fin_df$Order <- factor(fin_df$Order, levels = tree_order_color$Order)
fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)

p_pro <- ggplot(fin_df, aes(x = Order, y = Value, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = position_fill(reverse = F), width = 0.68) +
    labs(
        x = "",
        y = "Relative abundance (18S)"
    ) +
    
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_df$Color, name = "Protists") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 7, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 6,  color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "right"
    )

# 4. Merge figure
# Extract figure legend
legend_16s <- get_legend(p_16s)
legend_its <- get_legend(p_its)
legend_pro <- get_legend(p_pro)

# remove figure legend
p_16s <- p_16s + theme(legend.position = "none")
p_its <- p_its + theme(legend.position = "none")
p_pro <- p_pro + theme(legend.position = "none")

# patchwork merge main figure
p <- (p_16s + p_its + p_pro) + 
    plot_layout(nrow = 3) &
    theme(plot.margin = margin(0, 0, 0, 0))

# merge legend
legend <- plot_grid(
    legend_16s,
    legend_its,
    legend_pro,
    nrow = 3,
    align = "v"
)

# merge main figure and legend
p_all <- plot_grid(
    p,
    legend + theme(plot.margin = margin(-0.5, 0, 0, -0.5, unit = "cm")),
    ncol = 2,
    rel_widths = c(3:2)
)

name <- paste0(dir_name, "/Barplot_Order_in_plant_top10_order")
width <- 9
height <- 16
ggsave(paste0(name, ".pdf"), p_all, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------
