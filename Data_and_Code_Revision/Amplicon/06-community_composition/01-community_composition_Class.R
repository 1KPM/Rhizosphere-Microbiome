### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "01-community_composition_Class"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(reshape2)
library(cowplot)

### Define variable -----------------------------------------------------------
unassigned_tax <- c("Unassigned", "uncultured", "unidentified", "")
color_manual <- c(colorRampPalette(brewer.pal(9, "Set1"))(11), "#666666")
legend_position <- "right"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
tree_file <- read.tree("../00-rawdata/phylogeny/tree_metadata_merge_info_final_align_tree.nwk")
barplot_16s <- read.csv("../01-sort_data/04-barplot_data/16S_barplot_data_relative_sorted.csv")
barplot_its <- read.csv("../01-sort_data/04-barplot_data/ITS_barplot_data_relative_sorted.csv")
barplot_pro <- read.csv("../01-sort_data/04-barplot_data/Protist_barplot_data_relative_sorted.csv")
tree_order_color <- read.csv("../01-sort_data/01-tree_color/Tree_top10_order_color.csv")
metadata_rs <- read.csv("../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../00-rawdata/metadata/tree_metadata_merge_info.csv")
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
sample_list <- Reduce(intersect, list(names(barplot_16s)[-1:-3], names(barplot_its)[-1:-3], names(barplot_pro)[-1:-3]))
tree_list <- intersect(tree_file$tip.label, metadata_rs[metadata_rs$FileID %in% sample_list, "TreeID"])

drop_list <- tree_file$tip.label[-match(tree_list, tree_file$tip.label)]
fin_tree <- drop.tip(tree_file, drop_list)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. Tree
sub_metadata_tree <- metadata_tree[metadata_tree$TreeID %in% tree_list, c("TreeID", "Order")]
sub_metadata_tree$Order <- ifelse(sub_metadata_tree$Order %in% tree_order_color$Order, sub_metadata_tree$Order, "Others")
sub_metadata_tree$Order <- factor(sub_metadata_tree$Order, levels = tree_order_color$Order)

group_info <- split(sub_metadata_tree$TreeID, sub_metadata_tree$Order)
fin_tree <- groupOTU(fin_tree, group_info)

color_vector <- with(tree_order_color, setNames(Color, Order))
p_tree <- ggtree(fin_tree, aes(color = group), branch.length = "none", size = 0.3, ladderize = F) +
    scale_color_manual(values = color_vector, name = "Plant") +
    labs(
        x = "Plant"
    ) + 
    theme(
        plot.title = element_text(size = 7, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 6,  color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        legend.key.size = unit(0.25, "cm"),
        legend.margin = margin(0, -0.5, 0, 0, unit = "cm"),
        legend.justification = "top",
        legend.position = legend_position
    )

# 2. 16S
raw_df <- barplot_16s

res_df <- raw_df[raw_df$Level == "Class",]
res_df <- res_df[order(-res_df$Mean), -2]

unassigned_df <- res_df[res_df$Taxonomy %in% unassigned_tax,]
unassigned_df$Taxonomy <- "Unassigned"

fin_df <- res_df[!(res_df$Taxonomy %in% unassigned_tax),]
top_Taxonomy <- fin_df$Taxonomy[1:10]
fin_df$Taxonomy <- ifelse(fin_df$Taxonomy %in% top_Taxonomy, fin_df$Taxonomy, "Others")
fin_df <- rbind(fin_df, unassigned_df)

top_df <- aggregate(fin_df[-1], by = list(Taxonomy = fin_df$Taxonomy), FUN = sum)
row.names(top_df) <- top_df$Taxonomy
top_df$Taxonomy <- factor(top_df$Taxonomy, levels = c(top_Taxonomy, "Others", "Unassigned"))
top_df <- top_df[order(top_df$Taxonomy),]

color_df <- data.frame(top_df[c("Taxonomy", "Mean")], Color = color_manual)
write.csv(color_df, paste0(dir_name, "/16S_top10Class_color_meanValue.csv"), row.names = F)

tmp_df <- data.frame(t(top_df[-1:-2]), check.names = F)
tmp_df <- merge(tmp_df, metadata_rs[c("FileID", "TreeID")], by.x = "row.names", by.y = "FileID", all.x = T)
tmp_df <- tmp_df[-1]

aggregate_df <- aggregate(tmp_df[-ncol(tmp_df)], by = list(Tree = tmp_df$TreeID), FUN = mean)
melt_df <- melt(aggregate_df, id.vars = "Tree", variable.name = "Taxonomy", value.name = "Value")
melt_df <- melt_df[melt_df$Tree %in% tree_list,]
melt_df$Tree <- factor(melt_df$Tree, levels = fin_tree$tip.label)
write.csv(melt_df, paste0(dir_name, "/16S_top10Class_barplot_data.csv"), row.names = F)

p_16s <- ggplot(melt_df, aes(x = Tree, y = Value, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = position_fill(reverse = F), width = 1) +
    labs(
        x = "",
        y = "Bacteria"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_manual, name = "Bacteria") +
    coord_flip() +
    theme(plot.title = element_text(size = 7, color = "black", hjust = 0.5),
          plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
          axis.title = element_text(size = 7, color = "black"),
          axis.line = element_blank(),
          axis.text = element_text(size = 6, color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 7, color = "black"),
          legend.text = element_text(size = 6,  color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          legend.position = legend_position) 

# 3. ITS
raw_df <- barplot_its

res_df <- raw_df[raw_df$Level == "Class",]
res_df <- res_df[order(-res_df$Mean), -2]

unassigned_df <- res_df[res_df$Taxonomy %in% unassigned_tax,]
unassigned_df$Taxonomy <- "Unassigned"

fin_df <- res_df[!(res_df$Taxonomy %in% unassigned_tax),]
top_Taxonomy <- fin_df$Taxonomy[1:10]
fin_df$Taxonomy <- ifelse(fin_df$Taxonomy %in% top_Taxonomy, fin_df$Taxonomy, "Others")
fin_df <- rbind(fin_df, unassigned_df)

top_df <- aggregate(fin_df[-1], by = list(Taxonomy = fin_df$Taxonomy), FUN = sum)
row.names(top_df) <- top_df$Taxonomy
top_df$Taxonomy <- factor(top_df$Taxonomy, levels = c(top_Taxonomy, "Others", "Unassigned"))
top_df <- top_df[order(top_df$Taxonomy),]

color_df <- data.frame(top_df[c("Taxonomy", "Mean")], Color = color_manual)
write.csv(color_df, paste0(dir_name, "/ITS_top10Class_color_meanValue.csv"), row.names = F)

tmp_df <- data.frame(t(top_df[-1:-2]), check.names = F)
tmp_df <- merge(tmp_df, metadata_rs[c("FileID", "TreeID")], by.x = "row.names", by.y = "FileID", all.x = T)
tmp_df <- tmp_df[-1]

aggregate_df <- aggregate(tmp_df[-ncol(tmp_df)], by = list(Tree = tmp_df$TreeID), FUN = mean)
melt_df <- melt(aggregate_df, id.vars = "Tree", variable.name = "Taxonomy", value.name = "Value")
melt_df <- melt_df[melt_df$Tree %in% tree_list,]
melt_df$Tree <- factor(melt_df$Tree, levels = fin_tree$tip.label)
write.csv(melt_df, paste0(dir_name, "/ITS_top10Class_barplot_data.csv"), row.names = F)

p_its <- ggplot(melt_df, aes(x = Tree, y = Value, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = position_fill(reverse = F), width = 1) +
    labs(
        x = "",
        y = "Fungi"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_manual, name = "Fungi") +
    coord_flip() +
    theme(plot.title = element_text(size = 7, color = "black", hjust = 0.5),
          plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
          axis.title = element_text(size = 7, color = "black"),
          axis.line = element_blank(),
          axis.text = element_text(size = 6, color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 7, color = "black"),
          legend.text = element_text(size = 6,  color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          legend.position = legend_position) 

# 4. Protist
raw_df <- barplot_pro

res_df <- raw_df[raw_df$Level == "Class",]
res_df <- res_df[order(-res_df$Mean), -2]

unassigned_df <- res_df[res_df$Taxonomy %in% unassigned_tax,]
unassigned_df$Taxonomy <- "Unassigned"

fin_df <- res_df[!(res_df$Taxonomy %in% unassigned_tax),]
top_Taxonomy <- fin_df$Taxonomy[1:10]
fin_df$Taxonomy <- ifelse(fin_df$Taxonomy %in% top_Taxonomy, fin_df$Taxonomy, "Others")
fin_df <- rbind(fin_df, unassigned_df)

top_df <- aggregate(fin_df[-1], by = list(Taxonomy = fin_df$Taxonomy), FUN = sum)
row.names(top_df) <- top_df$Taxonomy
top_df$Taxonomy <- factor(top_df$Taxonomy, levels = c(top_Taxonomy, "Others", "Unassigned"))
top_df <- top_df[order(top_df$Taxonomy),]

color_df <- data.frame(top_df[c("Taxonomy", "Mean")], Color = color_manual)
write.csv(color_df, paste0(dir_name, "/Protist_top10Class_color_meanValue.csv"), row.names = F)

tmp_df <- data.frame(t(top_df[-1:-2]), check.names = F)
tmp_df <- merge(tmp_df, metadata_rs[c("FileID", "TreeID")], by.x = "row.names", by.y = "FileID", all.x = T)
tmp_df <- tmp_df[-1]

aggregate_df <- aggregate(tmp_df[-ncol(tmp_df)], by = list(Tree = tmp_df$TreeID), FUN = mean)
melt_df <- melt(aggregate_df, id.vars = "Tree", variable.name = "Taxonomy", value.name = "Value")
melt_df <- melt_df[melt_df$Tree %in% tree_list,]
melt_df$Tree <- factor(melt_df$Tree, levels = fin_tree$tip.label)
write.csv(melt_df, paste0(dir_name, "/Protist_top10Class_barplot_data.csv"), row.names = F)

p_pro <- ggplot(melt_df, aes(x = Tree, y = Value, fill = Taxonomy)) +
    geom_bar(stat = "identity", position = position_fill(reverse = F), width = 1) +
    labs(
        x = "",
        y = "Protists"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = color_manual, name = "Protists") +
    coord_flip() +
    theme(plot.title = element_text(size = 7, color = "black", hjust = 0.5),
          plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
          axis.title = element_text(size = 7, color = "black"),
          axis.line = element_blank(),
          axis.text = element_text(size = 6, color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 7, color = "black"),
          legend.text = element_text(size = 6,  color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          legend.position = legend_position) 

legend_p_tree <- cowplot::get_legend(p_tree)
legend_p_16s <- cowplot::get_legend(p_16s)
legend_p_its <- cowplot::get_legend(p_its)
legend_p_pro <- cowplot::get_legend(p_pro)

new_p_tree <- p_tree + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
new_p_16s <- p_16s + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
new_p_its <- p_its + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))
new_p_pro <- p_pro + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "cm"))

legend_p <- cowplot::plot_grid(
    legend_p_tree,
    legend_p_16s,
    legend_p_its,
    legend_p_pro,
    ncol = 1
)


p <- plot_grid(
    legend_p,
    new_p_tree + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    new_p_16s + theme(plot.margin = unit(c(0, 0, 0, -0.5), "cm")),
    new_p_its + theme(plot.margin = unit(c(0, 0, 0, -0.5), "cm")), 
    new_p_pro + theme(plot.margin = unit(c(0, 0, 0, -0.5), "cm")), 
    ncol = 5, rel_widths = c(1.5, 0.8, 1, 1, 1), align = "h"
)

name <- paste0(dir_name, "/Barplot_Class")
width <- 10
height <- 16
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
