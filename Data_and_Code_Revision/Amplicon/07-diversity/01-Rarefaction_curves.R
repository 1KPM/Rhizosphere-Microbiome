### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(1994)

# Import Packages
library(tidyverse)
library(reshape2)
library(RColorBrewer)


# Create Directory
dir_name <- "01-Rarefaction_curves"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Import Data ----------------------------------------------------------------
metadata_rs <- read.csv("../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv")
tree_order_color <- read.csv("../01-sort_data/01-tree_color/Tree_top10_order_color.csv")
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
data_df <- data.frame()
for (amp in amplicon) {
    tmp_df <- read.csv(paste0("../00-rawdata/diversity/",amp, "/shannon.csv"), row.names = 1)
    
    raw_df <- data.frame(t(tmp_df[-ncol(tmp_df)]))
    raw_df$Depth <- gsub('depth.', '', row.names(raw_df))
    raw_df$Depth <- paste0('Depth_', gsub("_iter\\.\\d+", '', raw_df$Depth))
    
    res_df <- aggregate(raw_df[-ncol(raw_df)], by = list(Depth = raw_df$Depth), FUN = mean)
    row.names(res_df) <- res_df$Depth
    res_df <- data.frame(t(res_df[-1]))
    
    group_df <- metadata_rs
    rownames(group_df) <- group_df$FileID

    merge_df <- merge(res_df, group_df["Order"], by = 'row.names')
    merge_df <- merge_df[-1]
    merge_df$Order <- ifelse(merge_df$Order %in% tree_order_color$Order, merge_df$Order, 'Others')
    
    aggregate_df <- aggregate(merge_df[-ncol(merge_df)], by = list(Order = merge_df$Order), FUN = function(x) mean(x, na.rm = TRUE))
    
    fin_df <- melt(aggregate_df, id.var = 'Order', variable.name = 'Depth', value.name = 'Shannon')
    fin_df$Depth <- as.numeric(gsub('Depth_', '', fin_df$Depth))
    fin_df$Order <- factor(fin_df$Order, levels = tree_order_color$Order)
    
    fin_df$Group <- ifelse(amp == "Protist", "18S", amp)
    data_df <- bind_rows(data_df, fin_df)
}
name <- paste0(dir_name, "/shannon")
write.csv(data_df, paste0(name, ".csv"), quote = F, row.names = F)
data_df$Group <- factor(data_df$Group, levels = c("16S", "ITS", "18S"))

p_shannon <- ggplot(data_df, aes(x = Depth, y = Shannon, color = Order), size = 3) + 
    geom_line() + 
    geom_vline(aes(xintercept = 10000), linetype = "dashed", linewidth = 0.5, color = "#666")+
    labs(
        x = 'Sampling depth',
        y = paste0('Shannon index')
    ) + 
    facet_grid(rows = vars(Group), scales = "free", space = "free_x") +
    scale_color_manual(values = tree_order_color$Color) +
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.3, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "none"
    )


data_df <- data.frame()
for (amp in amplicon) {
    tmp_df <- read.csv(paste0("../00-rawdata/diversity/",amp, "/observed_features.csv"), row.names = 1)
    
    raw_df <- data.frame(t(tmp_df[-ncol(tmp_df)]))
    raw_df$Depth <- gsub('depth.', '', row.names(raw_df))
    raw_df$Depth <- paste0('Depth_', gsub("_iter\\.\\d+", '', raw_df$Depth))
    
    res_df <- aggregate(raw_df[-ncol(raw_df)], by = list(Depth = raw_df$Depth), FUN = mean)
    row.names(res_df) <- res_df$Depth
    res_df <- data.frame(t(res_df[-1]))
    
    group_df <- metadata_rs
    rownames(group_df) <- group_df$FileID
    
    merge_df <- merge(res_df, group_df["Order"], by = 'row.names')
    merge_df <- merge_df[-1]
    merge_df$Order <- ifelse(merge_df$Order %in% tree_order_color$Order, merge_df$Order, 'Others')
    
    aggregate_df <- aggregate(merge_df[-ncol(merge_df)], by = list(Order = merge_df$Order), FUN = function(x) mean(x, na.rm = TRUE))
    
    fin_df <- melt(aggregate_df, id.var = 'Order', variable.name = 'Depth', value.name = 'Value')
    fin_df$Depth <- as.numeric(gsub('Depth_', '', fin_df$Depth))
    fin_df$Order <- factor(fin_df$Order, levels = tree_order_color$Order)
    
    fin_df$Group <- ifelse(amp == "Protist", "18S", amp)
    data_df <- bind_rows(data_df, fin_df)
}
name <- paste0(dir_name, "/observed_features")
write.csv(data_df, paste0(name, ".csv"), quote = F, row.names = F)
data_df$Group <- factor(data_df$Group, levels = c("16S", "ITS", "18S"))

p_obs <- ggplot(data_df, aes(x = Depth, y = Value, color = Order), size = 3) + 
    geom_line() + 
    labs(
        x = 'Sampling depth',
        y = paste0('Observed features')
    ) + 
    facet_grid(rows = vars(Group), scales = "free", space = "free_x") +
    scale_color_manual(values = tree_order_color$Color) +
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.3, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "none"
    )

p_legend <- ggplot(data_df, aes(x = Depth, y = Value, color = Order), size = 3) + 
    geom_line() + 
    labs(
        x = 'Sampling depth',
        y = paste0('Observed features')
    ) + 
    facet_grid(rows = vars(Group), scales = "free", space = "free_x") +
    scale_color_manual(values = tree_order_color$Color) +
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "right"
    )


legend_p <- cowplot::get_legend(p_legend)


p <- cowplot::plot_grid(
    p_shannon,
    p_obs + theme(plot.margin = margin(0, 0, 0, 0.5, "cm")),
    legend_p,
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 3, nrow = 1, rel_widths = c(3, 3, 1) 
)

width <- 17
height <- 15

name <- paste0(dir_name, "/Rarefaction_curves")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

# ------------------------------------------------------------------------------



