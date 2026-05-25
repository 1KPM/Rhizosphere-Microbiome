### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(reshape2)

### Define variable -----------------------------------------------------------
kingdom <- c('Bacteria', 'Fungi', 'Protist')
num <- 10
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
tax_all <- read.csv('data/All_core_ASV_taxonomy.csv', row.names = 1)
top_order <- read.csv('data/top_order.csv')
kegg_info <- read.csv('data/kegg_pathway_for_1KPM.csv')
# ------------------------------------------------------------------------------
### 1. Within
hub_info <- read.csv('data/Within_Microbiome_network_hub_info.csv', row.names = 1)
tmp_df <- merge(hub_info, tax_all, by = 'row.names')

for (k in kingdom) {
    raw_df <- tmp_df[tmp_df$Clade == k,]
    raw_df$Order <- ifelse(raw_df$Order %in% top_order[, k], raw_df$Order, 'Others')
    
    all_df <- data.frame(t(table(raw_df$Order)))
    hub_df <- data.frame(t(table(raw_df[raw_df$roles != 'Peripherals', 'Order'])))
    
    res_df <- merge(all_df[2:3], hub_df[2:3], by = 'Var2', all.x = T)
    res_df[is.na(res_df)] <- 0
    names(res_df) <- c('Taxonomy', paste0('All (', sum(all_df$Freq), ')'), paste0('Hub (', sum(hub_df$Freq), ')'))
    
    fin_df <- melt(res_df, id.vars = 'Taxonomy', variable.name = 'ASV', value.name = 'Value')
    color_df <- data.frame(
        Taxonomy = c(top_order[, k], 'Others'),
        Color = colorRampPalette(brewer.pal(9, 'Set1'))(num+1)
    )
    color_df <- color_df[color_df$Taxonomy %in% fin_df$Taxonomy,]
    fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)
    color_manual <- color_df$Color
    
    p <- ggplot(fin_df, aes(x = ASV, y = Value, fill = Taxonomy)) +
        geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
        labs(
            title = NULL,
            subtitle = NULL,
            x = NULL,
            y = "Proportion of ASVs (Within-kingdom)"
        ) + 
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = color_manual, name = k) +
        theme_bw() + theme(
            text = element_text(color = "black", size = 6),
            plot.title = element_text(size = 7, hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.title = element_text(size = 7),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
            panel.spacing = unit(0.1, "cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.key.size = unit(0.25, "cm"),
            legend.position = "right"
        )
    assign(paste0("p_within_", k), p)
}

p_within <- cowplot::plot_grid(
    p_within_Bacteria + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p_within_Fungi + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p_within_Protist + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 3, nrow = 1, rel_widths = c(1, 1, 1) 
)


### 2. Cross
hub_info <- read.csv('data/Cross_Microbiome_network_hub_info.csv', row.names = 1)
tmp_df <- merge(hub_info, tax_all, by = 'row.names')

for (k in kingdom) {
    raw_df <- tmp_df[tmp_df$Clade == k,]
    raw_df$Order <- ifelse(raw_df$Order %in% top_order[, k], raw_df$Order, 'Others')
    
    all_df <- data.frame(t(table(raw_df$Order)))
    hub_df <- data.frame(t(table(raw_df[raw_df$roles != 'Peripherals', 'Order'])))
    
    res_df <- merge(all_df[2:3], hub_df[2:3], by = 'Var2', all.x = T)
    res_df[is.na(res_df)] <- 0
    names(res_df) <- c('Taxonomy', paste0('All (', sum(all_df$Freq), ')'), paste0('Hub (', sum(hub_df$Freq), ')'))
    
    fin_df <- melt(res_df, id.vars = 'Taxonomy', variable.name = 'ASV', value.name = 'Value')
    color_df <- data.frame(
        Taxonomy = c(top_order[, k], 'Others'),
        Color = colorRampPalette(brewer.pal(9, 'Set1'))(num+1)
    )
    color_df <- color_df[color_df$Taxonomy %in% fin_df$Taxonomy,]
    fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)
    color_manual <- color_df$Color
    
    p <- ggplot(fin_df, aes(x = ASV, y = Value, fill = Taxonomy)) +
        geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
        labs(
            title = NULL,
            subtitle = NULL,
            x = NULL,
            y = "Proportion of ASVs (Cross-kingdom)"
        ) + 
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = color_manual, name = k) +
        theme_bw() + theme(
            text = element_text(color = "black", size = 6),
            plot.title = element_text(size = 7, hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.title = element_text(size = 7),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6, color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
            panel.spacing = unit(0.1, "cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.key.size = unit(0.25, "cm"),
            legend.position = "right"
        )
    assign(paste0("p_cross_", k), p)
}

p_cross <- cowplot::plot_grid(
    p_cross_Bacteria + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p_cross_Fungi + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p_cross_Protist + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 3, nrow = 1, rel_widths = c(1, 1, 1) 
)


# 3.meta
pathawy_color <- read.csv("data/top_superpathway.csv")
hub_info <- read.csv('data/Cross_meta_network_hub_info.csv', row.names = 1)
hub_info$KO <- gsub("^[bfp]", "", row.names(hub_info))
all_tmp_df <- merge(hub_info, kegg_info, by = 'KO')

top_pathawy <- data.frame(table(all_tmp_df$level2))
top_pathawy <- top_pathawy[order(-top_pathawy$Freq),]
pathawy_color <- pathawy_color[match(top_pathawy$Var1, pathawy_color$SuperPathway),]

all_data_df <- data.frame()
for (kin in kingdom) {
    raw_df <- all_tmp_df[all_tmp_df$Clade == kin,]
    all_df <- data.frame(t(table(raw_df$level2)))
    hub_df <- data.frame(t(table(raw_df[raw_df$roles != 'Peripherals', 'level2'])))
    
    res_df <- merge(all_df[2:3], hub_df[2:3], by = 'Var2', all.x = T)
    res_df[is.na(res_df)] <- 0
    names(res_df) <- c('Taxonomy', paste0('All (', sum(all_df$Freq), ')'), paste0('Hub (', sum(hub_df$Freq), ')'))
    
    fin_df <- melt(res_df, id.vars = 'Taxonomy', variable.name = 'ASV', value.name = 'Value')
    fin_df$Kingdom <- kin
    
    all_data_df <- bind_rows(all_data_df, fin_df)
}
all_data_df$Taxonomy <- factor(all_data_df$Taxonomy, levels = pathawy_color$SuperPathway)

p_meta <- ggplot(all_data_df, aes(x = ASV, y = Value, fill = Taxonomy)) +
    geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Proportion of ASVs (Cross-kingdom)"
    ) + 
    facet_grid(cols = vars(Kingdom), scales = "free", space = "free_x") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = pathawy_color$Color, name = "Super pathway") +
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )


p_amplicon <- cowplot::plot_grid(
    p_within + theme(plot.margin = margin(0, 0, 0.2, 0, "cm")),
    p_cross + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 1, nrow = 2
)

p <- cowplot::plot_grid(
    p_amplicon + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p_meta + theme(plot.margin = margin(0, 3, 0, 0, "cm")), 
    ncol = 1, nrow = 2, rel_heights = c(2, 1)
)

width <- 17
height <- 17
name <- paste0(dir_name, "/FigureS12")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
