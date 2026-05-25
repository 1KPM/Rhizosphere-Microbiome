### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(agricolae)

ord <- "Rosales"
### Import data ----------------------------------------------------------------
maptree_edge <- read.csv(paste0('data/cross_maptree_edge_', ord, '_absolute.csv'))
maptree_node <- read.csv(paste0('data/cross_maptree_node_', ord, '_absolute.csv'))
top_5percent_taxa <- read.csv('data/top_5percent_taxa_absolute.csv')

# ------------------------------------------------------------------------------


### get results ----------------------------------------------------------------
top_16s_phylum <- top_5percent_taxa$top_16s_phylum[top_5percent_taxa$top_16s_phylum != '']
top_its_phylum <- top_5percent_taxa$top_its_phylum[top_5percent_taxa$top_its_phylum != '']
top_pro_phylum <- top_5percent_taxa$top_pro_phylum[top_5percent_taxa$top_pro_phylum != '']

maptree_node[maptree_node$Clade == 'Bacteria', 'Phylum'] <- 
    ifelse(maptree_node[maptree_node$Clade == 'Bacteria', 'Phylum'] %in% top_16s_phylum, 
           maptree_node[maptree_node$Clade == 'Bacteria', 'Phylum'], 'Other bacteria')

maptree_node[maptree_node$Clade == 'Fungi', 'Phylum'] <- 
    ifelse(maptree_node[maptree_node$Clade == 'Fungi', 'Phylum'] %in% top_its_phylum, 
           maptree_node[maptree_node$Clade == 'Fungi', 'Phylum'], 'Other Fungi')

maptree_node[maptree_node$Clade == 'Protist', 'Phylum'] <- 
    ifelse(maptree_node[maptree_node$Clade == 'Protist', 'Phylum'] %in% top_pro_phylum, 
           maptree_node[maptree_node$Clade == 'Protist', 'Phylum'], 'Other Protist')

maptree_node$Phylum <- 
    factor(maptree_node$Phylum, levels = c(top_16s_phylum, 'Other bacteria', top_its_phylum, 'Other Fungi',
                                           top_pro_phylum, 'Other Protist'))


color_phylum <- colorRampPalette(brewer.pal(12, 'Paired'))(13)

p <- ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Correlation), 
                 alpha = 0.01, linewidth = 0.1, data = maptree_edge) + 
    geom_point(aes(X1, X2, fill = Phylum, size = Degree), pch = 21, data = maptree_node) +
    scale_fill_manual(values = color_phylum) +
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(breaks = NULL) +
    theme_classic() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'right')

width <- 24
height <- 18
name <- paste0(dir_name, '/FigureS15D')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p, width = width, height = height, dpi = 600, units = "cm", compression = "lzw")
# ------------------------------------------------------------------------------


all_info <- read.csv(paste0("data/cross_network_hub_info_", ord, "_absolute.csv"), row.names = 1)
all_info$Clade <- ifelse(grepl("bASV", row.names(all_info)), "Bacteria", 
                         ifelse(grepl("fASV", row.names(all_info)), "Fungi", "Protist"))

hub_info <- all_info[all_info$roles != "Peripherals",]

all_df <- data.frame(table(all_info$Clade))
names(all_df) <- c("Kingdom", "Value")
all_df$Ratio <- all_df$Value / sum(all_df$Value)

hub_df <- data.frame(table(hub_info$Clade))
names(hub_df) <- c("Kingdom", "Value")
hub_df$Ratio <- hub_df$Value / sum(hub_df$Value)


color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
all_df$Type <- paste0("All ASVs (", sum(all_df$Value), ")")
hub_df$Type <- paste0("Hub ASVs (", sum(hub_df$Value), ")")
data_df <- bind_rows(all_df, hub_df)

p <- ggplot(data_df, aes(x = '', y = Value, fill = Kingdom)) +
    geom_bar(stat = "identity") +
    geom_text(
        aes(label = scales::percent(Ratio, 0.1)),
        color = "black", position = position_stack(vjust = 0.5),
        size = 7 / 2.835
    ) +
    coord_polar(theta = "y") +
    facet_wrap(~ Type, strip.position = "bottom", scales = "free",) +
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
        panel.spacing = unit(0.5, "lines"),
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        legend.key.size = unit(0.25, 'cm'),
        legend.position = "none"
    ) 


width <- 7
height <- 2.5
name <- paste0(dir_name, "/FigureS15E")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

# ------------------------------------------------------------------------------
hub_info <- read.csv(paste0("data/cross_network_hub_info_", ord, "_absolute.csv"), row.names = 1)
tax_all <- read.csv('data/All_core_ASV_taxonomy.csv', row.names = 1)
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)


raw_df <- merge(hub_info, tax_all, by = 'row.names')
names(raw_df)[2:4] <- c('Degree', 'Closeness', 'Betweenness')
raw_df$Clade <- factor(raw_df$Clade, levels = c('Bacteria', 'Fungi', 'Protist'))


fin_df <- raw_df[c('Clade', "Degree")]
names(fin_df) <- c('Group', 'Value')

max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
mean_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = mean)

n <- nrow(mean_df)
dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
dunn_res_df <- data.frame(dunn_res[-1])

pvalue_df <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
    for(j in (i + 1):n){ 
        k <- k + 1
        pvalue_df[i,j] <- dunn_res_df$p.adj[k]
        pvalue_df[j,i] <- dunn_res_df$p.adj[k]
    }
}

letter_df <- orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = TRUE)
letter_df <- letter_df[levels(fin_df$Group),]
letter_vector <- letter_df$groups

p <- ggplot(fin_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = max_df, aes(y = Value * c(1.1, 1.3, 1.5), label = letter_vector), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(
        x = '',
        y = 'Degree centrality of within-kingdom networks',
    ) + 
    theme_bw() + 
    coord_flip() +
    scale_color_manual(values = color_manual) + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

width <- 7
height <- 3
name <- paste0(dir_name, '/FigureS15F')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p, width = width, height = height, dpi = 600, units = "cm", compression = "lzw")
write.csv(dunn_res_df, paste0(name, '.csv'), row.names = F)