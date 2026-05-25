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


### Import data ----------------------------------------------------------------
maptree_edge <- read.csv('data/within_maptree_edge_absolute.csv')
maptree_node <- read.csv('data/within_maptree_node_absolute.csv')
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
    geom_point(aes(X1, X2, fill = Phylum, size = Degree ** 2), pch = 21, data = maptree_node) +
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
name <- paste0(dir_name, '/Figure3A')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p, width = width, height = height, dpi = 600, units = "cm", compression = "lzw")
# ------------------------------------------------------------------------------

