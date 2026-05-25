### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Import package
library(tidyverse)
library(RColorBrewer)


### Import data ----------------------------------------------------------------
maptree_edge <- read.csv('all_edge_mt2.csv')
maptree_node <- read.csv('all_node_mt2.csv')
maptree_node$Clade <- ifelse(substr(maptree_node[, 1], 1, 1) == "b", "Bacteria", ifelse(substr(maptree_node[, 1], 1, 1) == "f", "Fungi", "Protist"))

# ------------------------------------------------------------------------------


### get results ----------------------------------------------------------------

color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)

p <- ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Correlation), 
               alpha = 0.01, linewidth = 0.01, data = maptree_edge) + 
  geom_point(aes(X1, X2, fill = Clade, size = Degree ** 0.6), pch = 21, data = maptree_node) +
  scale_color_manual(values = c("+" = "#00BFC4", "-" = "#F8766D"))+
  scale_fill_manual(values = color_manual) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  theme_classic() + 
  theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
        plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
        axis.title = element_text(size = 7, color = 'black'), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(size = 6, color = 'black'), 
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        legend.position = 'none')

width <- 50
height <- 50
name <- paste0('all_network(degree^0.6)-50-50')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p, width = width, height = height, dpi = 600, units = "cm", compression = "lzw")
# ------------------------------------------------------------------------------

