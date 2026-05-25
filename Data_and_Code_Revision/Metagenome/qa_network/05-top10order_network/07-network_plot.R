# ******************************************************************************
# @File: 08-network_plot.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-12 16:32:47
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "07-network_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
net <- "inter"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
# ------------------------------------------------------------------------------

order_dimensions <- data.frame(
  Order = c("Fabales", "Rosales", "Malpighiales", "Gentianales", "Sapindales", 
            "Lamiales", "Malvales", "Myrtales", "Arecales", "Asparagales"),
  Width = c(30, 30, 38, 35, 37, 40, 40, 35, 38, 40),
  Height = c(40, 45, 35, 35, 25, 32, 40, 35, 38, 20)
)

for (i in 1:nrow(order_dimensions)) {
  ord <- order_dimensions$Order[i]
  width <- order_dimensions$Width[i]
  height <- order_dimensions$Height[i]
  
  tmp_path <- paste0("03-get_top10order_network_property/", ord, "_", net)
  
  maptree_edge <- read.csv(paste0(tmp_path, "_edge_mt2.csv"), header = T)
  maptree_node <- read.csv(paste0(tmp_path, "_node_mt2.csv"), header = T)
  maptree_node$Clade <- ifelse(substr(maptree_node[, 1], 1, 1) == "b", "Bacteria", 
                               ifelse(substr(maptree_node[, 1], 1, 1) == "f", "Fungi", "Protist"))
  
  color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
  
  p <- ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Correlation), 
                 alpha = 0.01, linewidth = 0.01, data = maptree_edge) + 
    geom_point(aes(X1, X2, fill = Clade, size = Degree ** 0.6), pch = 21, data = maptree_node) +
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
  
  name <- paste0(dir_name, '/', ord, '_', net, '_network_maptree_kingdom')
  
  ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
  ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}
# ------------------------------------------------------------------------------

