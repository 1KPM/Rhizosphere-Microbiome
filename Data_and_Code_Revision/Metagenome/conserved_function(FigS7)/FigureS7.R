pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

library(ggplot2)
library(tidyverse)
library(cowplot)
library(patchwork)
library(RColorBrewer)

final_abundance <- read.csv("./data/All_level2_abundance.csv",row.names = 1)

top_superpathways <- final_abundance %>%
  rownames_to_column("superpathway") %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(total)) %>%
  pull(superpathway)

data_processed <- final_abundance %>%
  rownames_to_column("superpathway") %>%
  mutate(superpathway = factor(superpathway, levels = top_superpathways)) %>%
  group_by(superpathway) %>%
  summarise(across(where(is.numeric), sum)) %>%
  ungroup()

df_long <- data_processed %>%
  pivot_longer(-superpathway, names_to = "Order", values_to = "abundance")

fill_color <- read.csv("./data/superpathway_color.csv")
fill_colors <- fill_color$fill_color2
names(fill_colors) <- fill_color$superpathway

orders_list <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales","Others")
df_long$Order <- factor(df_long$Order, levels = orders_list)


p_all <- ggplot(df_long, aes(x = Order, y = abundance, fill = superpathway)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = fill_colors) +
  labs(title = "All", x = "", y = "Relative Abundance", fill = "Superpathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.key.size = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(size = 6)
  ) +
  guides(fill = guide_legend(ncol = 1))

all_superpathway_order <- levels(data_processed$superpathway)


kingdom <- c("Bacteria","Fungi","Protist")
plot_list <- list()
for (k in kingdom){
  final_abundance <- read.csv(paste0("./data/",k,"_level2_abundance.csv"),row.names = 1)
  data_processed <- final_abundance %>%
    rownames_to_column("superpathway") %>%
    mutate(superpathway = factor(superpathway, levels = all_superpathway_order)) %>%
    group_by(superpathway) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup()
  
  df_long <- data_processed %>%
    pivot_longer(-superpathway, names_to = "Order", values_to = "abundance")
  
  orders_list <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales","Others")
  df_long$Order <- factor(df_long$Order, levels = orders_list)
  
  p <- ggplot(df_long, aes(x = Order, y = abundance, fill = superpathway)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(title = k, x = "",  y = "",fill = "Superpathway") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5)
    ) + 
    guides(fill = guide_legend(ncol = 1))
  plot_list[[k]] <- p
}

kegg_plot <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")



family_levels <- c("Glycoside Hydrolases", "Glycosyl Transferases", "Polysaccharide Lyases",
                   "Carbohydrate Esterases", "Auxiliary Activities", "Carbohydrate-Binding Modules")
fill_colors <- colorRampPalette(brewer.pal(8, 'Set3'))(6)
names(fill_colors) <- family_levels
orders_list <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales","Others")

kingdom <- c("Bacteria", "Fungi", "Protist")
plot_list <- list()
for (k in kingdom){
  processed_data <- read.csv(paste0("./data/",k,"_cazyme_family_abundance.csv"),row.names = 1)
  data_processed <- processed_data %>%
    rownames_to_column("Family") %>%
    mutate(Family = factor(Family, levels = family_levels)) %>%
    complete(Family = family_levels, fill = as.list(setNames(rep(0, length(orders_list)), orders_list))) %>%
    distinct(Family, .keep_all = TRUE)
  
  df_long <- data_processed %>%
    pivot_longer(-Family, names_to = "Order", values_to = "abundance") %>%
    mutate(Order = factor(Order, levels = orders_list),
           Family = factor(Family, levels = family_levels))
  
  p <- ggplot(df_long, aes(x = Order, y = abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors, breaks = family_levels, drop = FALSE) +
    labs(title = k, x = "",  y = "",fill = "Cazyme Family") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(fill = guide_legend(ncol = 1))
  plot_list[[k]] <- p
}

cazy_plot <- wrap_plots(plot_list, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

combined_plot <- wrap_plots(kegg_plot,cazy_plot,ncol = 1)

width <-25
height <- 24
name <- paste0(dir_name,"/FigureS7")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")
ggsave(paste0(name, ".jpg"), combined_plot, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), combined_plot, width = width, height = height, units = "cm")