pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(ggplot2)
library(tidyverse)
library(cowplot)
library(patchwork)
library(RColorBrewer)

final_abundance <- read.csv("./data/FigureS5-All_level2_abundance.csv",row.names = 1)

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

# 
# p_all <- ggplot(df_long, aes(x = Order, y = abundance, fill = superpathway)) +
#   geom_bar(stat = "identity", position = "fill") +
#   scale_fill_manual(values = fill_colors) +
#   labs(title = "All", x = "", y = "Relative Abundance", fill = "Superpathway") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, color = "black"),
#     axis.text.y = element_text(color = "black"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
#     panel.background = element_blank(),
#     axis.line = element_line(color = "black"),
#     legend.position = "right",
#     legend.key.size = unit(0.3, "cm"),
#     legend.key.height = unit(0.3, "cm"),
#     legend.spacing.y = unit(0.1, "cm"),
#     legend.text = element_text(size = 7)
#   ) +
#   guides(fill = guide_legend(ncol = 1))

all_superpathway_order <- levels(data_processed$superpathway)


kingdom <- c("Bacteria","Fungi","Protist")

### 1. Bacteria 
k <- "Bacteria"
final_abundance <- read.csv(paste0("./data/FigureS5-",k,"_level2_abundance.csv"),row.names = 1)
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

p1 <- ggplot(df_long, aes(x = Order, y = abundance, fill = superpathway)) +
geom_bar(stat = "identity", position = "fill") +
scale_fill_manual(values = fill_colors) +
labs(title = k, x = "",  y = "Relative abundance",fill = "Superpathway") +
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
        legend.position = "none"
    )

### 2. Fungi 
k <- "Fungi"
final_abundance <- read.csv(paste0("./data/FigureS5-",k,"_level2_abundance.csv"),row.names = 1)
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

p2 <- ggplot(df_long, aes(x = Order, y = abundance, fill = superpathway)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(title = k, x = "",  y = "",fill = "Superpathway") +
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
        legend.position = "none"
    )

### 3. Protist 
k <- "Protist"
final_abundance <- read.csv(paste0("./data/FigureS5-",k,"_level2_abundance.csv"),row.names = 1)
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

p3 <- ggplot(df_long, aes(x = Order, y = abundance, fill = superpathway)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(title = k, x = "",  y = "",fill = "Superpathway") +
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
    ) +
    guides(fill = guide_legend(ncol = 1))

p3
legend_p <- cowplot::get_legend(p3)
p123 <- cowplot::plot_grid(
    p1 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p2 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p3 + theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "none"), 
    labels = c("a", "b", "c"), label_size = 10, label_fontface = "bold",
    legend_p,
    ncol = 4, nrow = 1
)
p123


family_levels <- c("Glycoside Hydrolases", "Glycosyl Transferases", "Polysaccharide Lyases",
                   "Carbohydrate Esterases", "Auxiliary Activities", "Carbohydrate-Binding Modules")
fill_colors <- colorRampPalette(brewer.pal(8, 'Set3'))(6)
names(fill_colors) <- family_levels
orders_list <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales","Others")

kingdom <- c("Bacteria", "Fungi", "Protist")

# 4. Bacteria
k <- "Bacteria"
processed_data <- read.csv(paste0("./data/FigureS5-",k,"_cazyme_family_abundance.csv"),row.names = 1)
data_processed <- processed_data %>%
rownames_to_column("Family") %>%
mutate(Family = factor(Family, levels = family_levels)) %>%
complete(Family = family_levels, fill = as.list(setNames(rep(0, length(orders_list)), orders_list))) %>%
distinct(Family, .keep_all = TRUE)

df_long <- data_processed %>%
pivot_longer(-Family, names_to = "Order", values_to = "abundance") %>%
mutate(Order = factor(Order, levels = orders_list),
       Family = factor(Family, levels = family_levels))

p4 <- ggplot(df_long, aes(x = Order, y = abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors, breaks = family_levels, drop = FALSE) +
    labs(title = k, x = "", y = "Relative abundance",fill = "Cazyme Family") +
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
        legend.position = "none"
    )

p4

# 5. Fungi
k <- "Fungi"
processed_data <- read.csv(paste0("./data/FigureS5-",k,"_cazyme_family_abundance.csv"),row.names = 1)
data_processed <- processed_data %>%
    rownames_to_column("Family") %>%
    mutate(Family = factor(Family, levels = family_levels)) %>%
    complete(Family = family_levels, fill = as.list(setNames(rep(0, length(orders_list)), orders_list))) %>%
    distinct(Family, .keep_all = TRUE)

df_long <- data_processed %>%
    pivot_longer(-Family, names_to = "Order", values_to = "abundance") %>%
    mutate(Order = factor(Order, levels = orders_list),
           Family = factor(Family, levels = family_levels))

p5 <- ggplot(df_long, aes(x = Order, y = abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors, breaks = family_levels, drop = FALSE) +
    labs(title = k, x = "",  y = "",fill = "Cazyme Family") +
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
        legend.position = "none"
    )

p5

# 6. Protist
k <- "Protist"
processed_data <- read.csv(paste0("./data/FigureS5-",k,"_cazyme_family_abundance.csv"),row.names = 1)
data_processed <- processed_data %>%
    rownames_to_column("Family") %>%
    mutate(Family = factor(Family, levels = family_levels)) %>%
    complete(Family = family_levels, fill = as.list(setNames(rep(0, length(orders_list)), orders_list))) %>%
    distinct(Family, .keep_all = TRUE)

df_long <- data_processed %>%
    pivot_longer(-Family, names_to = "Order", values_to = "abundance") %>%
    mutate(Order = factor(Order, levels = orders_list),
           Family = factor(Family, levels = family_levels))

p6 <- ggplot(df_long, aes(x = Order, y = abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors, breaks = family_levels, drop = FALSE) +
    labs(title = k, x = "",  y = "",fill = "Cazyme Family") +
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

p6

legend_p <- cowplot::get_legend(p6)
p456 <- cowplot::plot_grid(
    p4 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p5 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p6 + theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = "none"), 
    legend_p,
    labels = c("d", "e", "f"), label_size = 10, label_fontface = "bold",
    ncol = 4, nrow = 1
)
p456

p <- cowplot::plot_grid(
    p123 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p456 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    ncol = 1, nrow = 2
)
p
width <- 17.5
height <- 15
name <- paste0("FigureS5")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")