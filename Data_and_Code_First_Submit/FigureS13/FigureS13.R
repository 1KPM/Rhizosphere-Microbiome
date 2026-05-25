### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[['path']])
setwd(pwd)

# Set Seed
set.seed(2024)

# Import Packages
library(reshape2)
library(ggplot2)
# Create Directory
dir_name <- 'results'
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------



### Get Results ----------------------------------------------------------------
# 1. within
robustness_df <- read.csv('data/within_robustness_by_asv_level.csv')
robustness_kingdom <- read.csv('data/within_robustness_by_asv_group_kingdom.csv')
names(robustness_df) <- c('Number', 'Core ASV', 'Random ASV')
names(robustness_kingdom) <- c("Clade", "Number", "Target", "Random")
robustness_kingdom$Group <- paste0(robustness_kingdom$Clade, " (", robustness_kingdom$Number, ")")
fin_df <- melt(robustness_df, id.vars = 'Number', variable.name = 'ASV', value.name = 'Robustness')
p1 <- ggplot(fin_df, aes(x = Number, y = Robustness, color = ASV)) + 
    geom_line() +
    geom_point(size = 0.5) +
    labs(
      x = 'Removed ASV number',
      y = 'Within-kingdom network roubustness'
    ) +
    theme_bw() + 
    scale_color_manual(values = c('#c01f25','#57657c'), name = "ASV type") +
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
 
          legend.position = 'right')

kingdom_df <- melt(robustness_kingdom[c("Group", "Random", "Target")], id.vars = "Group", variable.name = 'Remove', value.name = 'Robustness')

p2 <- ggplot(
    data = kingdom_df, 
    mapping = aes(x = Group, y = Robustness, fill = Remove)
    ) +
    labs(
        x = "",
        y = "Within-kingdom network roubustness"
    ) + 
    geom_bar(
        stat = "identity",
        position = "dodge",
        width = 0.68
    ) +
    scale_fill_manual(values = c('#57657c', '#c01f25'), name = "Remove type") +
    theme_bw() +
    theme(
    plot.title = element_text(size = 7, color = "black", hjust = 0.5),
    plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
    legend.title = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 6,  color = "black"),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.25, "cm"),
    )

# 2. cross
robustness_df <- read.csv('data/cross_robustness_by_asv_level.csv')
robustness_kingdom <- read.csv('data/cross_robustness_by_asv_group_kingdom.csv')
names(robustness_df) <- c('Number', 'Core ASV', 'Random ASV')
names(robustness_kingdom) <- c("Clade", "Number", "Target", "Random")
robustness_kingdom$Group <- paste0(robustness_kingdom$Clade, " (", robustness_kingdom$Number, ")")
fin_df <- melt(robustness_df, id.vars = 'Number', variable.name = 'ASV', value.name = 'Robustness')
p3 <- ggplot(fin_df, aes(x = Number, y = Robustness, color = ASV)) + 
    geom_line() +
    geom_point(size = 0.5) +
    labs(
        x = 'Removed ASV number',
        y = 'Cross-kingdom network roubustness'
    ) +
    theme_bw() + 
    scale_color_manual(values = c('#c01f25','#57657c'), name = "ASV type") +
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          
          legend.position = 'right')

kingdom_df <- melt(robustness_kingdom[c("Group", "Random", "Target")], id.vars = "Group", variable.name = 'Remove', value.name = 'Robustness')

p4 <- ggplot(
    data = kingdom_df, 
    mapping = aes(x = Group, y = Robustness, fill = Remove)
) +
    labs(
        x = "",
        y = "Cross-kingdom network roubustness"
    ) + 
    geom_bar(
        stat = "identity",
        position = "dodge",
        width = 0.68
    ) +
    scale_fill_manual(values = c('#57657c', '#c01f25'), name = "Remove type") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 7, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 6,  color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.25, "cm"),
    )

p <- cowplot::plot_grid(
    p1 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p2 + theme(plot.margin = margin(0, 0, 0, 0.5, "cm")), 
    p3 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p4 + theme(plot.margin = margin(0, 0, 0, 0.5, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 2, nrow = 2
)
width <- 17
height <- 12
name <- paste0(dir_name, "/FigureS13")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

# ------------------------------------------------------------------------------

