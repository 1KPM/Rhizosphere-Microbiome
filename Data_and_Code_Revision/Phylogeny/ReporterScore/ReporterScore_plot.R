pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(patchwork)

Order <- c("Rosales", "Lamiales", "Malpighiales", "Sapindales", "Gentianales", 
           "Asparagales","Malvales", "Myrtales","Arecales", "Others")

data<-read.csv("./data/Bacteria_Fabales_vs_all_orders_pathway_list1.csv") 
kegg_info <- read.csv("./data/kegg_pathway_for_1KPM.csv")
kegg_info$pathway <- sub("\\s*\\[.*", "", kegg_info$pathway)
kegg_info <- kegg_info %>% select(pathway, L3) %>% distinct()

data_filter <- data %>% 
  filter(p.adjust < 0.05 & abs(ReporterScore) >1.96) %>%
  filter(ID %in% kegg_info$L3)


data_enrich <- data_filter %>% 
  filter(ReporterScore > 0) %>%
  group_by(Description) %>%
  filter(n() >= 3) %>% 
  ungroup() %>%
  separate(GroupComparison, into = c("Class", "Group"), sep = "\\s+vs\\s+")

data_enrich$Group <- factor(data_enrich$Group, levels = Order)

enriched_p <- ggplot(data_enrich, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) + 
  geom_point(aes(size = ReporterScore, fill = -log10(p.adjust)), 
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(name = "-log10(p.adjust)",
                      low = "#FED976", high = "#E31A1C") +
  scale_size_continuous(name = "Reporter Score",
                        range = c(2, 7),
                        breaks = c(2, 4, 6)) +
  labs(title = "Enriched pathways", x = "-log10(p.adjust)", y = "") +
  facet_grid(~ Class + Group, switch = "x") +
  theme_bw(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "#F0F0F0"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5)
  )

data_depleted <- data_filter %>% 
  filter(ReporterScore < 0) %>%
  group_by(Description) %>%
  filter(n() >= 3) %>% 
  ungroup() %>%
  mutate(ReporterScore = ReporterScore * -1) %>%
  separate(GroupComparison, into = c("Class", "Group"), sep = "\\s+vs\\s+")


data_depleted$Group <- factor(data_depleted$Group, levels = Order)

depleted_p <- ggplot(data_depleted, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) + 
  geom_point(aes(size = ReporterScore, fill = -log10(p.adjust)), 
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(name = "-log10(p.adjust)",
                      low = "lightblue", high = "#1F78B4",guide = "none") +
  scale_size_continuous(name = "Reporter Score",
                        range = c(2, 7),
                        breaks = c(2, 4, 6),guide = "none") +
  labs(title = "Depleted pathways", x = "-log10(p.adjust)", y = "") +
  facet_grid(~ Class + Group, switch = "x") +
  theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "#F0F0F0"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5)
  )

combined_plot <- enriched_p / depleted_p + 
  plot_layout(ncol = 1, heights = c(3, 5)) +
  plot_annotation(tag_levels = 'A')


width <- 20
height <- 26
name <- "./Bacteria_GRSA"

ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), combined_plot, width = width, height = height, units = "cm")



##Fungi
data<-read.csv("./data/Fungi_Fabales_vs_all_orders_pathway_list1.csv") 
kegg_info <- read.csv("./data/kegg_pathway_for_1KPM.csv")
kegg_info$pathway <- sub("\\s*\\[.*", "", kegg_info$pathway)
kegg_info <- kegg_info %>% select(pathway, L3) %>% distinct()

data_filter <- data %>% 
  filter(p.adjust < 0.05 & abs(ReporterScore) >1.96) %>%
  filter(ID %in% kegg_info$L3)


data_enrich <- data_filter %>% 
  filter(ReporterScore > 0) %>%
  group_by(Description) %>%
  filter(n() >= 3) %>% 
  ungroup() %>%
  separate(GroupComparison, into = c("Class", "Group"), sep = "\\s+vs\\s+")

data_enrich$Group <- factor(data_enrich$Group, levels = Order)

enriched_p <- ggplot(data_enrich, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) + 
  geom_point(aes(size = ReporterScore, fill = -log10(p.adjust)), 
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(name = "-log10(p.adjust)",
                      low = "#FED976", high = "#E31A1C") +
  scale_size_continuous(name = "Reporter Score",
                        range = c(2, 7),
                        breaks = c(2, 4, 6)) +
  labs(title = "Enriched pathways", x = "-log10(p.adjust)", y = "") +
  facet_grid(~ Class + Group, switch = "x") +
  theme_bw(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "#F0F0F0"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5)
  )

data_depleted <- data_filter %>% 
  filter(ReporterScore < 0) %>%
  group_by(Description) %>%
  filter(n() >= 3) %>% 
  ungroup() %>%
  mutate(ReporterScore = ReporterScore * -1) %>%
  separate(GroupComparison, into = c("Class", "Group"), sep = "\\s+vs\\s+")


data_depleted$Group <- factor(data_depleted$Group, levels = Order)

depleted_p <- ggplot(data_depleted, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) + 
  geom_point(aes(size = ReporterScore, fill = -log10(p.adjust)), 
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(name = "-log10(p.adjust)",
                      low = "lightblue", high = "#1F78B4",guide = "none") +
  scale_size_continuous(name = "Reporter Score",
                        range = c(2, 7),
                        breaks = c(2, 4, 6),guide = "none") +
  labs(title = "Depleted pathways", x = "-log10(p.adjust)", y = "") +
  facet_grid(~ Class + Group, switch = "x") +
  theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "#F0F0F0"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5)
  )

combined_plot <- enriched_p / depleted_p + 
  plot_layout(ncol = 1, heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')


width <- 20
height <- 26
name <- "./Fungi_GRSA"

ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), combined_plot, width = width, height = height, units = "cm")


###################################################################################
##Protist-无结果
data<-read.csv("./data/Protist_Fabales_vs_all_orders_pathway_list1.csv") 
kegg_info <- read.csv("./data/kegg_pathway_for_1KPM.csv")
kegg_info$pathway <- sub("\\s*\\[.*", "", kegg_info$pathway)
kegg_info <- kegg_info %>% select(pathway, L3) %>% distinct()

data_filter <- data %>% 
  filter(p.adjust < 0.05 & abs(ReporterScore) >1.64) %>%
  filter(ID %in% kegg_info$L3)

data_enrich <- data_filter %>% 
  filter(ReporterScore > 0) %>%
  group_by(Description) %>%
  filter(n() >= 3) %>% 
  ungroup() %>%
  separate(GroupComparison, into = c("Class", "Group"), sep = "\\s+vs\\s+")

data_enrich$Group <- factor(data_enrich$Group, levels = Order)

enriched_p <- ggplot(data_enrich, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) + 
  geom_point(aes(size = ReporterScore, fill = -log10(p.adjust)), 
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(name = "-log10(p.adjust)",
                      low = "#FED976", high = "#E31A1C") +
  scale_size_continuous(name = "Reporter Score",
                        range = c(2, 7),
                        breaks = c(2, 4, 6)) +
  labs(title = "Enriched pathways", x = "-log10(p.adjust)", y = "") +
  facet_grid(~ Class + Group, switch = "x") +
  theme_bw(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "#F0F0F0"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5)
  )

data_depleted <- data_filter %>% 
  filter(ReporterScore < 0) %>%
  group_by(Description) %>%
  filter(n() >= 3) %>% 
  ungroup() %>%
  mutate(ReporterScore = ReporterScore * -1) %>%
  separate(GroupComparison, into = c("Class", "Group"), sep = "\\s+vs\\s+")


data_depleted$Group <- factor(data_depleted$Group, levels = Order)

depleted_p <- ggplot(data_depleted, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) + 
  geom_point(aes(size = ReporterScore, fill = -log10(p.adjust)), 
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_gradient(name = "-log10(p.adjust)",
                      low = "lightblue", high = "#1F78B4",guide = "none") +
  scale_size_continuous(name = "Reporter Score",
                        range = c(2, 7),
                        breaks = c(2, 4, 6),guide = "none") +
  labs(title = "Depleted pathways", x = "-log10(p.adjust)", y = "") +
  facet_grid(~ Class + Group, switch = "x") +
  theme_bw(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    strip.text.x = element_text(size = 6),
    strip.background = element_rect(fill = "#F0F0F0"),
    legend.position = "none",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5)
  )

combined_plot <- enriched_p / depleted_p + 
  plot_layout(ncol = 1, heights = c(2, 1)) +
  plot_annotation(tag_levels = 'A')


width <- 20
height <- 26
name <- "./Protist_GRSA"

ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), combined_plot, width = width, height = height, units = "cm")