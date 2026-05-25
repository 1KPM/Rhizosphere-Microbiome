pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(VennDiagram)
library(ggplotify)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

## 1.NR Plot
nr_data <- read.csv("./data/nr-pie.csv")

nr_number <-nr_data %>%
  mutate(
    Group = factor(Name, levels = c("Bacteria","Fungi","Protist","Others","Unassign")),
    ratio = round(Number / sum(Number), 4)*100,
    count = ifelse(ratio < 1, 1, ratio)
  )
p1 <- ggpubr::ggpie(
  data = nr_number,
  x = "count",
  #  label = c("36.62%","0.13%","0.03%","0.17%","63.06%"),
  label = c("36.62%","","","","63.06%"),
  fill = "Name",
  lab.pos = "in",
  lab.font = 3,
  palette = c("Bacteria" = "#7FC97F","Fungi" = "#BEAED4", "Protist" = "#FDC086","Others" = "#80B1D3","Unassign" = "grey80")) +
  labs(title = "Gene Number") +
  theme(
    plot.background   = element_blank(),
    panel.background  = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(size = 8,hjust = 0.5, face = "bold"),
    legend.position = "none",
    aspect.ratio = 1,
    plot.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0), 
    legend.box.margin = margin(0,0,0,0)
    )
#width <- 5 height <- 5
#ggsave("result/nr_number.png", p1, width = width, height = height, dpi = 600, units = "cm")
#ggsave("result/nr_number.pdf", p1, width = width, height = height, units = "cm")

nr_abundance <-nr_data %>%
  mutate(
    Group = factor(Name, levels = c("Bacteria","Fungi","Protist","Others","Unassign")),
    ratio = round(Abundance / sum(Abundance), 6)*100,
    count = ifelse(ratio < 1, 1, ratio)
  )

p2 <- ggpubr::ggpie(
  data = nr_abundance,
  x = "count",
  #  label = c("43.65%","0.10%","0.03%","0.27%","55.94%"),
  label = c("43.65%","","","","55.94%"),
  fill = "Name",
  lab.pos = "in",
  lab.font = 3,
  palette = c("Bacteria" = "#7FC97F","Fungi" = "#BEAED4", "Protist" = "#FDC086","Others" = "#80B1D3","Unassign" = "grey80")) +
  labs(title = "Abundance") +
  theme(
    plot.background   = element_blank(),
    panel.background  = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(size = 8,hjust = 0.5, face = "bold"),
    legend.position = "right",
    aspect.ratio = 1,
    plot.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0), 
    legend.box.margin = margin(0,0,0,0)
  )

#ggsave("result/nr_abundance.png", p2, width = width, height = height, dpi = 600, units = "cm")
#ggsave("result/nr_abundance.pdf", p2, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------

## 2.KEGG vs eggNOG Plot
kegg_cog_data <-read.csv("data/kegg-cog.csv")

kegg_cog_number <- kegg_cog_data %>% select(region,Number)

number_regions <- list(
  n1 = kegg_cog_number$Number[kegg_cog_number$region == "eggnog"],
  n2 = kegg_cog_number$Number[kegg_cog_number$region == "kegg"],
  n12 = kegg_cog_number$Number[kegg_cog_number$region == "eggnog-kegg"]
)
p3 <- draw.pairwise.venn(
  area1 = number_regions$n1 + number_regions$n12,
  area2 = number_regions$n2 + number_regions$n12,
  cross.area = number_regions$n12,
  category = c("COG", "KEGG"),
  fill = c("#2E86AB", "#F18F01"),
  alpha = 0.6,
  lty = "solid",
  lwd = 2,
  col = "black",
  cex = 0.7,
  cat.cex = 0.8,
  cat.col = c("#2E86AB", "#F18F01"),
  cat.dist = c(0.05, 0.05),
  cat.pos = c(0, 0),
  print.mode = "percent",
  sigdigs = 3,
  margin = 0,
  ind = TRUE,
  euler.d = TRUE,
  scaled = FALSE,
  direct.area = FALSE
)

#png("result/COG-kegg_number.png", width = width, height = height, units = "cm", res = 600); print(p3); dev.off()
#pdf("result/COG-kegg_number.pdf", width = width / 2.54, height = height / 2.54); print(p3); dev.off()

kegg_cog_abun <- kegg_cog_data %>% select(region,Abundance)

abun_regions <- list(
  n1 = kegg_cog_abun$Abundance[kegg_cog_abun$region == "eggnog"],
  n2 = kegg_cog_abun$Abundance[kegg_cog_abun$region == "kegg"],
  n12 = kegg_cog_abun$Abundance[kegg_cog_abun$region == "eggnog-kegg"]
)
unassign_count <- kegg_cog_abun$Abundance[kegg_cog_abun$region == "unassign"]

p4 <- draw.pairwise.venn(
  area1 = abun_regions$n1 + abun_regions$n12,
  area2 = abun_regions$n2 + abun_regions$n12,
  cross.area = abun_regions$n12,
  category = c("COG", "KEGG"),
  fill = c("#2E86AB", "#F18F01"),
  alpha = 0.6,
  lty = "solid",
  lwd = 2,
  col = "black",
  cex = 0.7,
  cat.cex = 0.8,
  cat.col = c("#2E86AB", "#F18F01"),
  cat.dist = c(0.05, 0.05),
  cat.pos = c(0, 0),
  print.mode = "percent",
  sigdigs = 3,
  margin = 0,
  ind = TRUE,
  euler.d = TRUE,
  scaled = FALSE,
  direct.area = FALSE
)

#png("result/COG-kegg_abundance.png", width = width, height = height, units = "cm", res = 600); print(p4); dev.off()
#pdf("result/COG-kegg_abundance.pdf", width = width / 2.54, height = height / 2.54); print(p4); dev.off()
# ------------------------------------------------------------------------------
##core
core_data <- read.csv("./data/core-pie.csv")

core_number <-core_data %>%
  mutate(
    Group = factor(Name),
    ratio = round(Number / sum(Number), 4)*100,
    label = sprintf("%.2f%%", ratio)
  )
p5 <- ggpubr::ggpie(
  data = core_number,
  x = "ratio",
  label = "label",
  fill = "Name",
  lab.pos = "in",
  lab.font = 3,
  palette = c("#6EB9C3","#C98B88")) +
  theme(
    plot.background   = element_blank(),
    panel.background  = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(size = 8,hjust = 0.5, face = "bold"),
    legend.position = "none",
    aspect.ratio = 1,
    plot.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0), 
    legend.box.margin = margin(0,0,0,0)
  )

ggsave("result/core_number.png", p5, width = 5, height = 5, dpi = 600, units = "cm")
ggsave("result/core_number.pdf", p5, width = 5, height = 5, units = "cm")

core_abundance <-core_data %>%
  mutate(
    Group = factor(Name),
    ratio = round(Abundance / sum(Abundance), 6)*100,
    label = sprintf("%.2f%%", ratio)
  )
p6 <- ggpubr::ggpie(
  data = core_abundance,
  x = "ratio",
  label = "label",
  fill = "Name",
  lab.pos = "in",
  lab.font = 3,
  palette = c("#6EB9C3","#C98B88")) +
  labs(fill = "Core KO") + 
  theme(
    plot.background   = element_blank(),
    panel.background  = element_blank(),
    legend.background = element_blank(),
    plot.title = element_text(size = 8,hjust = 0.5, face = "bold"),
    legend.position = "right",
    aspect.ratio = 1,
    plot.margin = margin(0,0,0,0),
    legend.margin = margin(0,0,0,0), 
    legend.box.margin = margin(0,0,0,0)
  )

ggsave("result/core_abundance.png", p6, width = 5, height = 5, dpi = 600, units = "cm")
ggsave("result/core_abundance.pdf", p6, width = 5, height = 5, units = "cm")


p3_g <- as.ggplot(grid::grid.grabExpr(grid::grid.draw(p3))) + coord_equal()
p4_g <- as.ggplot(grid::grid.grabExpr(grid::grid.draw(p4))) + coord_equal()

final_plot <- p1 + p2 + p3_g + p4_g + p5 + p6 +
  plot_layout(ncol = 2, byrow = TRUE) +
  plot_annotation(tag_levels = "a")


width <- 12.5
height <- 15.5
name <- paste0("result/final_plot")
ggsave(paste0(name, ".png"), final_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), final_plot, width = width, height = height, units = "cm",paper = "a4")

