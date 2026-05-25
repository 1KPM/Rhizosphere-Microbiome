pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

##nr
nr_data <- read.csv("./data/FigureS4-nr-pie.csv")

nr_number <-nr_data %>%
  filter(Name != "Total") %>% 
  mutate(
    Group = factor(Name, levels = c("Bacteria","Fungi","Protist","Others","Unassign")),
    ratio = round(Number / sum(Number), 6)*100,
    count = ifelse(ratio < 1, 1, ratio)
  )
p1 <- ggpubr::ggpie(
  data = nr_number,
  x = "count",
  #  label = c("36.62%","0.13%","0.03%","0.17%","63.06%"),
  label = c("36.62%","","","","63.06%"),
  fill = "Name",
  lab.pos = "in",
  color = "white",
  lab.font = 2,
  palette = c("Bacteria" = "#7FC97F","Fungi" = "#BEAED4", "Protist" = "#FDC086","Others" = "#80B1D3","Unassign" = "grey80")) +
  labs(title = "Gene Number") +
  theme(
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 1)
p1
nr_abundance <-nr_data %>%
  filter(Name != "Total") %>% 
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
  lab.font = 2,
  color = "white",
  palette = c("Bacteria" = "#7FC97F","Fungi" = "#BEAED4", "Protist" = "#FDC086","Others" = "#80B1D3","Unassign" = "grey80")) +
  labs(title = "Gene Abundance") +
  theme(
     text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(0.25, "cm"),
    legend.title = element_text(size = 7),
    aspect.ratio = 1)
p2
##kegg
kegg_data <- read.csv("./data/FigureS4-kegg-pie.csv")

kegg_number <- kegg_data %>%
  filter(!Name %in% c("Total")) %>%
  mutate(
    Category = case_when(
      Name == "Unassign" ~ "Unassign",
      Name == "Others" ~ "Others",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Category) %>%
  mutate(
    Group = case_when(
      Category == "Unassign" ~ "Unassign",
      Category == "Others" ~ "Others",
      row_number() <= 10 ~ Name,
      TRUE ~ "Others_New"
    )
  ) %>%
  ungroup() %>%
  mutate(Group = if_else(Group == "Others_New", "Others", Group)) %>%
  group_by(Group) %>%
  summarise(Number = sum(Number)) %>%
  mutate(
    ratio = round(Number / sum(Number) * 100, 6),
    #    label = sprintf("%.2f%%", ratio)
  ) %>%
  mutate(
    sort_group = case_when(
      Group == "Unassign" ~ 3,
      Group == "Others" ~ 2,
      TRUE ~ 1
    )
  ) %>%
  arrange(sort_group, desc(Number)) %>%
  mutate(
    Group = factor(Group, levels = c(
      setdiff(unique(Group), c("Others", "Unassign")),
      "Others",
      "Unassign"
    ))
  ) %>%
  select(-sort_group)


named_palette <- c(
  "Carbohydrate metabolism" = "#E41A1C",
  "Amino acid metabolism" = "#A24057",
  "Energy metabolism" = "#606692",
  "Cellular community" = "#3A85A8",
  "Lipid metabolism" = "#9C509B",
  "Biosynthesis of other secondary metabolites" = "#E1C62F",
  "Cell growth and death" = "#CC6A6F",
  "Glycan biosynthesis and metabolism" = "#FFF830",
  "Folding, sorting and degradation" = "#BF862B",
  "Cell motility" = "#EB7AA9",
  "Others" = "#CCCCCC",  # 浅灰色
  "Unassign" = "#999999"  # 中灰色
)

p3 <- ggpubr::ggpie(
  data = kegg_number,
  x = "Number",
  #  label = c("7.92%","6.49%","3.31%","2.72%","2.06%","1.52%","1.34%","1.09%","0.77%","0.56%","33.78%","38.44%"),
  label = c("7.92%","6.49%","3.31%","2.72%","2.06%","","","","","","33.78%","38.44%"),
  fill = "Group",
  lab.pos = "in",
  color = "white",
  lab.font = 2,
)+
  scale_fill_manual(values = named_palette) +
  #  labs(title = "Gene Number") +
  theme(
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 1)
p3
kegg_abundance <- kegg_data %>%
  filter(!Name %in% c("Total")) %>%
  mutate(
    Category = case_when(
      Name == "Unassign" ~ "Unassign",
      Name == "Others" ~ "Others",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Category) %>%
  mutate(
    Group = case_when(
      Category == "Unassign" ~ "Unassign",
      Category == "Others" ~ "Others",
      row_number() <= 10 ~ Name,
      TRUE ~ "Others_New"
    )
  ) %>%
  ungroup() %>%
  mutate(Group = if_else(Group == "Others_New", "Others", Group)) %>%
  group_by(Group) %>%
  summarise(Abundance = sum(Abundance)) %>%
  mutate(
    ratio = round(Abundance / sum(Abundance) * 100, 6),
    #    label = sprintf("%.2f%%", ratio)
  ) %>%
  mutate(
    sort_group = case_when(
      Group == "Unassign" ~ 3,
      Group == "Others" ~ 2,
      TRUE ~ 1
    )
  ) %>%
  arrange(sort_group, desc(Abundance)) %>%
  mutate(
    Group = factor(Group, levels = c(
      setdiff(unique(Group), c("Others", "Unassign")),
      "Others",
      "Unassign"
    ))
  ) %>%
  select(-sort_group)


length <- length(kegg_abundance$Name)
my_palette <- c(colorRampPalette(brewer.pal(9, 'Set1'))(12)[1:11],'grey80')

p4 <- ggpubr::ggpie(
  data = kegg_abundance,
  x = "Abundance",
  #  label = c("7.88%","6.69%","3.80%","2.71%","1.89%","1.17%","1.09%","1.06%","0.88%","0.58%","32.65%","39.58%"),
  label = c("7.88%","6.69%","3.80%","2.71%","","","","","","","32.65%","39.58%"),
  fill = "Group",
  lab.pos = "in",
  lab.font = 2,
  color = "white",
)+
  scale_fill_manual(values = named_palette) +
  #  labs(title = "Abundance") +
    theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7,hjust = 0.5),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 7),
        aspect.ratio = 1)
p4

##core
core_data <- read.csv("./data/FigureS4-core-KO.csv")

core_number <-core_data %>%
  filter(Name != "Total") %>% 
  mutate(
    Group = factor(Name),
    ratio = round(Number / sum(Number), 6)*100,
    label = sprintf("%.2f%%", ratio)
  )
p5 <- ggpubr::ggpie(
  data = core_number,
  x = "ratio",
  label = "label",
  fill = "Name",
  lab.pos = "in",
  lab.font = 2,
  color = "white",
  palette = c("#6EB9C3","#C98B88")) +
  #  labs(title = "Gene Number") +
  theme(
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position = "none",
    aspect.ratio = 1)

p5

core_abundance <-core_data %>%
  filter(Name != "Total") %>% 
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
  color = "white",
  lab.font = 2,
  palette = c("#6EB9C3","#C98B88")) +
  #  labs(title = "Abundance") +
    theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7,hjust = 0.5),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 7),
        aspect.ratio = 1)
p6

p <- cowplot::plot_grid(
    p1 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p2 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p3 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p4 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p5 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    p6 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    labels = c("a", " ", "b", " ", "c", " "), label_size = 10, label_fontface = "bold",
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 2, nrow = 3, rel_heights = c(1, 1, 1) 
)
width <- 22
height <- 15
name <- paste0("FigureS4")

ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

