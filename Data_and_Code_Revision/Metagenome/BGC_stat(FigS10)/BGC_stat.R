pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}


library(tidyverse)
library(cowplot)
library(ggplot2)
library(patchwork)

data <- read.csv("./data/BGC.csv")
bgc_order <- c("Polyketide", "NRP", "Saccharide", "Terpene", "RiPP", "Alkaloid", "Others")
color <- c("#F7837E","#58ACA9","#A363A2","#C49BCB","#E8E276","#019ED2","#D3D3D3")

taxonomic_levels <- c("Order")
top_n_list <- c(
  Phylum = 10,
  Class = 10,
  Order = 20,
  Family = 20
  #  Genus = 30,
  #  Species = 30
)

plot_list <- list()
for (tax_level in taxonomic_levels) {
  current_top_n <- top_n_list[tax_level]
  tax_sum <- data %>% 
    filter(!!sym(tax_level) != "Unassign") %>%
    group_by(!!sym(tax_level)) %>% 
    summarise(Total = sum(c_across(Alkaloid:Terpene))) %>%
    arrange(desc(Total))
  top_tax <- tax_sum %>% 
    slice_head(n = current_top_n) %>% 
    pull(!!sym(tax_level))
  data_processed <- data %>%
    mutate(
      Group = case_when(
        !!sym(tax_level) %in% top_tax ~ !!sym(tax_level),
        !!sym(tax_level) == "Unassign" ~ "Unassign",
        TRUE ~ "Others"
      )
    )
  
  data_agg <- data_processed %>% 
    group_by(Group) %>% 
    summarise(across(Alkaloid:Terpene, sum, na.rm = TRUE)) %>% 
    ungroup()
  order_levels <- c(top_tax, "Others", "Unassign")
  data_agg$Group <- factor(data_agg$Group, levels = order_levels)
  data_agg_long <- data_agg %>% 
    pivot_longer(
      cols = -Group,
      names_to = "Function",
      values_to = "Count"
    ) %>% 
    mutate(
      Function = factor(Function, levels = bgc_order)
    )
  p <- ggplot(data_agg_long, 
              aes(x = Group, y = Count, fill = Function)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.3) +
    labs(
      title = tax_level,
      x = NULL, y = "BGCs Numbers", fill = NULL
    ) +
    scale_fill_manual(values = color) +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.title = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                 color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "none",
      legend.box.spacing = unit(0.2, "cm"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "mm")
    ) +
    scale_x_discrete(drop = FALSE, expand = expansion(add = 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  plot_list[[tax_level]] <- p
}
width <- 20
height <- 10
name <- paste0(dir_name,"/BGCs_number")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm",paper = "a4")
ggsave(paste0(name, ".jpg"), p, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p, width = width, height = height, units = "cm")


data <- read.csv("./data/BGC_abundance.csv",row.names = 1)
bgc_order <- c("Polyketide", "NRP", "Saccharide", "Terpene", "RiPP", "Alkaloid", "Others")
order_list <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales","Others")
color <- c("#F7837E","#58ACA9","#A363A2","#C49BCB","#E8E276","#019ED2","#D3D3D3")

data_ordered <- data[bgc_order, order_list]

data_long <- data_ordered %>% 
  rownames_to_column(var = "BGC") %>% 
  pivot_longer(
    cols = -BGC,
    names_to = "Orders",
    values_to = "value"
  )
data_long$BGC <- factor(data_long$BGC, levels = bgc_order)
data_long$Orders <- factor(data_long$Orders, levels = order_list)

p1<-ggplot(data_long, aes(x = Orders, y = value, fill = BGC)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = color) +
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.position = "right",
    legend.box.spacing = unit(0.2, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

width <- 20
height <- 10
name <- paste0(dir_name,"/BGCs_abundance")
ggsave(paste0(name, ".pdf"), p1, width = width, height = height, units = "cm",paper = "a4")
ggsave(paste0(name, ".jpg"), p1, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p1, width = width, height = height, units = "cm")