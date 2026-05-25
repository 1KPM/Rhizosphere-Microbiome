### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "10-inter_network_cor_info"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

library(dplyr)
library(ggplot2)
library(ggvenn)
library(patchwork)

thresholds <- c("0.4", "0.6", "0.8", "0.72(AOD)")
type <- c("inter")

hub_list <- list()

##Summary barplot
for (thr in thresholds) {
  for (typ in type) {
    file_path <- paste0("00-data/inter_network_hub_info(cor",thr,").csv")
    
    tmp_df <- read.csv(file_path, row.names = 1)
    
    tmp_df <- tmp_df %>%
      mutate(Clade = case_when(
        substr(rownames(.), 1, 1) == "b" ~ "Bacteria",
        substr(rownames(.), 1, 1) == "f" ~ "Fungi",
        substr(rownames(.), 1, 1) == "p" ~ "Protist",
        TRUE ~ NA_character_
      )) %>%
      mutate(Clade = case_when(
        Clade == "Protist" ~ "Protists",
        TRUE ~ Clade
      ),
      Clade = factor(Clade, levels = c("Bacteria", "Fungi", "Protists")))
    
    hub_info <- tmp_df %>%
      filter(roles != "Peripherals") %>%
      group_by(Clade, .drop = FALSE) %>%
      summarise(Count = n(), .groups = "drop") %>%
      mutate(Thresholds = thr, Type = typ)
    hub_list[[paste0(thr, "_", typ)]] <- hub_info
  }
}

all_hub_info <- bind_rows(hub_list)
all_hub_info$Thresholds <- factor(all_hub_info$Thresholds, levels = thresholds)
write.csv(all_hub_info, paste0(dir_name,"/inter_hub_info_summary.csv"), row.names = FALSE)


p <- ggplot(all_hub_info, aes(x = Thresholds, y = Count, fill = Thresholds)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.3, size = 2) +
  facet_grid(~ Clade, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_x_discrete(labels = function(x) paste0("r = ", x)) + 
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Numbers of hub KOs"
  ) + 
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
    panel.spacing.y = unit(0.5, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "none"
  )

name <- paste0(dir_name, "/cor_barplot")
width <- 11
height <- 6
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


##Venn
hub_inter <- list()
for (typ in type) {
  for (thr in thresholds) {
    file_path <- paste0("00-data/inter_network_hub_info(cor",thr,").csv")
    
    tmp_df <- read.csv(file_path, row.names = 1)
    
    tmp_df <- tmp_df %>%
      mutate(Clade = case_when(
        substr(rownames(.), 1, 1) == "b" ~ "Bacteria",
        substr(rownames(.), 1, 1) == "f" ~ "Fungi",
        substr(rownames(.), 1, 1) == "p" ~ "Protist",
        TRUE ~ NA_character_
      )) %>%
      mutate(Clade = case_when(
        Clade == "Protist" ~ "Protists",
        TRUE ~ Clade
      ),
      Clade = factor(Clade, levels = c("Bacteria", "Fungi", "Protists")))
    
      hub_list <- tmp_df[tmp_df$roles != 'Peripherals',]$name
      hub_inter[[paste0("r = ", thr)]] <- hub_list
  }
}


p_venn_inter <- ggvenn(
  hub_inter,
  show_elements = FALSE,
  label_sep = "\n",
  show_percentage = FALSE,
  digits = 1,
  fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"),
  fill_alpha = 0.5,
  stroke_color = "white",
  stroke_alpha = 0.5,
  stroke_size = 0.5,
  stroke_linetype = "solid",
  set_name_color = "black",
  text_color = "black",
  set_name_size = 7 / 2.835,
  text_size = 6 / 2.835
) + labs(
  title = NULL,
  subtitle = NULL,
  x = NULL,
  y = NULL
) + theme(
  plot.title = element_text(size = 7, hjust = 0.5),
  plot.margin = margin(10, 20, 10, 20)
) + coord_cartesian(clip = "off")

name <- paste0(dir_name, "/cor_venn")
width <- 7
height <- 6
ggsave(paste0(name, ".png"), p_venn_inter, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p_venn_inter, width = width, height = height, units = "cm")


p_all <- p + p_venn_inter +
  plot_layout(
    widths = c(2, 1),
  )
p_all
name <- paste0(dir_name, "/merge")
width <- 17
height <- 6
ggsave(paste0(name, ".png"), p_all, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p_all, width = width, height = height, units = "cm")