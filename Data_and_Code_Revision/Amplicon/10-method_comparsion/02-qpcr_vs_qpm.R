### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "2-qpcr_vs_qpm"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggtext)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
qpcr_data <- read.csv("../00-rawdata/qpcr/qPCR_results.csv", header = T)
absolute_abundance <- read.csv("../01-sort_data/03-absolute_abundance/all_sample_copies.csv", header = T)
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
qpcr_df <- qpcr_data %>%
    pivot_longer(
        cols = -SampleID,
        names_to = "Kingdom",
        values_to = "Value"
    ) %>%
    mutate(
        Value = round(Value) * 200 * 100,
        Method = "qPCR"
    )

absolute_df <- absolute_abundance %>%
    filter(FileID %in% qpcr_df$SampleID) %>%
    select(SampleID = FileID,
           Bacteria = Copies_16S,
           Fungi = Copies_ITS,
           Protists = Copies_Protist) %>%
    pivot_longer(
        cols = -SampleID,
        names_to = "Kingdom",
        values_to = "Value"
    ) %>%
    mutate(Method = "Sequencing")

data_df <- rbind(qpcr_df %>% filter(SampleID %in% absolute_df$SampleID), absolute_df)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
cor_df <- data_df %>%
    pivot_wider(names_from = Method, values_from = Value)

p <- ggplot(data = cor_df, mapping = aes(x = log10(qPCR), y = log10(Sequencing))) +
    geom_point(mapping = aes(fill = Kingdom), size = 1.5, shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
    geom_smooth(method = "lm", color = "firebrick", fill = "grey80", linewidth = 0.5, alpha = 0.4) +
    facet_wrap(vars(Kingdom), scales = "free", strip.position = "top") +
    stat_cor(method = "pearson", size = 6 / 2.835, color = "black", label.x.npc = "left", label.y.npc = "top") +
    scale_color_manual(values = color_manual)+
    labs(
        title = NULL,
        subtitle = NULL,
        x = "Absolute abundance by qPCR (log<sub>10</sub>copies g<sup>-1</sup> soil)",
        y = "Absolute abundance by sequencing <br> (log<sub>10</sub>copies g<sup>-1</sup> soil)"
    ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
        axis.title.x = element_markdown(size = 7),
        axis.title.y = element_markdown(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.2, "cm"),
        legend.position = "none" 
    )

name <- paste0(dir_name, "/compare_correlation")
width_cor <- 13.5
height_cor <- 5
ggsave(paste0(name, ".pdf"), p, width = width_cor, height = height_cor, units = "cm")


# ------------------------------------------------------------------------------
