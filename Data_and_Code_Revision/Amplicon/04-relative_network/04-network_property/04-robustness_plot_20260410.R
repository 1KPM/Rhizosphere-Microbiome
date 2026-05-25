### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "04-robustness_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(reshape2)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
color_manual <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(2)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. Core ASV 累积去除后的Robustness
for (typ in type) {
    file_path <- paste0("03-robustness/", typ, "_robustness_by_cumulative_asv.csv")
    
    tmp_df <- read.csv(file_path, header = T)
    names(tmp_df) <- c("Counts", "Target", "Random")
    
    data_df <- melt(tmp_df, id.vars = 'Counts', variable.name = 'Type', value.name = 'Robustness')
    data_df$Network <- ifelse(typ == "all", "Whole network", ifelse(typ == "inter", "Interkingdom network", "Intrakingdom network"))
    
    p <- ggplot(data_df, aes(x = Counts, y = Robustness, color = Type)) + 
        geom_line() +
        geom_point(size = 0.5) +
        labs(
            x = 'Removed ASV number',
            y = 'Roubustness'
        ) +
        theme_bw() +
        scale_color_manual(values = c("#1B9E77", "#999999"), name = NULL) +
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.text = element_text(size = 6, color = 'black'), 
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.margin = margin(t = 0, r = 0, b = -5, l = 0, unit = "pt"),
              legend.box.margin = margin(b = -10, unit = "pt"),
              legend.position = 'top')
    
    name <- paste0(dir_name, "/robustness_by_cumulative_asv_",typ)
    width <- 5.3
    height <- 3.8
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}


# 2. Core ASV 按Kingdom去除后的Robustness
for (typ in type) {
    file_path <- paste0("03-robustness/", typ, "_robustness_by_kingdom.csv")
    
    tmp_df <- read.csv(file_path, header = T)
    names(tmp_df) <- c("Kingdom", "Counts", "Target", "Random")
    tmp_df <- tmp_df[tmp_df$Counts != 0,]
    tmp_df <- tmp_df %>%
        mutate(Kingdom = case_when(Kingdom == "Protist" ~ "Protists", TRUE ~ Kingdom))
    
    
    tmp_df$Group <- paste0(tmp_df$Kingdom, " (", tmp_df$Counts, ")")
    
    data_df <- melt(tmp_df[c("Group", "Random", "Target")], id.vars = "Group", variable.name = 'Type', value.name = 'Robustness')
    data_df$Network <- ifelse(typ == "all", "Whole network", ifelse(typ == "inter", "Interkingdom network", "Intrakingdom network"))
    
    
    name <- paste0(dir_name, "/robustness_by_kingdom_", typ)
    max_value <- ifelse(typ == "inter", 0.96, 1.01)
    min_value <- ifelse(typ == "inter", 0.88, 0.9)

    p <- ggplot(
        data = data_df, 
        mapping = aes(x = Group, y = Robustness, fill = Type)
    ) +
        labs(
            x = "",
            y = "Roubustness"
        ) + 
        geom_bar(
            stat = "identity",
            position = "dodge",
            width = 0.68
        ) +
        geom_text(
            data = data_df,
            mapping = aes(x = Group, y = Robustness * 1.01, label = round(Robustness, 3)),
            position = position_dodge(0.9),
            size = 6 / 2.835
        ) +
        coord_cartesian(ylim = c(min_value, max_value))+
        scale_fill_manual(values = c("#999999", "#1B9E77"), name = NULL) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 7, color = "black", hjust = 0.5),
            plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5),
            legend.title = element_text(size = 7, color = "black"),
            legend.text = element_text(size = 6,  color = "black"),
            axis.title = element_text(size = 7, color = "black"),
            axis.text = element_text(size = 6, color = "black"),
            # axis.text.x = element_text(angle = 45, hjust = 1),
            # axis.text.x = element_blank(),
            legend.key.size = unit(0.25, "cm"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.margin = margin(t = 0, r = 0, b = -5, l = 0, unit = "pt"),
            legend.box.margin = margin(b = -10, unit = "pt"),
            legend.position = "none"
        )
    p
    width <- 5.3
    height <- 3.8
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}


# ------------------------------------------------------------------------------