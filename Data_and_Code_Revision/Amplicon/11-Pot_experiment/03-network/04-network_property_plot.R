### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "04-network_property_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

treatment <- c("B", "BF", "BFP", "All")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
data_df <- NULL
for (tre in treatment) {
    network_property <- read.csv(paste0("02-get_network_property/", tre, "_netwrok_property.csv"), row.names = 1)
    network_property$Treatment <- tre
    
    data_df <- bind_rows(data_df, network_property)
}

data_df <- data_df %>%
    mutate(
        Property = case_when(
            Property == "num.edges"                  ~ "Number of edges",
            Property == "num.pos.edges"              ~ "Number of positive edges",
            Property == "num.neg.edges"              ~ "Number of negative edges",
            Property == "num.vertices"               ~ "Number of nodes",
            Property == "connectance"                ~ "Connectance",
            Property == "average.degree"             ~ "Average degree",
            Property == "average.path.length"        ~ "Average path length",
            Property == "diameter"                   ~ "Network diameter",
            Property == "edge.connectivity"          ~ "Edge connectivity",
            Property == "clustering.coefficient"     ~ "Clustering coefficient",
            Property == "no.clusters"                ~ "Number of modules",
            Property == "centralization.degree"      ~ "Degree centralization",
            Property == "centralization.betweenness" ~ "Betweenness centralization",
            Property == "centralization.closeness"   ~ "Closeness centralization",
            Property == "robustness"                 ~ "Robustness",
            TRUE ~ Property
        ),
        Treatment = factor(Treatment, levels = treatment)
    ) %>%
    filter(Property != "Edge connectivity")



p <- ggplot(data = data_df, mapping = aes(x = Treatment, y = Value)) +
    geom_bar(mapping = aes(fill = Treatment), stat = "identity", position = "stack", width = 0.68) +
    geom_text(
        mapping = aes(x = Treatment, y = Value * 1.2, label = round(Value, 3)),
        position = position_dodge(0.9),
        size = 7 / 2.835
    ) +
    facet_wrap(vars(Property), ncol = 4, scales = "free", strip.position = "top") +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top"
    )
p
name <- paste0(dir_name, "/network_property")
width <- 17
height <- 16
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------