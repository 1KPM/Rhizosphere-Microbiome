### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(2024)

# Create directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Import data ----------------------------------------------------------------
all_info <- read.csv("data/cross_network_hub_info_absolute.csv", row.names = 1)
all_info$Clade <- ifelse(grepl("bASV", row.names(all_info)), "Bacteria", 
                         ifelse(grepl("fASV", row.names(all_info)), "Fungi", "Protist"))
                         
hub_info <- all_info[all_info$roles != "Peripherals",]
# ------------------------------------------------------------------------------
all_df <- data.frame(table(all_info$Clade))
names(all_df) <- c("Kingdom", "Value")
all_df$Ratio <- all_df$Value / sum(all_df$Value)

hub_df <- data.frame(table(hub_info$Clade))
names(hub_df) <- c("Kingdom", "Value")
hub_df$Ratio <- hub_df$Value / sum(hub_df$Value)



### Get results ----------------------------------------------------------------
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
all_df$Type <- paste0("All ASVs (", sum(all_df$Value), ")")
hub_df$Type <- paste0("Hub ASVs (", sum(hub_df$Value), ")")
data_df <- bind_rows(all_df, hub_df)

p <- ggplot(data_df, aes(x = '', y = Value, fill = Kingdom)) +
    geom_bar(stat = "identity") +
    geom_text(
        aes(label = scales::percent(Ratio, 0.1)),
        color = "black", position = position_stack(vjust = 0.5),
        size = 7 / 2.835
    ) +
    coord_polar(theta = "y") +
    facet_wrap(~ Type, strip.position = "bottom", scales = "free",) +
    theme_void() +
    scale_fill_manual(values = color_manual) +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(size = 7, hjust = 0.5),
        strip.text = element_text(size = 7),
        panel.spacing = unit(-1, "lines"),
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        legend.key.size = unit(0.25, 'cm'),
        legend.position = "right"
    ) 


width <- 7.5
height <- 4
name <- paste0(dir_name, "/Figure3E")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

