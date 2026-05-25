### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "09-asv_proportion"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(patchwork)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
kingdom <- c('Bacteria', 'Fungi', 'Protist')

num <- 10
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
top_order <- read.csv("../../01-sort_data/07-top_order/top_order.csv", header = T)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
p_list <- NULL
for (typ in type) {
    first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
    second_prefix <- ifelse(typ == "all", "02", "01")
    
    dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
    hub_info <- read.csv(paste0(dir_path, '/', typ, '_network_hub_info.csv'),row.names = 1)
    
    raw_df <- merge(hub_info, core_taxonomy, by = 'row.names')
    
    for (kin in kingdom) {
        order_list <- top_order[,kin]
        tmp_df <- raw_df %>%
            filter(Clade == kin) %>%
            mutate(Order = case_when(Order %in% order_list ~ Order, TRUE ~ "Others"))
        
        n_all <- nrow(tmp_df)
        n_hub <- sum(tmp_df$roles != 'Peripherals', na.rm = TRUE)
        
        fin_df <- tmp_df %>%
            group_by(Taxonomy = Order) %>%
            summarise(
                !!paste0('All (', n_all, ')') := n(),                                       # 统计总数
                !!paste0('Hub (', n_hub, ')') := sum(roles != 'Peripherals', na.rm = TRUE)  # 统计条件数 (自动避开NA)
            ) %>%
            pivot_longer(cols = -Taxonomy, names_to = "ASV", values_to = "Value")
        
        color_df <- tibble(
            Taxonomy = c(top_order[, kin], 'Others'),
            Color = colorRampPalette(brewer.pal(9, 'Set1'))(num + 1)
        ) %>% 
            filter(Taxonomy %in% fin_df$Taxonomy)
        
        fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)
        color_manual <- color_df$Color
        
        network <- ifelse(typ == "all", "Whole network", paste0(str_to_title(typ), "kingdom network"))
        clade <- ifelse(kin == "Protist", "Protists", kin)
        
        p <- ggplot(fin_df, aes(x = ASV, y = Value, fill = Taxonomy)) +
            geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
            labs(
                title = NULL,
                subtitle = NULL,
                x = NULL,
                y = paste0("Proportion of ASVs (", network, ")")
            ) + 
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = color_manual, name = clade) +
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
                legend.box.spacing = unit(0.1,"cm"),
                legend.key.size = unit(0.25, "cm"),
                legend.position = "right"
            )
        p_list[[paste0(typ, "_", kin)]] <- p
    }
}

p <- p_list[["all_Bacteria"]] + p_list[["all_Fungi"]] + p_list[["all_Protist"]] + 
    p_list[["inter_Bacteria"]] + p_list[["inter_Fungi"]] + p_list[["inter_Protist"]] + 
    p_list[["intra_Bacteria"]] + p_list[["intra_Fungi"]] + p_list[["intra_Protist"]] + 
    plot_layout(ncol = 3, nrow = 3)

name <- paste0(dir_name, "/asv_proportion")
width <- 16
height <- 16
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------