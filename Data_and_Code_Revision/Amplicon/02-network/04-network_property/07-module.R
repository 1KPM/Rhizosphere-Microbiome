### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "07-module"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggh4x)

### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
kingdom <- c('Bacteria', 'Fungi', 'Protist')
num <- 10
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. Module Correlation
all_data_df <- all_main_module <- all_main_edge <- NULL

for (typ in type) {
    first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
    second_prefix <- ifelse(typ == "all", "02", "01")
    network <- ifelse(typ == "all", "Whole", ifelse(typ == "inter", "Interkingdom", "Intrakingdom"))
    
    dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
    
    node <- read.csv(paste0(dir_path, '/', typ, '_maptree_node.csv'))
    edge <- read.csv(paste0(dir_path, '/', typ, '_maptree_edge.csv'))
    module <- read.csv(paste0(dir_path, '/', typ, '_maptree_module.csv'))
    
    module_freq <- as.data.frame(table(module$group)) %>%
        arrange(-Freq) %>%
        filter(Var1 != "mother_no" & Freq >= nrow(module) * 0.1)
        
    main_module <- module %>%
        select(-degree) %>%
        filter(group %in% module_freq$Var1) %>%
        mutate(
            Network = paste0(network, " network"),
            group = gsub("model_", "M", group)
            ) 
    
    main_edge <- main_module %>%
        rename(OTU_1 = ID, Module_OTU_1 = group) %>%
        inner_join(edge, by = "OTU_1") %>%
        inner_join(
            main_module %>%
                rename(OTU_2 = ID, Module_OTU_2 = group), by = "OTU_2"
            ) %>%
        select(X1, Y1, OTU_1, Module_OTU_1, X2, Y2, OTU_2, Module_OTU_2, weight, Correlation)
    
    data_df <- main_edge %>%
        mutate(
            Module_1 = pmin(Module_OTU_1, Module_OTU_2),
            Module_2 = pmax(Module_OTU_1, Module_OTU_2)
        ) %>%
        count(Module_1, Module_2, Correlation, name = "Count") %>%
        rename(
            Module_OTU_1 = Module_1, 
            Module_OTU_2 = Module_2
        ) %>%
        mutate(
            Network = paste0(network, " network"),
            Correlation = case_when(Correlation == "-" ~ "Negative", TRUE ~ "Positive")
        )
    
    all_data_df <- bind_rows(all_data_df, data_df)
    all_main_module <- bind_rows(all_main_module, main_module)
    all_main_edge <- bind_rows(all_main_edge, main_edge %>% mutate(Network = paste0(network, " network")))
}

data_df$Network <- factor(data_df$Network, levels = c("Whole network", "Interkingdom network", "Intrakingdom network"))

p <- all_data_df %>%
    mutate(Interaction_Pair = paste(Module_OTU_1, Module_OTU_2, sep = "-")) %>%
    
    ggplot(aes(x = Interaction_Pair, y = Count, fill = Correlation)) +
    geom_text(
        mapping = aes(x = Interaction_Pair, y = Count + 10000, label = Count),
        position = position_dodge(0.9),
        size = 6 / 2.835
    ) +
    facet_grid(cols = vars(Network), scales = "free", space = "free_x") +
    geom_col(position = "dodge", color = "black", alpha = 0.8) +
    labs(
        title = NULL,
        x = NULL,
        y = "Number of edges",
        fill = "Correlation"
    )+
    scale_fill_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F8766D")) +
    
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )
p
name <- paste0(dir_name, "/module_correlation")
width <- 15
height <- 5.5
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

write.csv(all_data_df, paste0(name, "_data.csv"), quote = F, row.names = F)
write.csv(all_main_module, paste0(name, "_main_module.csv"), quote = F, row.names = F)
write.csv(all_main_edge, paste0(name, "_main_edge.csv"), quote = F, row.names = F)


# 2. Module Composition
top_order <- read.csv("../../01-sort_data/07-top_order/top_order.csv", header = T)
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
merge_p <- NULL
for (typ in type) {
    network <- ifelse(typ == "all", "Whole", ifelse(typ == "inter", "Interkingdom", "Intrakingdom"))
    main_module <- all_main_module %>%
        filter(Network == paste0(network, " network"))
    
    tmp_df <- merge(main_module, core_taxonomy, by.x = 'ID', by.y = "ASVID")
    
    p_list <- NULL
    for (kin in kingdom) {
        raw_df <- tmp_df %>%
            filter(Clade == kin) %>%
            mutate(Order = case_when(Order %in% top_order[, kin] ~ Order, TRUE ~ "Others"))
        
        all_df <- NULL
        for (mod in unique(raw_df$group)) {
            res_df <- raw_df %>%
                filter(group == mod)
            fin_df <- data.frame(t(table(res_df$Order))) %>%
                select(-Var1) %>%
                rename(Taxonomy = Var2, Count = Freq) %>%
                mutate(Module = paste0(mod, " (", sum(Count), ")"))
            
            all_df <- bind_rows(all_df, fin_df)
        }
        
        
        color_df <- data.frame(
            Taxonomy = c(top_order[, kin], 'Others'),
            Color = colorRampPalette(brewer.pal(9, 'Set1'))(num+1)
        )
        color_df <- color_df[color_df$Taxonomy %in% all_df$Taxonomy,]
        all_df$Taxonomy <- factor(all_df$Taxonomy, levels = color_df$Taxonomy)
        color_manual <- color_df$Color
        
        p <- ggplot(all_df, aes(x = Module, y = Count, fill = Taxonomy)) +
            geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
            labs(
                title = NULL,
                subtitle = NULL,
                x = NULL,
                y = paste0("Proportion of ASVs (", network, " network)")
            ) + 
            scale_y_continuous(labels = scales::percent) +
            scale_fill_manual(values = color_manual, name = ifelse(kin == "Protist", "Protists", kin)) +
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
        p_list[[kin]] <- p
        
        name <- paste0(dir_name, "/", typ, "_", kin, "_module_feature_distribution")
        write.csv(all_df, paste0(name, "_data.csv"), quote = F, row.names = F)
    }
    
    rel_widths <- switch(
        typ,
        "all"   = c(1, 0.95, 1.15),
        "inter" = c(1, 0.95, 1.15),
        c(1, 0.8, 0.98)
    )
    
    p <- cowplot::plot_grid(
        p_list[["Bacteria"]],
        p_list[["Fungi"]], 
        p_list[["Protist"]], 
        ncol = 3, nrow = 1, rel_widths = rel_widths
    )
    merge_p[[typ]] <- p
    name <- paste0(dir_name, "/", typ, "_module_feature_distribution")
    width <- 17
    height <- 5.5
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}
p <- merge_p[["all"]] / 
    merge_p[["inter"]] / 
    merge_p[["intra"]]
name <- paste0(dir_name, "/module_feature_distribution")
width <- 17
height <- 16.5
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------
all <- all_main_edge %>%
    filter(Network == "Whole network") %>%
    group_by(Module_OTU_1) %>%
    summarise(across(c(X1, Y1), mean, na.rm = TRUE))
all
inter <- all_main_edge %>%
    filter(Network == "Interkingdom network") %>%
    group_by(Module_OTU_1) %>%
    summarise(across(c(X1, Y1), mean, na.rm = TRUE))
inter
intra <- all_main_edge %>%
    filter(Network == "Intrakingdom network") %>%
    group_by(Module_OTU_1) %>%
    summarise(across(c(X1, Y1), mean, na.rm = TRUE))
intra
