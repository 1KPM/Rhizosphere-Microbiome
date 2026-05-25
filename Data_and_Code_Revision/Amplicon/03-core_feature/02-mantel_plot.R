### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "02-mantel_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(vegan)
library(ggExtra)
library(scales)
library(ggtext)
library(patchwork)

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
occurrence_cutoff <- seq(0.1, 0.9, 0.1)
kingdom <- c("16S", "ITS", "Protist")

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (kin in kingdom) {
    otu_path <- paste0("../00-rawdata/feature_table/", kin, "/fin_feature_table_relative.tsv")
    otu <- read.delim(otu_path, row.names = 1)
    
    otu_PA <- 1 * ((otu > 0) == 1)
    otu_occ <- rowSums(otu_PA) / ncol(otu_PA)
    
    data_df <- data.frame(Proportion = seq(0.05, 1, 0.05), Numbers = NA)
    for (i in seq(0.05, 1, 0.05)) {
        data_df[data_df$Proportion == i, "Numbers"] <- length(otu_occ[otu_occ >= i])
    }
    
    name <- paste0(dir_name, "/", kin, "_sample_proportion")
    write.csv(data_df, paste0(name, ".csv"), quote = F, row.names = F)
    
    y_label <- ifelse(kin == "16S", "Bacteria", ifelse(kin == "ITS", "Fungi", "Protists"))
    p <- ggplot(data = data_df, mapping = aes(x = Proportion, y = Numbers)) +
        geom_point(
            size = 1, color = "#4c92c3"
        ) + 
        geom_vline(aes(xintercept = 0.2), linetype = "dashed", linewidth = 0.5, color = "#666")+
        
        labs(
            title = NULL,
            subtitle = NULL,
            x = "Proportion of samples",
            y = paste0("Number of features (", y_label, ")")
        ) + 
        theme_bw() + theme(
            text = element_text(color = "black", size = 6),
            plot.title = element_text(size = 7, hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.title = element_text(size = 7),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6, color = "black"),
            panel.spacing = unit(0.1, "cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.key.size = unit(0.25, "cm"),
            legend.position = "right"
        )
    
    assign(paste0("p_", kin), p)
}

p <- cowplot::plot_grid(
    get("p_16S") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
    get("p_ITS") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")), 
    get("p_Protist") + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 3, nrow = 1, rel_widths = c(1, 1, 1) 
)

width <- 17
height <- 5
name <- paste0(dir_name, "/sample_proportion")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


cutoff <- 0.2
plot_list <- list() 

mantel_df <- read.csv("01-mantel_results/core_feature_mantel_results.csv")
for (kin in kingdom) {
    bc_full <- read.csv(paste0("01-mantel_results/bray_curtis_", kin, "_", cutoff, "_full.csv"), row.names = 1)
    bc_sub <- read.csv(paste0("01-mantel_results/bray_curtis_", kin, "_", cutoff, "_sub.csv"), row.names = 1)
    mantel <- subset(mantel_df, Kingdom == kin)
    
    data_df <- data.frame(Core = as.dist(bc_sub), All = as.dist(bc_full))
    
    p_base <- ggplot(data_df, aes(x = All, y = Core)) +
        geom_point(alpha = 0.4, shape = 18, size = 2) + 
        geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 1.5) +
        geom_abline(slope = 1,  color = "red", linewidth = 1.5) +
        xlim(0, 1) + ylim(0, 1) +
        labs(
            x = "Bray-Curtis distance (All ASV)", 
            y = "Bray-Curtis distance (Core ASV)"
        ) +
        theme_bw() + theme(
            text = element_text(color = "black", size = 6),
            legend.title = element_text(size = 7),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6, color = "black")
        )
    
    p_marg <- ggMarginal(
        p_base,
        type = "histogram",
        fill = "orange",    
        xparams = list(fill = "mediumseagreen"), 
        bins = 30            
    )
    
    p_final <- wrap_elements(p_marg) + 
        labs(
            title = paste0("Mantel's: r = ", round(mantel$mantel_r, 2), ",  p = ", mantel$mantel_p),
            subtitle = paste0("lm(y ~ x): r<sup>2</sup> = ", round(mantel$linear_r2, 2))
        ) +
        theme(
            plot.title = element_text(size = 7, hjust = 0.5, color = "black"),
            plot.subtitle = element_markdown(size = 7, hjust = 0.5, color = "black")
        )
    
    plot_list[[kin]] <- p_final
}

p_combined <- plot_list[["16S"]] | plot_list[["ITS"]] | plot_list[["Protist"]]

width <- 17.5
height <- 6
name <- paste0(dir_name, "/mantel_plot")
ggsave(paste0(name, ".pdf"), p_combined, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------


