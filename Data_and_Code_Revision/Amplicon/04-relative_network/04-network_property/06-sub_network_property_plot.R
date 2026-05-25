# ******************************************************************************
# @File: 06-sub_network_property_plot.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-18 10:58:46
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "06-sub_network_property_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(ggh4x)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
type <- c("all", "inter", "intra")

network <- c("Whole network", "Interkingdom network", "Intrakingdom Network")
main_property <- c(
    "Average degree", "Connectance", "Degree centralization", "Betweenness centralization", 
    "Clustering coefficient", "Robustness"
) 
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
rs_metadata <- read.csv("../../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv")
top_tree <- read.csv("../../01-sort_data/01-tree_color/Tree_top10_order_color.csv", header = T) 
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
all_df <- NULL
for (typ in type) {
    tmp_df <- read.csv(paste0("05-sub_network/", typ, "_network_property.csv"), header = T)
    tmp_df$Type <- ifelse(typ == "all", "Whole network", 
                          ifelse(typ == "inter", "Interkingdom network", "Intrakingdom Network"))
    all_df <- bind_rows(all_df, tmp_df)
}

fin_df <- all_df %>%
    left_join(rs_metadata[c("FileID", "Order")], by = c("Sample" = "FileID")) %>%
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
        Order = case_when(Order %in% top_tree$Order ~ Order, TRUE ~ "Others"),
        Order = factor(Order, levels = top_tree$Order)
    )

property <- unique(fin_df$Property)[! unique(fin_df$Property) %in% c("Edge connectivity", "Closeness centralization")]

all_data_df <- all_summarise_df<- all_pvalue_df <- data.frame()
for (net in network) {
    for (pro in property) {
        data_df <- fin_df %>%
            filter(Type == net, Property == pro) %>%
            rename( "Group" = Order)
        
        summarise_df <- data_df %>%
            group_by(Group) %>%
            summarise(
                n = n(),
                mean = mean(Value, na.rm = TRUE),
                min = min(Value, na.rm = TRUE),
                max = max(Value, na.rm = TRUE),
                sd = sd(Value, na.rm = TRUE),
                se = sd / sqrt(n),
                ci = qt(0.975, df = n - 1) * se
            ) %>%
            mutate(allmax = max(max))
   
        # Non-parametric test (Kruskal-Wallis rank-sum test with Dunn's post hoc test)
        kruskal_test <- kruskal.test(Value ~ Group, data = data_df)
        sig_value <- kruskal_test$p.value
        dunn_test <- rstatix::dunn_test(data_df, Value ~ Group, p.adjust.method = p_adjust_method)
        pvalue_df <- data.frame(KruskalWallis = sig_value, dunn_test[c("group1", "group2", "p.adj")])
        
        pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                                  ifelse(pvalue_df$p.adj < 0.01, "**", 
                                         ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
        n <- nrow(summarise_df)
        pvalue_matrix <- matrix(1, ncol = n, nrow = n)
        k <- 0
        for(i in 1:(n - 1)) { 
            for(j in (i + 1):n) { 
                k <- k + 1
                pvalue_matrix[i, j] <- pvalue_df$p.adj[k]
                pvalue_matrix[j, i] <- pvalue_df$p.adj[k]
            }
        }
        letter_df <- agricolae::orderPvalue(summarise_df$Group, summarise_df$mean, 0.05, pvalue_matrix, console = TRUE)
        letter_df <- letter_df[levels(data_df$Group),]
        summarise_df$label <- letter_df$groups
        
        data_df$Type <- net
        data_df$Property <- pro
        
        summarise_df$Type <- net
        summarise_df$Property <- pro
        
        pvalue_df$Type <- net
        pvalue_df$Property <- pro
        
        all_data_df <- bind_rows(all_data_df, data_df)
        all_summarise_df <- bind_rows(all_summarise_df, summarise_df)
        all_pvalue_df <- bind_rows(all_pvalue_df, pvalue_df)
    }
}

name <- paste0(dir_name, "/network_property")
write.csv(all_data_df, paste0(name, "_data.csv"), row.names = F)
write.csv(all_summarise_df, paste0(name, "_summarise.csv"), row.names = F)
write.csv(all_pvalue_df, paste0(name, "_pvalue.csv"), row.names = F)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. 主要网络特征
tmp_data_df <- all_data_df %>%
    filter(Property %in% main_property) %>%
    mutate(Type = factor(Type, levels = network))

tmp_summarise_df <- all_summarise_df %>%
    filter(Property %in% main_property) %>%
    mutate(Type = factor(Type, levels = network))

p <- ggplot(tmp_data_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +  
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(
        data = tmp_summarise_df,
        mapping = aes(x = Group, y = max + allmax * 0.2, label = label),
        position = position_dodge(0.9),
        size = 7 / 2.835
    ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) + 
    facet_grid2(
        rows = vars(Property), 
        cols = vars(Type), 
        scales = "free", 
        independent = "y",  # 这是解除同一行共享 Y 轴的关键！
        space = "free_x"
    ) + 
    theme_bw() + 
    scale_color_manual(values = top_tree$Color) + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(color = "black", size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
          panel.spacing = unit(0.1, "cm"),
          legend.position = 'none')

width <- 17.5
height <- 22
name <- paste0(dir_name, "/network_property_main")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


# 2. 其他网络特征
tmp_data_df <- all_data_df %>%
    filter(! Property %in% main_property) %>%
    mutate(Type = factor(Type, levels = network))

tmp_summarise_df <- all_summarise_df %>%
    filter(! Property %in% main_property) %>%
    mutate(Type = factor(Type, levels = network))

p <- ggplot(tmp_data_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +  
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(
        data = tmp_summarise_df,
        mapping = aes(x = Group, y = max + allmax * 0.2, label = label),
        position = position_dodge(0.9),
        size = 7 / 2.835
    ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = NULL
    ) + 
    facet_grid2(
        rows = vars(Property), 
        cols = vars(Type), 
        scales = "free", 
        independent = "y",  # 这是解除同一行共享 Y 轴的关键！
        space = "free_x"
    ) + 
    theme_bw() + 
    scale_color_manual(values = top_tree$Color) + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(color = "black", size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
          panel.spacing = unit(0.1, "cm"),
          legend.position = 'none')

width <- 17.5
height <- 22
name <- paste0(dir_name, "/network_property_others")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------