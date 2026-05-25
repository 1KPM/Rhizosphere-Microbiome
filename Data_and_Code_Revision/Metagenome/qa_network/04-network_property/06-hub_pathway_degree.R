### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "06-hub_pathway_degree"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(agricolae)
library(patchwork)

data <- read.csv("../02-inter_network/01-get_inter_network_property/inter_network_hub_info.csv",row.names = 1)
names(data)[1:4] <- c("KO", "Degree", "Betweenness", "Closeness") 
pathway2ko <- read.csv("00-data/pathway2ko.csv")
pathway2superpathway <- read.csv("00-data/pathway2superpathway.csv")
superpathway2ko <- pathway2ko %>%
  left_join(pathway2superpathway, by = "Pathway") %>%
  select(KO, Pathway, Superpathway)

superpathway_color <- read.csv("00-data/superpathway_color.csv")
superpathway_levels <- levels(factor(superpathway_color$superpathway, 
                                     levels = unique(superpathway_color$superpathway)))
superpathway_color_manual <-setNames(superpathway_color$fill_color2, superpathway_color$superpathway)


pathway_color <- read.csv("00-data/pathway_color.csv")
pathway_levels <- levels(factor(pathway_color$Pathway, 
                                     levels = unique(pathway_color$Pathway)))
pathway_color_manual <-setNames(pathway_color$fill_color, pathway_color$Pathway)

categories <- c(f = "Fungi", p = "Protist",b = "Bacteria")
pathway_result_list <- list()
pathway_plot_list <- list()
# Pathway Plot
for (prefix in names(categories)) {
  sub_data <- data %>%
    filter(grepl(paste0("^", prefix), KO)) %>%
    mutate(
      KO_clean = sub(paste0("^", prefix), "", KO),
      Category = categories[prefix]
    )
  
  merged_data <- left_join(sub_data, 
                           pathway2ko[, c("KO", "Pathway")], 
                           by = c("KO_clean" = "KO"))
  
  pathway_result_list[[categories[prefix]]] <- merged_data %>%
    select(Pathway, Degree, roles) %>%
    na.omit()
  
  fin_df <- pathway_result_list[[categories[prefix]]] %>%
    filter(Pathway %in% pathway_levels)
  
  names(fin_df) <- c('Group', 'Value' ,'Roles')
  max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
  mean_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = mean)
  
  n <- nrow(mean_df)
  dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
  dunn_res_df <- data.frame(dunn_res[-1])
  
  pvalue_df <- matrix(1, ncol = n, nrow = n)
  k <- 0
  for(i in 1:(n - 1)) { 
    for(j in (i + 1):n){ 
      k <- k + 1
      pvalue_df[i,j] <- dunn_res_df$p.adj[k]
      pvalue_df[j,i] <- dunn_res_df$p.adj[k]
    }
  }
  
  letter_df <- orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = TRUE)
  letter_df <- as.data.frame(letter_df)
  letter_df$Group <- rownames(letter_df)
  
  max_df <- left_join(max_df, letter_df[, c("Group", "groups")], by = "Group") %>%
    rename(label = groups)
  
  fin_df <- fin_df %>%
    mutate(
      point_color = ifelse(
        Roles == "Peripherals",
        "grey80",
        pathway_color_manual[as.character(Group)]
      )
    )
  fin_df$Group <- factor(fin_df$Group, levels = pathway_levels)
  max_df$Group <- factor(max_df$Group, levels = pathway_levels)
  p <- ggplot(fin_df, aes(x = Group, y = Value)) + 
    geom_point(
      aes(color = point_color, group = Group),
      position = position_jitterdodge(dodge.width = 0.6),
      alpha = 0.8, size = 2
    ) +
    scale_color_identity() + 
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE,drop = FALSE) +  
    geom_text(data = max_df, aes(y = Value * 1.2, label = label), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(title = if (prefix == "b") "Pathway" else "", x = '', y = paste0("Degree centrality of networks(",categories[prefix],")")) + 
    theme_bw() + 
    theme(
      plot.title = element_text(size = 10, color = 'black', hjust = 0.5), 
      axis.title.y = element_text(size = 6, color = 'black'), 
      axis.text = element_text(size = 6, color = 'black'), 
      panel.grid = element_blank(),
      legend.position = 'none',
      axis.title.x = element_blank(),
      axis.text.x = if (prefix == "p") element_text(size = 6,angle = 45,hjust = 1) else element_blank(),,
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  pathway_plot_list[[categories[prefix]]] <- p
}


# Superpathway Plot
superpathway_result_list <- list()
superpathway_plot_list <- list()
for (prefix in names(categories)) {
  sub_data <- data %>%
    filter(grepl(paste0("^", prefix), KO)) %>%
    mutate(
      KO_clean = sub(paste0("^", prefix), "", KO),
      Category = categories[prefix]
    )
  
  merged_data <- left_join(sub_data, 
                           superpathway2ko[, c("KO", "Superpathway")], 
                           by = c("KO_clean" = "KO"))
  
  superpathway_result_list[[categories[prefix]]] <- merged_data %>%
    select(Superpathway, Degree, roles) %>%
    na.omit()
  
  fin_df <- superpathway_result_list[[categories[prefix]]] %>%
    filter(Superpathway %in% superpathway_levels)
  
  names(fin_df) <- c('Group', 'Value' ,'Roles')
  max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
  mean_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = mean)
  
  n <- nrow(mean_df)
  dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
  dunn_res_df <- data.frame(dunn_res[-1])
  
  pvalue_df <- matrix(1, ncol = n, nrow = n)
  k <- 0
  for(i in 1:(n - 1)) { 
    for(j in (i + 1):n){ 
      k <- k + 1
      pvalue_df[i,j] <- dunn_res_df$p.adj[k]
      pvalue_df[j,i] <- dunn_res_df$p.adj[k]
    }
  }
  
  letter_df <- orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = TRUE)
  letter_df <- as.data.frame(letter_df)
  letter_df$Group <- rownames(letter_df)
  
  max_df <- left_join(max_df, letter_df[, c("Group", "groups")], by = "Group") %>%
    rename(label = groups)
  
  fin_df <- fin_df %>%
    mutate(
      point_color = ifelse(
        Roles == "Peripherals",
        "grey80",
        superpathway_color_manual[as.character(Group)]
      )
    )
  fin_df$Group <- factor(fin_df$Group, levels = superpathway_levels)
  max_df$Group <- factor(max_df$Group, levels = superpathway_levels)
  p <- ggplot(fin_df, aes(x = Group, y = Value)) + 
    geom_point(
      aes(color = point_color, group = Group),
      position = position_jitterdodge(dodge.width = 0.6),
      alpha = 0.8, size = 2
    ) +
    scale_color_identity() + 
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE,drop = FALSE) +  
    geom_text(data = max_df, aes(y = Value * 1.2, label = label), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(title = if (prefix == "b") "Superpathway" else "", x = '', y = '') + 
    theme_bw() + 
    theme(
      plot.title = element_text(size = 10, color = 'black', hjust = 0.5), 
      axis.title.y = element_text(size = 6, color = 'black'), 
      axis.text = element_text(size = 6, color = 'black'), 
      panel.grid = element_blank(),
      legend.position = 'none',
      axis.title.x = element_blank(),
      axis.text.x = if (prefix == "p") element_text(size = 6,angle = 45,hjust = 1) else element_blank(),,
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  superpathway_plot_list[[categories[prefix]]] <- p
}

final_plot <- wrap_plots(
  pathway_plot_list[['Bacteria']],superpathway_plot_list[['Bacteria']],
  pathway_plot_list[['Fungi']],superpathway_plot_list[['Fungi']],
  pathway_plot_list[['Protist']],superpathway_plot_list[['Protist']],
  ncol = 2) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(0, 3, 0, 3))


width <- 17.5
height <- 20.5
name <- paste0(dir_name,"/hub_pathway_degree")
ggsave(paste0(name, ".png"), final_plot, width = width, height = height, units = "cm",dpi = 900,bg = "white")
ggsave(paste0(name, ".pdf"), final_plot, width = width, height = height, units = "cm",dpi = 900,bg = "white")
ggsave(paste0(name, ".tiff"), final_plot, width = width, height = height, units = "cm",dpi = 900,bg = "white")
