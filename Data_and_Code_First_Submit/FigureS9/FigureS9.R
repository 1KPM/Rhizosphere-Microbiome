### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(1994)

# Import Packages
library(dplyr)
library(rstatix)
library(agricolae)
library(ggplot2)
library(RColorBrewer)

# Create Directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
type <- c("Closeness", "Betweenness")
network <- c("Within_Microbiome", "Cross_Microbiome", "Cross_KEGG")
# ------------------------------------------------------------------------------



### Get Results ----------------------------------------------------------------
# 1 Betweenness
typ <- "Betweenness"
all_data_df <- all_max_df <- all_sig_df <- data.frame()
for (net in network) {
    raw_df <- read.csv(paste0("data/", net, "_network_hub_info.csv"), row.names = 1)
    names(raw_df)[2:3] <- type
    
    res_df <- raw_df[typ]
    res_df$Clade <- ifelse(grepl("^b.*", row.names(res_df)), "Bacteria", 
                           ifelse(grepl("^f.*", row.names(res_df)), "Fungi", "Protist"))
    
    fin_df <- res_df
    names(fin_df) <- c('Value', 'Group')
    fin_df$Group <- factor(fin_df$Group, levels = c('Bacteria', 'Fungi', 'Protist'))
    
    max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
    max_df$GroupMax <- max(max_df$Value)
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
    letter_df <- letter_df[levels(fin_df$Group),]
    max_df$Label <- letter_df$groups
    
    fin_df$Network <- sub("_", "-Kingdom (", net) |> sub("$", ")", x = _)
    max_df$Network <- sub("_", "-Kingdom (", net) |> sub("$", ")", x = _)
    dunn_res_df$Network <- sub("_", "-Kingdom (", net) |> sub("$", ")", x = _)
    
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, max_df)
    all_sig_df <- rbind(all_sig_df, dunn_res_df)
}
all_data_df$Network <- factor(
    all_data_df$Network, levels = c("Within-Kingdom (Microbiome)", "Cross-Kingdom (Microbiome)", "Cross-Kingdom (KEGG)"))
all_max_df$Network <- factor(
    all_max_df$Network, levels = c("Within-Kingdom (Microbiome)", "Cross-Kingdom (Microbiome)", "Cross-Kingdom (KEGG)"))

p1 <- ggplot(all_data_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = all_max_df, aes(y = Value + GroupMax * 0.2, label = Label), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(
        x = '',
        y = paste0(typ, ' centrality of networks'),
    ) + 
    facet_grid(
        Network ~ .,
        scales = "free_y",
        switch = NULL,
        space = "fixed",
    ) +
    scale_color_manual(values = color_manual) + 
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
        legend.position = "none"
    )



name <- paste0(dir_name, '/FigureS9ABC')

write.csv(all_sig_df, paste0(name, '_sig_fino.csv'), row.names = F)   
write.csv(all_data_df, paste0(name, '_data.csv'), row.names = F)
write.csv(all_max_df, paste0(name, '_data_mean.csv'), row.names = F)
# ------------------------------------------------------------------------------

# 2 Closeness
typ <- "Closeness"
all_data_df <- all_max_df <- all_sig_df <- data.frame()
for (net in network) {
    raw_df <- read.csv(paste0("data/", net, "_network_hub_info.csv"), row.names = 1)
    names(raw_df)[2:3] <- type
    
    res_df <- raw_df[typ]
    res_df$Clade <- ifelse(grepl("^b.*", row.names(res_df)), "Bacteria", 
                           ifelse(grepl("^f.*", row.names(res_df)), "Fungi", "Protist"))
    
    fin_df <- res_df
    names(fin_df) <- c('Value', 'Group')
    fin_df$Group <- factor(fin_df$Group, levels = c('Bacteria', 'Fungi', 'Protist'))
    
    max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
    max_df$GroupMax <- max(max_df$Value)
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
    letter_df <- letter_df[levels(fin_df$Group),]
    max_df$Label <- letter_df$groups
    
    fin_df$Network <- sub("_", "-Kingdom (", net) |> sub("$", ")", x = _)
    max_df$Network <- sub("_", "-Kingdom (", net) |> sub("$", ")", x = _)
    dunn_res_df$Network <- sub("_", "-Kingdom (", net) |> sub("$", ")", x = _)
    
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, max_df)
    all_sig_df <- rbind(all_sig_df, dunn_res_df)
}
all_data_df$Network <- factor(
    all_data_df$Network, levels = c("Within-Kingdom (Microbiome)", "Cross-Kingdom (Microbiome)", "Cross-Kingdom (KEGG)"))
all_max_df$Network <- factor(
    all_max_df$Network, levels = c("Within-Kingdom (Microbiome)", "Cross-Kingdom (Microbiome)", "Cross-Kingdom (KEGG)"))

p2 <- ggplot(all_data_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = all_max_df, aes(y = Value + GroupMax * 0.2, label = Label), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(
        x = '',
        y = paste0(typ, ' centrality of networks'),
    ) + 
    facet_grid(
        Network ~ .,
        scales = "free_y",
        switch = NULL,
        space = "fixed",
    ) +

    scale_color_manual(values = color_manual) + 
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
        legend.position = "none"
    )

name <- paste0(dir_name, '/FigureS9DEF')

write.csv(all_sig_df, paste0(name, '_sig_fino.csv'), row.names = F)   
write.csv(all_data_df, paste0(name, '_data.csv'), row.names = F)
write.csv(all_max_df, paste0(name, '_data_mean.csv'), row.names = F)

p <- cowplot::plot_grid(
    p1 + theme(plot.margin = margin(0, 0.5, 0, 0.5, "cm")),
    p2, 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    labels = c("a", "b"), label_size = 10, label_fontface = "bold",
    ncol = 2, nrow = 1
)

width <- 14
height <- 12
name <- paste0(dir_name, '/FigureS9')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
