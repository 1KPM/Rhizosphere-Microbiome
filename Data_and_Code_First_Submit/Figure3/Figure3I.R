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
library(rstatix)
library(agricolae)


### Import data ----------------------------------------------------------------
hub_info <- read.csv('data/Meta_cross_network_hub_info.csv', row.names = 1)
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
raw_df <- hub_info
names(raw_df)[1:3] <- c('Degree', 'Closeness', 'Betweenness')
raw_df$Clade <- factor(raw_df$Clade, levels = c('Bacteria', 'Fungi', 'Protist'))
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
fin_df <- raw_df[c('Clade', "Degree")]
names(fin_df) <- c('Group', 'Value')

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
letter_df <- letter_df[levels(fin_df$Group),]
letter_vector <- letter_df$groups

p <- ggplot(fin_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = max_df, aes(y = Value * c(1.2, 1.1, 1.1), label = letter_vector), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(
        x = '',
        y = 'Degree centrality of cross-kingdom networks',
    ) + 
    theme_bw() + 
    coord_flip() +
    scale_color_manual(values = color_manual) + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')
p
width <- 7
height <- 3.5
name <- paste0(dir_name, '/Figure3I')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
ggsave(paste0(name, ".tiff"), p, width = width, height = height, dpi = 600, units = "cm", compression = "lzw")
write.csv(dunn_res_df, paste0(name, '.csv'), row.names = F)
# ------------------------------------------------------------------------------
