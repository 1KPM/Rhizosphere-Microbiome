pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(agricolae)
library(patchwork)


node_mt2 <- read.csv("./data/FigureS17-Gentianales_node_mt2.csv",row.names = 1)
edge_mt2 <- read.csv("./data/FigureS17-Gentianales_edge_mt2.csv")
cross_mt2_p <- ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Correlation), 
               alpha = 0.01, linewidth = 0.1, data = edge_mt2) + 
  geom_point(aes(X1, X2, fill = Clade, size = Degree^0.5), pch = 21, data = node_mt2) +
  scale_fill_manual(values = c("Bacteria"="#7FC97F","Fungi"="#BEAED4","Protist"="#FDC086"))  +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  coord_flip() +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = "cm"))

width <- 20
height <- 14
ggsave('./FigureS17D.jpg', cross_mt2_p, width = width, height = height, dpi = 900, bg = "white", units = 'cm')
ggsave('./FigureS17D.tiff', cross_mt2_p, width = width, height = height, units = 'cm',dpi = 900, bg = "white")
ggsave('./FigureS17D.pdf', cross_mt2_p, width = width, height = height, units = 'cm',dpi = 900, bg = "white")

df<-read.csv('./data/FigureS17-Gentianales_Cross_kingdom_network_hub_info.csv')
all_df <- df %>%
  mutate(
    Taxonomy = case_when(
      startsWith(X, "b") ~ "Bacteria",
      startsWith(X, "f") ~ "Fungi",
      startsWith(X, "p") ~ "Protist"
    )
  ) %>% 
  count(Taxonomy, name = "count") %>%
  mutate(ratio = sprintf("%.1f%%", 100 * count / sum(count)))
hub_df <- df %>%
  filter(roles != "Peripherals") %>%
  mutate(
    Taxonomy = case_when(
      startsWith(label, "b") ~ "Bacteria",
      startsWith(label, "f") ~ "Fungi",
      startsWith(label, "p") ~ "Protist"
    )
  ) %>% 
  count(Taxonomy, name = "count") %>%
  mutate(ratio = sprintf("%.1f%%", 100 * count / sum(count)))


p_all <- ggpubr::ggpie(
  data = all_df,
  x = "count",
  label ="ratio",
  fill = "Taxonomy",
  lab.pos = "in",
  lab.font = 3,
  palette = c("Bacteria" = "#7FC97F","Fungi" = "#BEAED4", "Protist" = "#FDC086")) +
  labs(title = "All KO(5258)") +
  theme(
    plot.title = element_text(size = 8,hjust = 0.5, face = "bold"),
    legend.position = "none",
    aspect.ratio = 1)

p_hub <- ggpubr::ggpie(
  data = hub_df,
  x = "count",
  label ="ratio",
  fill = "Taxonomy",
  lab.pos = "in",
  lab.font = 3,
  palette = c("Bacteria" = "#7FC97F","Fungi" = "#BEAED4", "Protist" = "#FDC086")) +
  labs(title = "Hub KO(130)") +
  theme(
    plot.title = element_text(size = 8,hjust = 0.5, face = "bold"),
    legend.position = "none",
    aspect.ratio = 1)


raw_df<-read.csv('./data/FigureS17-Gentianales_Cross_kingdom_network_hub_info.csv')
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
names(raw_df)[2:4] <- c('Degree','Closeness','Betweenness')

fin_df <- raw_df[c('Clade','Degree')]
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
fin_df$Group <- as.factor(fin_df$Group)
letter_df <- letter_df[levels(fin_df$Group),]
letter_vector <- letter_df$groups

p_degree <- ggplot(fin_df, aes(x = Group, y = Value)) + 
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
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text = element_text(size = 7, color = 'black'), 
        legend.title = element_text(size = 7, color = 'black'), 
        legend.text = element_text(size = 6,  color = 'black'), 
        axis.text.x = element_text(angle = 0, hjust = 0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        legend.position = 'none',
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"))



#final_plot1 <- cross_mt2_p + (p_all | p_hub) / p_degree
final_plot <- (p_all | p_hub) / p_degree
width <- 10
height <- 10
ggsave('./FigureS17EF.jpg', final_plot, width = width, height = height, units = 'cm',dpi = 900, bg = "white")
ggsave('./FigureS17EF.tiff', final_plot, width = width, height = height, units = 'cm',dpi = 900, bg = "white")
ggsave('./FigureS17EF.pdf', final_plot, width = width, height = height, units = 'cm',dpi = 900, bg = "white")

p4 <- cross_mt2_p + final_plot +plot_layout(widths = c(3,2))
width <- 30
height <- 15
ggsave('./FigureS17-Gentianales.pdf', p4, width = width, height = height, units = 'cm')
ggsave('./FigureS17-Gentianales.jpg', p4, width = width, height = height, units = 'cm')
ggsave('./FigureS17-Gentianales.tiff', p4, width = width, height = height, units = 'cm')