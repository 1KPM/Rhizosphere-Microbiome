### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "12-each_order_hub_enrich"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(clusterProfiler)
library(patchwork)
library(ggh4x)
library(vegan)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

legend_position <- "top"
# ------------------------------------------------------------------------------

### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
kegg_name <- read.csv("00-data/K_gene_name.csv",row.names = 1)

pathway2ko <- read.csv('00-data/pathway2ko.csv')
pathway2superpathway <- read.csv('00-data/pathway2superpathway.csv')

# ------------------------------------------------------------------------------



##各个目hub富集—————————————————————————————————————————————————————————————————
tmp_list <- NULL
for (ord in names(top10order_list)) {
  hub_info <- read.csv(paste0('./03-get_top10order_network_property/', ord, '_inter_network_hub_info.csv'))
  hub_info <- hub_info %>% filter(roles != "Peripherals") %>% select(name)
  
  Bacteria_list <- gsub("^b", "", hub_info$name[startsWith(hub_info$name, "b")])
  Fungi_list <- gsub("^f", "", hub_info$name[startsWith(hub_info$name, "f")])
  Protist_list <-  gsub("^p", "", hub_info$name[startsWith(hub_info$name, "p")])
  max_len <- max(length(Bacteria_list), length(Fungi_list), length(Protist_list))
  
  df <- data.frame(
    Bacteria = c(Bacteria_list, rep(NA, max_len - length(Bacteria_list))),
    Fungi = c(Fungi_list, rep(NA, max_len - length(Fungi_list))),
    Protist = c(Protist_list, rep(NA, max_len - length(Protist_list)))
  )
  
  ko_res <- compareCluster(df, fun = 'enricher', TERM2GENE = pathway2ko)
  ko_res@compareClusterResult <- merge(ko_res@compareClusterResult, pathway2superpathway, by.x = 'ID', by.y = 'Pathway', all.x = T)
  ko_res@compareClusterResult$Order <- ord
  tmp_list[[ord]] <- ko_res@compareClusterResult
}
result <- do.call(rbind, tmp_list)
write.csv(result,paste0(dir_name,'/each_order_hub_KO_enrichment_results.csv'), row.names = F)
data_ko_df <- read.csv(paste0(dir_name,'/each_order_hub_KO_enrichment_results.csv'))

data_ko_df <- data_ko_df %>%
  filter(ID != "Carbon fixation pathways in prokaryotes") %>%
  mutate(
    GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )
data_ko_df$Order <- factor(data_ko_df$Order, levels = names(top10order_list))

p <- ggplot(data_ko_df, aes(x = Order, y = Description)) +
  geom_point(
    aes(size = GeneRatio * 1.4),
    position = position_nudge(x = 0.1, y = -0.1),
    color = "gray20",
    alpha = 0.4,
    shape = 19,
    show.legend = FALSE
  ) +
  geom_point(
    aes(size = GeneRatio, fill = p.adjust),
    color = "black",
    shape = 21,
    stroke = 0.8,
    alpha = 1 
  )+
  facet_wrap(~Cluster, nrow = 1, scales = "fixed",drop = FALSE) +
  scale_x_discrete(drop = FALSE) +  
  scale_fill_gradient(low = "#e17674", high = "#3b7eb8", name = "p.adjust") +
  scale_size_continuous(range = c(1, 5), name = "KO Ratio",breaks = c(0.1,0.2,0.3,0.4,0.5)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
    plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
    axis.title = element_text(size = 7, color = 'black'),
    axis.text = element_text(size = 6, color = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(size = 7, color = 'black', face = "bold"),
    strip.background = element_rect(fill = "gray80", color = "gray80"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    legend.title = element_text(size = 7, color = 'black'),
    legend.text = element_text(size = 6, color = 'black'),
    legend.box.margin = margin(0, 0, 0, -30),
    legend.key.size = unit(0.25, 'cm'),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    legend.position = legend_position
  ) +
  guides(
    fill = guide_colorbar(label.theme = element_text(angle = -90, size = 6, color = 'black', vjust = 0.5)),
    size = guide_legend(nrow = 2)
  ) +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL)

width = 18
height = 15
name <- paste0(dir_name, '/each_order_hub_ko')
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
#———————————————————————————————————————————————————————————————————————————————