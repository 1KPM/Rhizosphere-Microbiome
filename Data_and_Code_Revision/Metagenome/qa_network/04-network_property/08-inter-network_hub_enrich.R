### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "08-inter_network_enrich"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(patchwork)
library(cowplot) 
library(clusterProfiler)

### Define variable -----------------------------------------------------------
kingdom <- c('Bacteria', 'Fungi', 'Protists')
type <- c("inter")

### Import data ----------------------------------------------------------------
kegg_info <- read.csv('00-data/kegg_pathway_for_1KPM.csv')
pathway_color <- read.csv("00-data/superpathway_color.csv")
pathway2ko <- read.csv('00-data/pathway2ko.csv')
pathway2superpathway <- read.csv('00-data/pathway2superpathway.csv')
# ------------------------------------------------------------------------------

plot_list <- list()
# 1.summary
hub_info <- read.csv('../02-inter_network/01-get_inter_network_property/inter_network_hub_info.csv', row.names = 1)
  
hub_info <- hub_info %>%
  mutate(Clade = case_when(
    substr(.[[1]], 1, 1) == "b" ~ "Bacteria",
    substr(.[[1]], 1, 1) == "f" ~ "Fungi",
    substr(.[[1]], 1, 1) == "p" ~ "Protist",
    TRUE ~ NA_character_
  )) %>%
  mutate(Clade = ifelse(Clade == "Protist", "Protists", Clade),
         Clade = factor(Clade, levels = c("Bacteria", "Fungi", "Protists")))
hub_info$KO <- gsub("^[bfp]", "", row.names(hub_info))
all_tmp_df <- merge(hub_info, kegg_info, by = 'KO',all.x = T)

all_data_df <- data.frame()
for (kin in kingdom) {
  raw_df <- all_tmp_df[all_tmp_df$Clade == kin,]
  all_df <- data.frame(t(table(raw_df$level2)))
  hub_df <- data.frame(t(table(raw_df[raw_df$roles != 'Peripherals', 'level2'])))
    
  res_df <- merge(all_df[2:3], hub_df[2:3], by = 'Var2', all.x = T)
  res_df[is.na(res_df)] <- 0
  names(res_df) <- c('Taxonomy', paste0('All (', sum(hub_info$Clade == kin), ')'), paste0('Hub (', sum(hub_info$Clade == kin & hub_info$roles != "Peripherals"), ')'))
  
  fin_df <- melt(res_df, id.vars = 'Taxonomy', variable.name = 'ASV', value.name = 'Value')
  fin_df$Kingdom <- kin
  
  all_data_df <- bind_rows(all_data_df, fin_df)
}
all_data_df$Taxonomy <- factor(all_data_df$Taxonomy, levels = pathway_color$superpathway)

p_inter <- ggplot(all_data_df, aes(x = ASV, y = Value, fill = Taxonomy)) +
  geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Proportion of KOs (Inter network)"
  ) + 
  facet_grid(cols = vars(Kingdom), scales = "free", space = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = pathway_color$fill_color2, name = "Superpathway") +
  theme_bw() + 
  theme(
    text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(color = "black", size = 10, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.title = element_text(size = 10, color = 'black'), 
    legend.text = element_text(size = 9,  color = 'black'), 
    legend.box.margin = margin(0, 0, 0, -20),
    legend.key.size = unit(0.4, 'cm'),
    legend.box.spacing = unit(1,"cm"),
    legend.spacing = unit(0.2, "cm"),
    legend.position = "right",
  ) + 
  guides(fill = guide_legend(ncol = 1))


# 2.enrich
hub_info <- read.csv('../02-inter_network/01-get_inter_network_property/inter_network_hub_info.csv', row.names = 1)

hub_info <- hub_info %>%
  mutate(Clade = case_when(
    substr(.[[1]], 1, 1) == "b" ~ "Bacteria",
    substr(.[[1]], 1, 1) == "f" ~ "Fungi",
    substr(.[[1]], 1, 1) == "p" ~ "Protist",
    TRUE ~ NA_character_
  )) %>%
  mutate(Clade = ifelse(Clade == "Protist", "Protists", Clade),
         Clade = factor(Clade, levels = c("Bacteria", "Fungi", "Protists")))
hub_df <- hub_info %>% 
  filter(roles != "Peripherals") %>%
  group_by(Clade) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = row_id,
    names_from = Clade,
    values_from = name
  ) %>%
  select(-row_id) %>%
  mutate(across(everything(), ~ if_else(!is.na(.), str_sub(., 2), NA_character_)))

ko_res <- compareCluster(hub_df, fun = 'enricher', TERM2GENE = pathway2ko)
ko_res@compareClusterResult <- merge(ko_res@compareClusterResult, pathway2superpathway, by.x = 'ID', by.y = 'Pathway', all.x = T)
write.csv(ko_res@compareClusterResult,paste0(dir_name,'/each_kingdom_hub_KO_enrichment_results.csv'), row.names = F)

data_ko_df <- read.csv(paste0(dir_name,'/each_kingdom_hub_KO_enrichment_results.csv'))

data_ko_df <- data_ko_df %>%
  filter(ID != "Carbon fixation pathways in prokaryotes") %>%
  mutate(
    GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )
levels <- c("Two-component system",
            "Glyoxylate and dicarboxylate metabolism",
            "Valine, leucine and isoleucine degradation",
            "Fatty acid degradation",
            "Propanoate metabolism",
            "Citrate cycle (TCA cycle)",
            "Flagellar assembly",
            "Glycolysis / Gluconeogenesis",
            "Oxidative phosphorylation",
            "Selenocompound metabolism",
            "Pyrimidine metabolism",
            "Pyruvate metabolism",
            "Tryptophan metabolism"
            )

data_ko_df$Description <- factor(data_ko_df$Description, levels =levels)
p_enrich <- ggplot(data_ko_df, aes(x = Cluster, y = Description)) +
  geom_point(
    aes(size = GeneRatio * 1.01),
    position = position_nudge(x = 0.04, y = -0.04),
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
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) + 
  scale_fill_gradient(low = "#e17674", high = "#3b7eb8", name = "p.adjust") +
  scale_size_continuous(range = c(4, 10), name = "KO Ratio") +
  scale_y_discrete(limits = rev,position = "right") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 10, color = 'black'), 
        axis.text = element_text(size = 9, color = 'black'), 
        strip.text = element_text(size = 10, color = 'black'), 
        panel.spacing = unit(0.1, "cm"),
        legend.title = element_text(size = 10, color = 'black'), 
        legend.text = element_text(size = 9,  color = 'black'), 
        legend.box.margin = margin(0, 0, 0, -20),
        legend.key.size = unit(0.25, 'cm'),
        legend.box.spacing = unit(1,"cm"),
        legend.position = "left",
        plot.margin = margin(0, 0, 0, 20, "pt")) +
  guides(
    fill = guide_colorbar(
      label.theme = element_text(angle = 0, size = 8, color = 'black', vjust = 0.5)
    ),
    size = guide_legend(ncol = 1) 
  )

combined_plot <- plot_grid(p_inter, p_enrich, ncol = 1)

width = 17
height = 21
name <- paste0(dir_name, '/inter_kingdom_enrich')
ggsave(paste0(name, ".png"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), combined_plot, width = width, height = height, units = "cm")
