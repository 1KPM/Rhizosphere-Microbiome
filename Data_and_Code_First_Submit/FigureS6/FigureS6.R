pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(rstatix)
library(agricolae)

df <- read.csv("./data/FigureS6-MAG_abundace.csv",row.names = 1,check.names = F)
tax <- read.csv("./data/FigureS6-MAG_taxonomy.csv",row.names = 1)
metadata <- read.csv("./data/metadata.csv")
order_list1 <- c("Fabales","Rosales","Lamiales","Malpighiales","Sapindales","Gentianales","Asparagales","Malvales","Myrtales","Arecales")
metadata <- metadata %>%
  filter(TreeID %in% colnames(df)) %>%
  mutate(Group = ifelse(Order %in% order_list1, Order, "Others"))

df_long <- df %>%
  rownames_to_column("MAG") %>%
  pivot_longer(-MAG, names_to = "TreeID", values_to = "Abundance") %>%
  inner_join(metadata, by = "TreeID")

order_df <- df_long %>%
  group_by(MAG, Group) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = 'drop')

order_abundance <- order_df %>%
  pivot_wider(names_from = Group, values_from = TotalAbundance, values_fill = 0) %>%
  column_to_rownames("MAG")

tax_levels <- c("Phylum", "Class", "Order")
plot_list1 <- list()
for (tax_level in tax_levels) {
  tax_data <- tax %>%
    rownames_to_column("MAG") %>%
    select(MAG, !!tax_level)
  
  classified_abundance <- order_abundance %>%
    rownames_to_column("MAG") %>%
    left_join(tax_data, by = "MAG") %>% 
    select(-MAG) %>%
    rename(Classification = !!tax_level) %>%
    group_by(Classification) %>%
    summarise(across(everything(), sum)) %>% 
    column_to_rownames("Classification")
  
  classification_total <- classified_abundance %>% 
    rownames_to_column("Classification") %>%
    mutate(Total = rowSums(across(-Classification))) 
  
  top_classification <- classification_total %>% 
    filter(Classification != "Unassign") %>%
    slice_max(Total, n = 10) %>%
    pull(Classification) %>% 
    as.character()
  
  classification_processed <- classification_total %>% 
    mutate(
      Classification = case_when(
        Classification == "Unassign" ~ "Unassign",
        Classification %in% top_classification ~ Classification,
        TRUE ~ "Others"
      )
    ) %>% 
    group_by(Classification) %>% 
    summarise(across(everything(), sum))
  
  classification_order <- c(
    classification_total %>%
      filter(Classification %in% top_classification) %>% 
      arrange(desc(Total)) %>% 
      pull(Classification) %>% 
      as.character(),
    "Others",
    "Unassign"
  )
  
  classification_long <- classification_processed %>% 
    select(-Total) %>% 
    pivot_longer(-Classification, names_to = "Group", values_to = "Abundance") %>% 
    mutate(
      Group = factor(Group, levels = c(order_list1, "Others")),
      Classification = factor(Classification, levels = classification_order)
    )
  
  n_colors <- length(classification_order) - 1
  #fill_colors <- c(colorRampPalette(brewer.pal(11, "Paired"))(n_colors),"#BCBEC0")
  fill_colors <- c("#E41A1C","#596A98","#449B75","#6B886D","#AC5782","#FF7F00","#FFE528","#C9992C","#C66764","#E485B7","#999999","#666666")
  p <- ggplot(classification_long, aes(x = Group, y = Abundance, fill = Classification)) +
    geom_bar(stat = "identity", position = "fill", width = 0.85, color = NA) +
    scale_fill_manual(values = fill_colors, name = "Taxonomy") +
    scale_y_continuous(
      expand = c(0, 0),
      labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.2)
    ) +
    labs(
      x = "", 
      y = paste0("Relative abundance (", tax_level, ")")
    ) +
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
  
  plot_list1[[tax_level]] <- p
}

combined_plot1 <- wrap_plots(plot_list1, ncol = 1)



data <- read.csv("./data/FigureS6-MAG_info.csv")
data2 <- read.csv("./data/FigureS6-group.csv")
data <- left_join(data, data2 %>% select(MAG, size), by = "MAG")


top_Phylum <-c("Proteobacteria","Actinobacteriota","Acidobacteriota","Verrucomicrobiota","Chloroflexota","Bacteroidota","Patescibacteria","Desulfobacterota_B","Firmicutes","Planctomycetota")
top_Class <-c("Gammaproteobacteria","Alphaproteobacteria","Actinomycetia","Thermoleophilia","Acidobacteriae","Verrucomicrobiae","Vicinamibacteria","Acidimicrobiia","Bacteroidia","Saccharimonadia")
top_Order <- c("Rhizobiales","Pseudomonadales","Burkholderiales","Enterobacterales","Sphingomonadales","Solirubrobacterales","Gaiellales","Steroidobacterales","Chthoniobacterales","Vicinamibacterales")
top_Family <- c("Xanthobacteraceae","Moraxellaceae","Pseudomonadaceae","Enterobacteriaceae","Sphingomonadaceae","Burkholderiaceae","Steroidobacteraceae","Methyloligellaceae","Gaiellaceae","Rhizobiaceae")
top_Genus <-c("Acinetobacter","Pseudomonas_E","Methyloceanibacter","Pseudolabrys","Enterobacter","Bradyrhizobium","VAZQ01","Mycobacterium","Udaeobacter","Rhizobium")
plot_list2 <- list()

tax_levels <- c("Phylum","Class","Order")

for (tax in tax_levels) {
  top_tax <- get(paste0("top_",tax))
  fin_df <- data %>% 
    select(Group = !!sym(tax), Value = size) %>%
    filter(Group %in% top_tax) %>%
    mutate(Group = factor(Group, levels = top_tax))
  
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
  letter_df <- letter_df[levels(fin_df$Group), ] 
  letter_vector <- letter_df$groups
  color_manual1 <-  c("#E41A1C","#596A98","#449B75","#6B886D","#AC5782","#FF7F00","#FFE528","#C9992C","#C66764","#E485B7")
  group_levels1 <- levels(fin_df$Group)
  color_mapping1 <- setNames(color_manual1, group_levels1)
  
  
  p <- ggplot(fin_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = max_df, aes(y = Value * 1.2, label = letter_vector), position = position_dodge(0.9), size = 2.5) + 
    scale_color_manual(values = color_manual1) +
    labs(
      title = '',
      x = '',
      y = paste0(tax, ' genome size (MB)'),
    ) + 
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
          legend.position = "none"
      )
  
  plot_list2[[tax]] <- p
}

p <- 
  (plot_list1[[1]] | plot_list2[[1]])+
  plot_annotation(tag_levels = "a")

width <- 17.5
height <- 8
name <- "FigureS6"
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")

