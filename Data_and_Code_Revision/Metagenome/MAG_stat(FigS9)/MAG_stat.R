pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

library(tidyverse)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(rstatix)
library(agricolae)


##Figure A

mag <- read.csv(file = "./data/MAG_info.csv")
mag$Category <- with(mag, ifelse(Completeness > 90 & Contamination <= 5, "High-quality draft",
                                 ifelse(Completeness >= 50 & Contamination <= 10, "Medium-quality draft", "low-quality drafts")))

p1 <- ggplot(data = mag, aes(x = Completeness, y = Contamination, color = Category)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("low-quality drafts" = "#4DAF4A", "Medium-quality draft" = "#FF6347", "High-quality draft" = "dodgerblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(title = "", x = "Completeness", y = "Contamination")

mag_summary <- mag %>%
  mutate(ContaminationGroup = cut(Contamination, breaks = seq(0, 10, by = 0.5), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(ContaminationGroup, Category) %>%
  summarise(Count = n())
mag_summary$Category <- factor(mag_summary$Category, levels = c("High-quality draft", "Medium-quality draft", "low-quality drafts"))

p2 <- ggplot(data = mag_summary, aes(x = as.numeric(ContaminationGroup)/2, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("low-quality drafts" = "#4DAF4A", "Medium-quality draft" = "#FF6347", "High-quality draft" = "dodgerblue"))+
  scale_x_continuous(breaks = seq(0, 10, by = 0.5), labels = seq(0, 10, by = 0.5)) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "")

mag_summary <- mag %>%
  mutate(CompletenessGroup = cut(Completeness, breaks = seq(50, 100, by = 2), include.lowest = TRUE, labels = FALSE)) %>%
  group_by(CompletenessGroup, Category) %>%
  summarise(Count = n())
mag_summary$Category <- factor(mag_summary$Category, levels = c("High-quality draft", "Medium-quality draft", "low-quality drafts"))

p3 <- ggplot(data = mag_summary, aes(x = as.numeric(CompletenessGroup)/2, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("low-quality drafts" = "#4DAF4A", "Medium-quality draft" = "#FF6347", "High-quality draft" = "dodgerblue"),labels = c("low-quality drafts" = "low-quality drafts","Medium-quality draft" = "Medium-quality draft(7011 MAGs,79.68%)","High-quality draft" = "High-quality draft(1788 MAGs,20.32%)"))+
  scale_x_continuous(breaks = seq(50, 100, by = 2), labels = seq(50, 100, by = 2)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.direction = "vertical",
        legend.key.size = unit(0.5, "lines")) +
  labs(x = "", y = "",fill="")

p_a <- p3 + plot_spacer()+p1 + p2 + plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 3))

width <- 10
height <- 10
name <- paste0(dir_name,"/MAG_stat_A")
ggsave(paste0(name, ".pdf"), p_a, width = width, height = height, dpi = 600,units = "cm")
ggsave(paste0(name, ".jpg"), p_a, width = width, height = height, dpi = 600,units = "cm")
ggsave(paste0(name, ".tiff"), p_a, width = width, height = height, dpi = 600,units = "cm")

##Figure B
df <- read.csv("./data/MAG_abundace.csv",row.names = 1,check.names = F)
tax <- read.csv("./data/MAG_taxonomy.csv",row.names = 1)
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

tax_levels <- c("Phylum")
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
  fill_colors <- c("#E41A1C","#596A98","#449B75","#6B886D","#AC5782","#FF7F00","#FFE528","#C9992C","#C66764","#E485B7","#999999","#666666")
  p <- ggplot(classification_long, aes(x = Group, y = Abundance, fill = Classification)) +
    geom_bar(stat = "identity", position = "fill", width = 0.85, color = NA) +
    scale_fill_manual(values = fill_colors,
                      guide = guide_legend(
                        ncol = 1,
                        title.position = "top",
                        title.hjust = 0,
                        keyheight = unit(3, "mm"),
                        keywidth = unit(3, "mm"),
                      ),
                      name = "Phylum"
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.2)
    ) +
    labs(
      x = "", 
      y = "Relative abundance",
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.title = element_text(colour = "black", size = rel(1.1)),
      axis.text.x = element_text(
        colour = "black",
        angle = 45, 
        hjust = 1, 
        vjust = 1,
      ),
      legend.key.size = unit(0.3, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.text = element_text(size = 8),
      legend.box.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      plot.margin = unit(c(0, 5, 0, 2), "mm")
    )
  
  plot_list1[[tax_level]] <- p
}


width <- 15
height <- 10
name <- paste0(dir_name,"/MAG_stat_B")
ggsave(paste0(name, ".pdf"), plot_list1[[1]], width = width, height = height, dpi = 600,units = "cm")
ggsave(paste0(name, ".jpg"), plot_list1[[1]], width = width, height = height, dpi = 600,units = "cm")
ggsave(paste0(name, ".tiff"), plot_list1[[1]], width = width, height = height, dpi = 600,units = "cm")

##Figure C
data <- read.csv("./data/MAG_info.csv")
data2 <- read.csv("./data/group.csv")
data <- left_join(data, data2 %>% select(MAG, size), by = "MAG")

top_Phylum <-c("Proteobacteria","Actinobacteriota","Acidobacteriota","Verrucomicrobiota","Chloroflexota","Bacteroidota","Patescibacteria","Desulfobacterota_B","Firmicutes","Planctomycetota")
plot_list2 <- list()

tax_levels <- c("Phylum")

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
  color_manual1 <- colorRampPalette(brewer.pal(9, 'Set1'))(11)[1:10]
  group_levels1 <- levels(fin_df$Group)
  color_mapping1 <- setNames(color_manual1, group_levels1)
  
  
  p <- ggplot(fin_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = max_df, aes(y = Value * 1.2, label = letter_vector), position = position_dodge(0.9), size = 2.5) + 
    labs(
      title = '',
      x = '',
      y = 'Genomn size(MB)',
    ) + 
    scale_color_manual(values = color_mapping1) +
    theme_bw() + 
    theme(
      plot.title = element_blank(), 
      axis.title = element_text(colour = "black", size = rel(1.1)),
      axis.text.x = element_text(
        colour = "black",
        angle = 45, 
        hjust = 1, 
        vjust = 1,
      ),
      panel.grid = element_blank(),
      legend.position = 'none',
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  plot_list2[[tax]] <- p
}


width <- 10
height <- 8
name <- paste0(dir_name,"/MAG_stat_C")
ggsave(paste0(name, ".pdf"), plot_list2[[1]], width = width, height = height, dpi = 600,units = "cm")
ggsave(paste0(name, ".jpg"), plot_list2[[1]], width = width, height = height, dpi = 600,units = "cm")
ggsave(paste0(name, ".tiff"), plot_list2[[1]], width = width, height = height, dpi = 600,units = "cm")
