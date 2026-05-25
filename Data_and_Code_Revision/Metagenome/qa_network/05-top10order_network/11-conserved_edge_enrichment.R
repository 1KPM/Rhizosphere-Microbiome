### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Create directory
dir_name <- "11-conserved_edge_enrichment"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Set seed
set.seed(1994)

# Import package
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(clusterProfiler)
library(patchwork)
library(ggh4x)
library(data.table)
library(UpSetR)
library(zoo)

### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
pathway2ko <- read.csv('00-data/pathway2ko.csv')
pathway2superpathway <- read.csv('00-data/pathway2superpathway.csv')
# ------------------------------------------------------------------------------
net <- "inter"

# Data Prepare -----------------------------------------------------------------
tmp_list <- NULL
for (ord in names(top10order_list)) {
  file_path <- paste0("03-get_top10order_network_property/", ord, "_", net, "_edge_property_raw.csv")
  edge_mt2 <- read.csv(file_path, header = T)
  edge <- edge_mt2 %>% select(node1,node2,label)
  tmp_list[[ord]] <- edge
}

all_pairs <- do.call(rbind, lapply(tmp_list, function(df) df[, c("node1", "node2")]))
all_pairs <- as.data.frame(t(apply(all_pairs, 1, sort)))
colnames(all_pairs) <- c("node1", "node2")
unique_edges <- unique(all_pairs)

edge_names <- names(tmp_list)

edge_maps <- map(tmp_list, function(df) {
  df %>%
    rowwise() %>%
    mutate(edge_key = paste(sort(c(node1, node2)), collapse = "_")) %>%
    ungroup() %>%
    select(edge_key, label) %>%
    deframe()
})

unique_edges <- unique_edges %>%
  rowwise() %>%
  mutate(edge_key = paste(sort(c(node1, node2)), collapse = "_")) %>%
  ungroup()

setDT(unique_edges)
unique_edges[, edge_key := paste(pmin(node1, node2), pmax(node1, node2), sep = "_")]

edge_names <- names(tmp_list)

for (nm in edge_names) {
  dt <- as.data.table(tmp_list[[nm]])[, .(node1, node2, label)]
  dt[, edge_key := paste(pmin(node1, node2), pmax(node1, node2), sep = "_")]
  dt[, label := fcase(label == "+", 1L,
                      label == "-", -1L,
                      default = 0L)]
  
  setkey(unique_edges, edge_key)
  setkey(dt, edge_key)
  unique_edges <- merge(unique_edges, dt[, .(edge_key, label)], 
                        by = "edge_key", all.x = TRUE, sort = FALSE)
  setnames(unique_edges, "label", nm)
  unique_edges[is.na(get(nm)), (nm) := 0L]
}

unique_edges[, edge_key := NULL]

write.csv(unique_edges,paste0(dir_name,'/unique_edges.csv'),row.names = F)
#-------------------------------------------------------------------------------
all_edge_result <- read.csv(paste0(dir_name,'/unique_edges.csv'))

all_order_edge <- all_edge_result %>%
  mutate(across(3:12, ~ case_when(
    . == -1 ~ "N",
    . == 0  ~ "U",
    . == 1  ~ "P"
  ))) %>%
  unite("label", all_of(3:12), sep = "", remove = FALSE) %>%
  select(node1,node2,label)

otu_grouped <- all_order_edge %>%
  group_by(label) %>%
  summarise(unique_otu = list(unique(c(node1, node2))), .groups = "drop") %>%
  deframe()

final_df <- as.data.frame(
  do.call(cbind, lapply(otu_grouped, `length<-`, max(lengths(otu_grouped))))
)

tmp_df <- final_df %>%
  select(
    PPPPPPPPPP,
    UPPPPPPPPP,
    PUPPPPPPPP,
    PPUPPPPPPP,
    PPPUPPPPPP,
    PPPPUPPPPP,
    PPPPPUPPPP,
    PPPPPPUPPP,
    PPPPPPPUPP,
    PPPPPPPPUP,
    PPPPPPPPPU
  )
prefix_map <- c("b" = "Bacteria", "f" = "Fungi", "p" = "Protists")


pathway2ko <- read.csv('00-data/pathway2ko.csv')
pathway2superpathway <- read.csv('00-data/pathway2superpathway.csv')
for (prefix in names(prefix_map)) {
  final_df_filter <- as.data.frame(lapply(tmp_df, function(x) ifelse(grepl(paste0("^", prefix), x, perl = TRUE), x, NA)))
  final_df_filter <- as.data.frame(lapply(final_df_filter, function(x) {
    sub(paste0("^", prefix), "", x, perl = TRUE)
  }))
  final_df_filter <- final_df_filter[rowSums(is.na(final_df_filter)) < ncol(final_df_filter), ]
  final_df_filter <- as.data.frame(
    lapply(final_df_filter, function(x) na.locf(x, fromLast = TRUE, na.rm = FALSE))
  )
  col_unique <- lapply(final_df_filter, unique)
  max_length <- max(sapply(col_unique, length))
  final_df_col_unique <- as.data.frame(
    lapply(col_unique, function(col) c(col, rep(NA, max_length - length(col))))
  )
  ko_res <- compareCluster(final_df_col_unique, fun = 'enricher', TERM2GENE = pathway2ko)
  ko_res@compareClusterResult <- merge(ko_res@compareClusterResult, pathway2superpathway, by.x = 'ID', by.y = 'Pathway', all.x = T)
  write.csv(ko_res@compareClusterResult, paste0(dir_name, '/', prefix_map[prefix], '_conversed_positive_edge_enrichment_results.csv'), row.names = F)
}

data_list <- list()
plot_list <- list()
kingdoms <- c("Bacteria", "Fungi", "Protists")
#kingdoms <- c("Bacteria", "Fungi")
all_clusters <- c("PPPPPPPPPP","UPPPPPPPPP","PUPPPPPPPP","PPUPPPPPPP","PPPUPPPPPP",
                  "PPPPUPPPPP","PPPPPUPPPP","PPPPPPUPPP","PPPPPPPUPP","PPPPPPPPUP",
                  "PPPPPPPPPU")
for (current_kingdom in kingdoms) {
  data_df <- read.csv(paste0(dir_name, '/', current_kingdom, '_conversed_positive_edge_enrichment_results.csv'))
  data_df <- data_df %>%
    mutate(
      GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
  #data_df <- data_df %>% filter(p.adjust < 0.01 & GeneRatio > 0.01 & FoldEnrichment > 2 )
  data_df$Kingdom <- current_kingdom
  
  data_list[[current_kingdom]] <- data_df
}


all_data <- bind_rows(data_list)
all_data <- all_data %>%
  group_by(Description) %>%
  filter(n_distinct(Kingdom) >= 2) %>%
  ungroup()
################################################################################
data_ko_df <- all_data %>% filter(Kingdom == "Bacteria")
data_ko_df <- setDT(data_ko_df)
cluster_mat <- dcast(data_ko_df, Cluster ~ Description, value.var = "Count", fill = 0)
rownames(cluster_mat) <- cluster_mat$Cluster
cluster_mat <- cluster_mat[, -1]
ser_col <- seriate(as.matrix(t(cluster_mat)), method = "PCA")
x_order <- get_order(ser_col)
x_order <- colnames(cluster_mat)[get_order(ser_col)]

all_data <- all_data %>%
  mutate(
    Description = factor(Description, 
                         levels = c(x_order, setdiff(unique(Description), x_order)))
  )
################################################################################
all_data <- all_data %>% filter(!Description %in% c("Drug metabolism - other enzymes", 
                                                    "Carbon fixation in photosynthetic organisms",""))
all_data$p.adjust <- -log10(all_data$p.adjust)
all_data$Cluster <- factor(all_data$Cluster, levels = rev(all_clusters))

combined_plot <- ggplot(all_data, aes(x = Description, y = Cluster)) +
  geom_point(
    aes(size = Count * 1.15),
    position = position_nudge(x = 0.08, y = -0.08),
    color = "gray20",
    alpha = 0.4,
    shape = 19,
    show.legend = FALSE
  ) +
  geom_point(
    aes(size = Count, fill = p.adjust),
    color = "black",
    shape = 21,
    stroke = 0.8,
    alpha = 1 
  )+
  facet_wrap(~ Kingdom, ncol = 1, strip.position = "right") +
  scale_fill_gradientn(
    name = expression(-log[10](p.adjust)),
    colors = c("#FFFFE0", "#E41A1C"),
    values = scales::rescale(c(0, 5, 10, 20, 50)),
    limits = c(0, 50)
  ) + 
  scale_y_discrete(drop = FALSE) +
  scale_size_continuous(
    range = c(1, 4),
    breaks = c(50, 100, 150,200),
    name = "KO Count"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 8, color = 'black', hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = 'black', hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8, color = 'black'), 
    strip.background = element_rect(fill = "lightgray", color = "black"),
    panel.grid.major = element_line(color = "grey90",linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 0.6),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
    legend.title = element_text(size = 8, color = 'black'), 
    legend.text = element_text(size = 8,  color = 'black'), 
    legend.box.margin = margin(0, 0, 0, -20),
    legend.key.size = unit(0.25, 'cm'),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.position = "bottom"
  )+
  guides(
    fill = guide_colorbar(label.theme = element_text(size = 6, color = 'black', vjust = 0.5)),
    size = guide_legend(nrow = 2)
  ) +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL)

width = 12
height = 24
ggsave(paste0(dir_name,"/conversed_positive_edge_enrichment.pdf"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(dir_name,"/conversed_positive_edge_enrichment.png"), combined_plot, width = width, height = height, units = "cm")
#-------------------------------------------------------------------------------



###########################
tmp_df2 <- final_df %>%
  select(
    PPPPPPPPPP,
    PUUUUUUUUU,
    UPUUUUUUUU,
    UUPUUUUUUU,
    UUUPUUUUUU,
    UUUUPUUUUU,
    UUUUUPUUUU,
    UUUUUUPUUU,
    UUUUUUUPUU,
    UUUUUUUUPU,
    UUUUUUUUUP
  )


prefix_map <- c("b" = "Bacteria", "f" = "Fungi", "p" = "Protists")


pathway2ko <- read.csv('00-data/pathway2ko.csv')
pathway2superpathway <- read.csv('00-data/pathway2superpathway.csv')
for (prefix in names(prefix_map)) {
  final_df_filter <- as.data.frame(lapply(tmp_df2, function(x) ifelse(grepl(paste0("^", prefix), x, perl = TRUE), x, NA)))
  final_df_filter <- as.data.frame(lapply(final_df_filter, function(x) {
    sub(paste0("^", prefix), "", x, perl = TRUE)
  }))
  final_df_filter <- final_df_filter[rowSums(is.na(final_df_filter)) < ncol(final_df_filter), ]
  final_df_filter <- as.data.frame(
    lapply(final_df_filter, function(x) na.locf(x, fromLast = TRUE, na.rm = FALSE))
  )
  col_unique <- lapply(final_df_filter, unique)
  max_length <- max(sapply(col_unique, length))
  final_df_col_unique <- as.data.frame(
    lapply(col_unique, function(col) c(col, rep(NA, max_length - length(col))))
  )
  ko_res <- compareCluster(final_df_col_unique, fun = 'enricher', TERM2GENE = pathway2ko)
  ko_res@compareClusterResult <- merge(ko_res@compareClusterResult, pathway2superpathway, by.x = 'ID', by.y = 'Pathway', all.x = T)
  write.csv(ko_res@compareClusterResult, paste0(dir_name, '/', prefix_map[prefix], '_conversed_and_unique_positive_edge_enrichment_results.csv'), row.names = F)
}

data_list <- list()
plot_list <- list()
kingdoms <- c("Bacteria", "Fungi", "Protists")
#kingdoms <- c("Bacteria", "Fungi")

all_clusters2 <- c("PPPPPPPPPP","PUUUUUUUUU","UPUUUUUUUU","UUPUUUUUUU","UUUPUUUUUU",
                  "UUUUPUUUUU","UUUUUPUUUU","UUUUUUPUUU","UUUUUUUPUU","UUUUUUUUPU",
                  "UUUUUUUUUP")



for (current_kingdom in kingdoms) {
  data_df <- read.csv(paste0(dir_name, '/', current_kingdom, '_conversed_and_unique_positive_edge_enrichment_results.csv'))
  data_df <- data_df %>%
    mutate(
      GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
  #data_df <- data_df %>% filter(p.adjust < 0.01 & GeneRatio > 0.01 & FoldEnrichment > 2 )
  data_df$Kingdom <- current_kingdom
  
  data_list[[current_kingdom]] <- data_df
}


all_data <- bind_rows(data_list)
all_data <- all_data %>%
  group_by(Description) %>%
  filter(n_distinct(Kingdom) >= 2) %>%
  ungroup()
################################################################################
data_ko_df <- all_data %>% filter(Kingdom == "Bacteria")
data_ko_df <- setDT(data_ko_df)
cluster_mat <- dcast(data_ko_df, Cluster ~ Description, value.var = "Count", fill = 0)
rownames(cluster_mat) <- cluster_mat$Cluster
cluster_mat <- cluster_mat[, -1]
ser_col <- seriate(as.matrix(t(cluster_mat)), method = "PCA")
x_order <- get_order(ser_col)
x_order <- colnames(cluster_mat)[get_order(ser_col)]

all_data <- all_data %>%
  mutate(
    Description = factor(Description, 
                         levels = c(x_order, setdiff(unique(Description), x_order)))
  )
################################################################################
all_data <- all_data %>% filter(!Description %in% c("Drug metabolism - other enzymes", 
                                                    "Carbon fixation in photosynthetic organisms",""))
all_data$p.adjust <- -log10(all_data$p.adjust)
all_data$Cluster <- factor(all_data$Cluster, levels = rev(all_clusters2))

combined_plot <- ggplot(all_data, aes(x = Description, y = Cluster)) +
  geom_point(
    aes(size = Count * 1.15),
    position = position_nudge(x = 0.08, y = -0.08),
    color = "gray20",
    alpha = 0.4,
    shape = 19,
    show.legend = FALSE
  ) +
  geom_point(
    aes(size = Count, fill = p.adjust),
    color = "black",
    shape = 21,
    stroke = 0.8,
    alpha = 1 
  )+
  facet_wrap(~ Kingdom, ncol = 1, strip.position = "right") +
  scale_fill_gradientn(
    name = expression(-log[10](p.adjust)),
    colors = c("#FFFFE0", "#E41A1C"),
    values = scales::rescale(c(0, 5, 10, 20, 50)),
    limits = c(0, 50)
  ) + 
  scale_y_discrete(drop = FALSE) +
  scale_size_continuous(
    range = c(1, 4),
    breaks = c(50, 100, 150,200),
    name = "KO Count"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 8, color = 'black', hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = 'black', hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8, color = 'black'), 
    strip.background = element_rect(fill = "lightgray", color = "black"),
    panel.grid.major = element_line(color = "grey90",linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 0.6),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
    legend.title = element_text(size = 8, color = 'black'), 
    legend.text = element_text(size = 8,  color = 'black'), 
    legend.box.margin = margin(0, 0, 0, -20),
    legend.key.size = unit(0.25, 'cm'),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.position = "bottom"
  )+
  guides(
    fill = guide_colorbar(label.theme = element_text(size = 6, color = 'black', vjust = 0.5)),
    size = guide_legend(nrow = 2)
  ) +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL)

width = 15
height = 24
ggsave(paste0(dir_name,"/conversed_and_unique_positive_edge_enrichment.pdf"), combined_plot, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(dir_name,"/conversed_and_unique_positive_edge_enrichment.png"), combined_plot, width = width, height = height, units = "cm")
#-------------------------------------------------------------------------------






