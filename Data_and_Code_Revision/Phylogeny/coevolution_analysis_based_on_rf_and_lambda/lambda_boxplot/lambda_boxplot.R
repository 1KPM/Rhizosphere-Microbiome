### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(dplyr)
library(rstatix)
library(agricolae)
library(ggplot2)
library(purrr)
library(patchwork)

### KO数据
superpathway_color <- read.csv("./data/All_superpathway_color.csv")
kingdoms <- c("Bacteria", "Fungi", "Protists")
all_sp_data <- list()

ko_df <- read.csv("./data/kegg_pathway_for_1KPM.csv")
ko_df$pathway <- sub("\\[.*", "", ko_df$pathway)

for (kingdom in kingdoms) {
  file_path <- paste0("./data/", kingdom, "_core0.2_qa_log_phylogenetic_signal.csv")
  data <- read.csv(file_path)
  names(data)[1] <- "KO"
  ko <- ko_df %>% filter(level2 %in% superpathway_color$superpathway)
  
  fin_df <- data %>%
    mutate(Phylogenetically = ifelse(padj < 0.05, "conserved", "labile")) %>%
    inner_join(ko %>% select(KO, level2), by = "KO") %>%
    select(level2, lambda, Phylogenetically) %>%
    rename(Group = level2, Value = lambda, Roles = Phylogenetically) %>%
    left_join(superpathway_color %>% select(superpathway, Color),
              by = c("Group" = "superpathway")) %>%
    mutate(
      Kingdom = kingdom,
      point_color = ifelse(Roles == "conserved",
                           coalesce(Color, "red"),
                           "grey80")
    ) %>%
    select(-Color)
  
  fin_df$Group <- factor(fin_df$Group, levels = sort(unique(fin_df$Group)))
  
  mean_df <- aggregate(Value ~ Group, data = fin_df, FUN = mean)
  max_df <- aggregate(Value ~ Group, data = fin_df, FUN = max)
  
  dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
  dunn_res_df <- data.frame(dunn_res[-1])
  write.csv(dunn_res_df, paste0("./data/KO_", kingdom, "_dunn_res.csv"))
  
  n <- nrow(mean_df)
  pvalue_df <- matrix(1, ncol = n, nrow = n)
  k <- 0
  for (i in 1:(n - 1)) { 
    for (j in (i + 1):n) { 
      k <- k + 1
      pvalue_df[i, j] <- dunn_res_df$p.adj[k]
      pvalue_df[j, i] <- dunn_res_df$p.adj[k]
    }
  }
  
  letter_df <- orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = FALSE)
  letter_df <- letter_df[levels(fin_df$Group), ]
  
  max_df$groups <- letter_df$groups
  max_df$Kingdom <- kingdom
  
  all_sp_data[[kingdom]] <- list(fin_df = fin_df, max_df = max_df)
}

### ASV数据
###############################################################################
bac_asv <- read.csv("./data/Bacteria_absolute_format_core_ASV_phylogenetic_signal.csv")
fun_asv <- read.csv("./data/Fungi_absolute_format_core_ASV_phylogenetic_signal.csv")
pro_asv <- read.csv("./data/Protists_absolute_format_core_ASV_phylogenetic_signal.csv")
all_df <- rbind.data.frame(bac_asv, fun_asv, pro_asv)
names(all_df)[1] <- "FeatureID"
  
all_order_data <- list()
asv_df <- read.csv("./data/All_core_ASV_taxonomy.csv")

for (kingdom in kingdoms) {
  data <- all_df %>% filter(Kingdom == kingdom)
  order_color <- read.csv(paste0("./data/", kingdom, "_top10order_colors.csv"))
  
  asv <- asv_df %>% filter(Order %in% order_color$Taxa)
  
  fin_df <- data %>%
    filter(FeatureID %in% asv$FeatureID) %>%
    mutate(Phylogenetically = ifelse(padj < 0.05, "conserved", "labile")) %>%
    inner_join(asv %>% select(FeatureID, Order), by = "FeatureID") %>%
    select(Order, lambda, Phylogenetically) %>%
    rename(Group = Order, Value = lambda, Roles = Phylogenetically) %>%
    left_join(order_color %>% select(Taxa, Color), by = c("Group" = "Taxa")) %>%
    mutate(
      Kingdom = kingdom,
      point_color = ifelse(Roles == "conserved", coalesce(Color, "red"), "grey80")
    ) %>%
    select(-Color)
  
  fin_df$Group <- factor(fin_df$Group, levels = order_color$Taxa)
  
  
  mean_df <- aggregate(Value ~ Group, data = fin_df, FUN = mean)
  max_df <- aggregate(Value ~ Group, data = fin_df, FUN = max)
  
  dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
  dunn_res_df <- data.frame(dunn_res[-1])
  write.csv(dunn_res_df, paste0("./data/ASV_", kingdom, "_dunn_res.csv"))
  
  n <- nrow(mean_df)
  pvalue_df <- matrix(1, ncol = n, nrow = n)
  k <- 0
  for (i in 1:(n - 1)) { 
    for (j in (i + 1):n) { 
      k <- k + 1
      pvalue_df[i, j] <- dunn_res_df$p.adj[k]
      pvalue_df[j, i] <- dunn_res_df$p.adj[k]
    }
  }
  
  letter_df <- orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = FALSE)
  letter_df <- letter_df[levels(fin_df$Group), ]
  
  max_df$groups <- letter_df$groups
  max_df$Kingdom <- kingdom
  
  all_order_data[[kingdom]] <- list(fin_df = fin_df, max_df = max_df)
}



###############################################################################
common_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 6, 
      color = 'black',
      angle = 45, 
      hjust = 1, 
      vjust = 1
    ),
    axis.text.y = element_text(size = 6, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6, color = 'black'),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey80"),
    strip.text = element_text(size = 6),
    plot.margin = margin(l = 5, r = 5, unit = "mm"),
    plot.tag = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8) 
  )

plots <- list()

# Bacteria
plots[["Bacteria_Order"]] <- ggplot(all_order_data[["Bacteria"]]$fin_df, aes(x = Group, y = Value)) +
  geom_point(aes(color = point_color), position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = all_order_data[["Bacteria"]]$max_df,
    aes(y = Value * 1.1, label = groups),
    size = 1.7
  ) +
  scale_color_identity() +
  labs(
    x = "Order",
    y = expression("Phylogenetic Signal ("*lambda*")"),
    title = "Bacterial ASVs",
    tag = "a"
  ) +
  common_theme

plots[["Bacteria_Superpathway"]] <- ggplot(all_sp_data[["Bacteria"]]$fin_df, aes(x = Group, y = Value)) +
  geom_point(aes(color = point_color), position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = all_sp_data[["Bacteria"]]$max_df,
    aes(y = Value * 1.1, label = groups),
    size = 1.7
  ) +
  scale_color_identity() +
  labs(
    x = "Superpathway",
    y = expression("Phylogenetic Signal ("*lambda*")"),
    title = "Bacterial KOs",
    tag = "d"
  ) +
  common_theme

# Fungi
plots[["Fungi_Order"]] <- ggplot(all_order_data[["Fungi"]]$fin_df, aes(x = Group, y = Value)) +
  geom_point(aes(color = point_color), position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = all_order_data[["Fungi"]]$max_df,
    aes(y = Value * 1.1, label = groups),
    size = 1.7
  ) +
  scale_color_identity() +
  labs(
    x = "Order",
    y = expression("Phylogenetic Signal ("*lambda*")"),
    title = "Fungal  ASVs",
    tag = "b"
  ) +
  common_theme

plots[["Fungi_Superpathway"]] <- ggplot(all_sp_data[["Fungi"]]$fin_df, aes(x = Group, y = Value)) +
  geom_point(aes(color = point_color), position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = all_sp_data[["Fungi"]]$max_df,
    aes(y = Value * 1.1, label = groups),
    size = 1.7
  ) +
  scale_color_identity() +
  labs(
    x = "Superpathway",
    y = expression("Phylogenetic Signal ("*lambda*")"),
    title = "Fungal KOs",
    tag = "e"
  ) +
  common_theme

# Protists
plots[["Protists_Order"]] <- ggplot(all_order_data[["Protists"]]$fin_df, aes(x = Group, y = Value)) +
  geom_point(aes(color = point_color), position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = all_order_data[["Protists"]]$max_df,
    aes(y = Value * 1.1, label = groups),
    size = 1.7
  ) +
  scale_color_identity() +
  labs(
    x = "Order",
    y = expression("Phylogenetic Signal ("*lambda*")"),
    title = "Protistan ASVs",
    tag = "c"
  ) +
  common_theme

plots[["Protists_Superpathway"]] <- ggplot(all_sp_data[["Protists"]]$fin_df, aes(x = Group, y = Value)) +
  geom_point(aes(color = point_color), position = position_jitter(width = 0.1), alpha = 0.6, size = 1) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = all_sp_data[["Protists"]]$max_df,
    aes(y = Value * 1.1, label = groups),
    size = 1.7
  ) +
  scale_color_identity() +
  labs(
    x = "Superpathway",
    y = expression("Phylogenetic Signal ("*lambda*")"),
    title = "Protistan KOs",
    tag = "f"
  ) +
  common_theme

plot_list <- list(
  plots[["Bacteria_Order"]],
  plots[["Bacteria_Superpathway"]],
  plots[["Fungi_Order"]],
  plots[["Fungi_Superpathway"]],
  plots[["Protists_Order"]],
  plots[["Protists_Superpathway"]]
)

combined_plot <- wrap_plots(
  plot_list,
  nrow = 3,
  ncol = 2,
  widths = c(1, 2)
)


ggsave("./lambda_boxplot.pdf", combined_plot, width = 20, height = 25, units = "cm",  dpi = 600)
ggsave("./lambda_boxplot.png", combined_plot, width = 20, height = 25, units = "cm",  dpi = 600)


