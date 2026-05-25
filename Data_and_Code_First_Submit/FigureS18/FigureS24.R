pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(dplyr)
library(rstatix)
library(agricolae)
library(ggplot2)
library(purrr)
library(patchwork)

##Total
data_bac <- read.csv("./data/Bacteria_KO_phylogenetic_signal_qa(log).csv") %>% mutate(Kingdom = "Bacteria")
data_fun <- read.csv("./data/Fungi_KO_phylogenetic_signal_qa(log).csv") %>% mutate(Kingdom = "Fungi")
data_pro <- read.csv("./data/Protist_KO_phylogenetic_signal_qa(log).csv") %>% mutate(Kingdom = "Protist")
all_data <- rbind(data_bac,data_fun,data_pro)

fin_df <- all_data %>%
  mutate(Phylogenetically = ifelse(P < 0.05, "conserved", "labile")) %>%
  select(Kingdom, lambda, Phylogenetically) %>%
  rename(Group = Kingdom, Value = lambda, Roles = Phylogenetically) %>%
  mutate(point_color = case_when(
    # 当 Roles 为 "conserved" 时，根据 Group 分配颜色
    Roles == "conserved" & Group == "Bacteria" ~ "#7FC97F",
    Roles == "conserved" & Group == "Fungi" ~ "#BEAED4",
    Roles == "conserved" & Group == "Protist" ~ "#FDC086",
    Roles == "labile" ~ "grey80",
  ))

fin_df$Group <- factor(fin_df$Group, levels = sort(unique(fin_df$Group)))
mean_df <- aggregate(Value ~ Group, data = fin_df, FUN = mean)
max_df <- aggregate(Value ~ Group, data = fin_df, FUN = max)

dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
dunn_res_df <- data.frame(dunn_res[-1])
write.csv(dunn_res_df,"./data/KO_kingdom_dunn_res.csv",row.names = F)
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

p2 <- ggplot(fin_df, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color),
    position = position_jitter(width = 0.1),
    alpha = 0.6, 
    size = 1
  ) +
  geom_boxplot(width = 0.25, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.4, alpha = 0.2, na.rm = TRUE) +  
  geom_text(
    data = max_df,
    aes(y = Value * 1.1, label = groups),
    size = 4
  ) +
  scale_color_identity() +
  labs(
    x = "",
    y = expression("Phylogenetic Signal ("*lambda*")")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 10, 
      color = 'black',
      angle = 45, 
      hjust = 1, 
      vjust = 1
    ),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey80"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = margin(l = 15, r = 15,unit = "mm"),
  )



##Kingdom-superpathway
superpathway_color <- read.csv("./data/All_superpathway_color.csv")
kingdoms <- c("Bacteria", "Fungi", "Protist")
all_data <- list()

ko_df <- read.csv("./data/kegg_pathway_for_1KPM.csv")
ko_df$pathway <- sub("\\[.*", "", ko_df$pathway)

all_order_colors <- list()

for (kingdom in kingdoms) {
  file_path <- paste0("./data/", kingdom, "_KO_phylogenetic_signal_qa(log).csv")
  data <- read.csv(file_path)
  all_order_colors[[kingdom]] <- superpathway_color$superpathway
  ko <- ko_df %>% filter(level2 %in% superpathway_color$superpathway)
  fin_df <- data %>%
    mutate(Phylogenetically = ifelse(P < 0.05, "conserved", "labile")) %>%
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
  
  all_data[[kingdom]] <- list(fin_df = fin_df, max_df = max_df)
}

combined_fin <- map_dfr(all_data, ~ .x$fin_df)
combined_max <- map_dfr(all_data, ~ .x$max_df)

combined_fin <- combined_fin %>%
  mutate(
    Group = factor(Group, levels = unique(unlist(all_order_colors)))
  )

p4 <- ggplot(combined_fin, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color),
    position = position_jitter(width = 0.1),
    alpha = 0.6, 
    size = 1
  ) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = combined_max,
    aes(y = Value * 1.1, label = groups),
    size = 4
  ) +
  scale_color_identity() +
  facet_wrap(~ Kingdom, scales = "free_x", nrow = 1,strip.position = "bottom") +
  labs(
    x = "Superpathway",
    y = expression("Phylogenetic Signal ("*lambda*")")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 10, 
      color = 'black',
      angle = 45, 
      hjust = 1, 
      vjust = 1
    ),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey80"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = margin(l = 15, r = 15,unit = "mm"),
  )


##Total
all_data <- read.csv("./data/absolute_format_core_ASV_phylogenetic_signal.csv")
data_bac <- all_data %>% filter(Kingdom == "Bacteria")
data_fun <- all_data %>% filter(Kingdom == "Fungi")
data_pro <- all_data %>% filter(Kingdom == "Protist")

fin_df <- all_data %>%
  mutate(Phylogenetically = ifelse(P < 0.05, "conserved", "labile")) %>%
  select(Kingdom, lambda, Phylogenetically) %>%
  rename(Group = Kingdom, Value = lambda, Roles = Phylogenetically) %>%
  mutate(point_color = case_when(
    # 当 Roles 为 "conserved" 时，根据 Group 分配颜色
    Roles == "conserved" & Group == "Bacteria" ~ "#7FC97F",
    Roles == "conserved" & Group == "Fungi" ~ "#BEAED4",
    Roles == "conserved" & Group == "Protist" ~ "#FDC086",
    Roles == "labile" ~ "grey80",
  ))

fin_df$Group <- factor(fin_df$Group, levels = sort(unique(fin_df$Group)))
mean_df <- aggregate(Value ~ Group, data = fin_df, FUN = mean)
max_df <- aggregate(Value ~ Group, data = fin_df, FUN = max)

dunn_res <- dunn_test(fin_df, Value ~ Group, p.adjust.method = 'fdr')
dunn_res_df <- data.frame(dunn_res[-1])
write.csv(dunn_res_df,"./data/ASV_kingdom_dunn_res.csv",row.names = F)
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

p1 <- ggplot(fin_df, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color),
    position = position_jitter(width = 0.1),
    alpha = 0.6, 
    size = 1
  ) +
  geom_boxplot(width = 0.25, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.4, alpha = 0.2, na.rm = TRUE) +  
  geom_text(
    data = max_df,
    aes(y = Value * 1.1, label = groups),
    size = 4
  ) +
  scale_color_identity() +
  labs(
    x = "",
    y = expression("Phylogenetic Signal ("*lambda*")")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 10, 
      color = 'black',
      angle = 45, 
      hjust = 1, 
      vjust = 1
    ),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey80"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = margin(l = 15, r = 15,unit = "mm"),
  )



##Kingdom-Order
all_df <- read.csv("./data/absolute_format_core_ASV_phylogenetic_signal.csv")
kingdoms <- c("Bacteria", "Fungi", "Protist")
all_data <- list()
asv_df <- read.csv("./data/All_core_ASV_taxonomy.csv")
all_order_colors <- list()

for (kingdom in kingdoms) {
  data <- all_df %>% filter(Kingdom == kingdom)
  order_color <- read.csv(paste0("./data/", kingdom, "_top10order_colors.csv"))
  all_order_colors[[kingdom]] <- order_color$Taxa
  asv <- asv_df %>% filter(Order %in% order_color$Taxa)
  
  fin_df <- data %>%
    filter(FeatureID %in% asv$FeatureID) %>%
    mutate(Phylogenetically = ifelse(P < 0.05, "conserved", "labile")) %>%
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
  
  all_data[[kingdom]] <- list(fin_df = fin_df, max_df = max_df)
}

combined_fin <- map_dfr(all_data, ~ .x$fin_df)
combined_max <- map_dfr(all_data, ~ .x$max_df)

# 确保每个分面的 Group 顺序正确
combined_fin <- combined_fin %>%
  mutate(
    Group = factor(Group, levels = unique(unlist(all_order_colors)))
  )

p3 <- ggplot(combined_fin, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color),
    position = position_jitter(width = 0.1),
    alpha = 0.6, 
    size = 1
  ) +
  geom_boxplot(width = 0.4, alpha = 0.2, na.rm = TRUE) +
  geom_text(
    data = combined_max,
    aes(y = Value * 1.1, label = groups),
    size = 4
  ) +
  scale_color_identity() +
  facet_wrap(~ Kingdom, scales = "free_x", nrow = 1, strip.position = "bottom") +
  labs(
    x = "Order",
    y = expression("Phylogenetic Signal ("*lambda*")")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 10, 
      color = 'black',
      angle = 45, 
      hjust = 1, 
      vjust = 1
    ),
    axis.text.y = element_text(size = 10, color = 'black'),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey80"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = margin(l = 15, r = 15, unit = "mm"),
  )


combined <- 
  (plot_spacer() + p1 + p2 + plot_spacer() + 
     plot_layout(widths = c(1, 8, 8, 1))) /
  p3 / 
  p4 +
  plot_layout(
    ncol = 1,
    heights = c(1, 1, 1)
  ) + 
  plot_annotation(tag_levels = 'A')
ggsave("./FigureS24.pdf", combined, width = 14, height = 16, dpi = 600)
ggsave("./FigureS24.png", combined, width = 14, height = 16, dpi = 600)
ggsave("./FigureS24.tiff", combined, width = 14, height = 16, dpi = 600)