pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(agricolae)
library(patchwork)

data <- read.csv("./data/FigureS14-Cross_kingdom_network_hub_info.csv")
names(data)[1:4] <- c("KO", "Degree", "Closeness", "Betweenness") 
ko_df <- read.csv("./data/kegg_pathway_for_1KPM.csv")
ko_df$pathway <- sub("\\[.*", "", ko_df$pathway)
color_manual1 <- colorRampPalette(brewer.pal(9, 'Set1'))(25)[1:24]


categories <- c(b = "Bacteria", f = "Fungi", p = "Protist")
result_list <- list()

for (prefix in names(categories)) {
  sub_data <- data %>%
    filter(grepl(paste0("^", prefix), KO)) %>%
    mutate(
      KO_clean = sub(paste0("^", prefix), "", KO),
      Category = categories[prefix]
    )
  
  merged_data <- left_join(sub_data, 
                           ko_df[, c("KO", "pathway")], 
                           by = c("KO_clean" = "KO"))
  
  result_list[[categories[prefix]]] <- merged_data %>%
    select(pathway, Degree, roles)
}
group_counts <- result_list[["Fungi"]] %>%
  group_by(pathway) %>%
  summarise(count = n())
filtered_df <- result_list[["Fungi"]] %>%
  filter(pathway %in% group_counts$pathway[group_counts$count >= 10])
Fungi_fin_df1 <- filtered_df %>%
  filter(!is.na(pathway))

names(Fungi_fin_df1) <- c('Group', 'Value' ,'Roles')
Fungi_max_df1 <- aggregate(Fungi_fin_df1['Value'], by = list(Group = Fungi_fin_df1$Group), FUN = max)
Fungi_mean_df1 <- aggregate(Fungi_fin_df1['Value'], by = list(Group = Fungi_fin_df1$Group), FUN = mean)

n <- nrow(Fungi_mean_df1)
Fungi_dunn_res1 <- dunn_test(Fungi_fin_df1, Value ~ Group, p.adjust.method = 'fdr')
Fungi_dunn_res_df1 <- data.frame(Fungi_dunn_res1[-1])

Fungi_pvalue_df1 <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
  for(j in (i + 1):n){ 
    k <- k + 1
    Fungi_pvalue_df1[i,j] <- Fungi_dunn_res_df1$p.adj[k]
    Fungi_pvalue_df1[j,i] <- Fungi_dunn_res_df1$p.adj[k]
  }
}

Fungi_letter_df1 <- orderPvalue(Fungi_mean_df1$Group, Fungi_mean_df1$Value, 0.05, Fungi_pvalue_df1, console = TRUE)
Fungi_fin_df1$Group <- as.factor(Fungi_fin_df1$Group)
Fungi_letter_df1 <- Fungi_letter_df1[levels(Fungi_fin_df1$Group),]
Fungi_letter_vector1 <- Fungi_letter_df1$groups

group_levels1 <- levels(Fungi_fin_df1$Group)
color_mapping1 <- setNames(color_manual1, group_levels1)

Fungi_fin_df1 <- Fungi_fin_df1 %>%
  mutate(
    point_color = ifelse(
      Roles == "Peripherals",
      "grey80",
      color_mapping1[as.character(Group)]
    )
  )

p2 <- ggplot(Fungi_fin_df1, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color, group = Group),
    position = position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8, size = 1
  ) +
  scale_color_identity() + 
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
  geom_text(data = Fungi_max_df1, aes(y = Value * 1.2, label = Fungi_letter_vector1), position = position_dodge(0.9), size = 7 / 2.835) + 
  labs(
    title = '',
    x = '',
    y = 'Fungi',
  ) + 
theme_bw() + theme(
    text = element_text(color = "black", size = 6),
    plot.title = element_text(size = 7, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 7),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6, color = "black"),
    axis.text.x = element_blank(),
    strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1,"cm"),
    legend.key.size = unit(0.25, "cm"),
    legend.position = "right"
)


Bacteria_fin_df1 <- result_list[["Bacteria"]] %>%
  filter(pathway %in% group_levels1)

names(Bacteria_fin_df1) <- c('Group', 'Value' ,'Roles')
Bacteria_max_df1 <- aggregate(Bacteria_fin_df1['Value'], by = list(Group = Bacteria_fin_df1$Group), FUN = max)
Bacteria_mean_df1 <- aggregate(Bacteria_fin_df1['Value'], by = list(Group = Bacteria_fin_df1$Group), FUN = mean)

n <- nrow(Bacteria_mean_df1)
Bacteria_dunn_res1 <- dunn_test(Bacteria_fin_df1, Value ~ Group, p.adjust.method = 'fdr')
Bacteria_dunn_res_df1 <- data.frame(Bacteria_dunn_res1[-1])

Bacteria_pvalue_df1 <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
  for(j in (i + 1):n){ 
    k <- k + 1
    Bacteria_pvalue_df1[i,j] <- Bacteria_dunn_res_df1$p.adj[k]
    Bacteria_pvalue_df1[j,i] <- Bacteria_dunn_res_df1$p.adj[k]
  }
}

Bacteria_letter_df1 <- orderPvalue(Bacteria_mean_df1$Group, Bacteria_mean_df1$Value, 0.05, Bacteria_pvalue_df1, console = TRUE)
Bacteria_fin_df1$Group <- as.factor(Bacteria_fin_df1$Group)
Bacteria_letter_df1 <- Bacteria_letter_df1[levels(Bacteria_fin_df1$Group),]
Bacteria_letter_vector1 <- Bacteria_letter_df1$groups

Bacteria_fin_df1 <- Bacteria_fin_df1 %>%
  mutate(
    point_color = ifelse(
      Roles == "Peripherals",
      "grey80",
      color_mapping1[as.character(Group)]
    )
  )

p1 <- ggplot(Bacteria_fin_df1, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color, group = Group),
    position = position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8, size = 1
  ) +
  scale_color_identity() + 
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
  geom_text(data = Bacteria_max_df1, aes(y = Value * 1.2, label = Bacteria_letter_vector1), 
            position = position_dodge(0.9), size = 7 / 2.835) + 
  labs(
    title = 'Pathway',
    x = '',
    y = 'Bacteria',
  ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )


full_df <- data.frame(pathway = group_levels1)

Protist_fin_df1 <- full_df %>%
  left_join(result_list[["Protist"]], by = "pathway") %>%
  mutate(
    Degree = ifelse(is.na(Degree), 0, Degree),
    roles = ifelse(is.na(roles), "Peripherals", roles)
  )

names(Protist_fin_df1) <- c('Group', 'Value' ,'Roles')
Protist_max_df1 <- aggregate(Protist_fin_df1['Value'], by = list(Group = Protist_fin_df1$Group), FUN = max)
Protist_mean_df1 <- aggregate(Protist_fin_df1['Value'], by = list(Group = Protist_fin_df1$Group), FUN = mean)

n <- nrow(Protist_mean_df1)
Protist_dunn_res1 <- dunn_test(Protist_fin_df1, Value ~ Group, p.adjust.method = 'fdr')
Protist_dunn_res_df1 <- data.frame(Protist_dunn_res1[-1])

Protist_pvalue_df1 <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
  for(j in (i + 1):n){ 
    k <- k + 1
    Protist_pvalue_df1[i,j] <- Protist_dunn_res_df1$p.adj[k]
    Protist_pvalue_df1[j,i] <- Protist_dunn_res_df1$p.adj[k]
  }
}

Protist_letter_df1 <- orderPvalue(Protist_mean_df1$Group, Protist_mean_df1$Value, 0.05, Protist_pvalue_df1, console = TRUE)
Protist_fin_df1$Group <- as.factor(Protist_fin_df1$Group)
Protist_letter_df1 <- Protist_letter_df1[levels(Protist_fin_df1$Group),]
Protist_letter_vector1 <- Protist_letter_df1$groups


Protist_fin_df1 <- Protist_fin_df1 %>%
  mutate(
    point_color = ifelse(
      Roles == "Peripherals",
      "grey80",
      color_mapping1[as.character(Group)]
    )
  )

p3 <- ggplot(Protist_fin_df1, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color, group = Group),
    position = position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8, size = 1
  ) +
  scale_color_identity() + 
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
  geom_text(data = Protist_max_df1, aes(y = Value * 1.2, label = Protist_letter_vector1), 
            position = position_dodge(0.9), size = 7 / 2.835) + 
  labs(
    title = '',
    x = '',
    y = 'Protist',
  ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )


categories <- c(b = "Bacteria", f = "Fungi", p = "Protist")
result_list <- list()
fill_color <- read.csv("./data/superpathway_color.csv")
color_manual2 <- fill_color$fill_color2
names(color_manual2) <- fill_color$superpathway
for (prefix in names(categories)) {
  sub_data <- data %>%
    filter(grepl(paste0("^", prefix), KO)) %>%
    mutate(
      KO_clean = sub(paste0("^", prefix), "", KO),
      Category = categories[prefix]
    )
  
  merged_data <- left_join(sub_data, 
                           ko_df[, c("KO", "level2")], 
                           by = c("KO_clean" = "KO"))
  
  result_list[[categories[prefix]]] <- merged_data %>%
    select(level2, Degree, roles)
  
}


Bacteria_fin_df <- result_list[["Bacteria"]] %>%
  mutate(
    level2 = case_when(
      level2 == "Cellular community - prokaryotes" ~ "Cellular community",
      TRUE ~ level2
    ),
    level2 = replace_na(level2, "Others")
  ) %>%
  rename(Group = level2, Value = Degree, Roles = roles)

Bacteria_fin_df <- Bacteria_fin_df %>% filter(Group != "Others")
names(Bacteria_fin_df) <- c('Group', 'Value' ,'Roles')
#Bacteria_fin_df$Group <- factor(Bacteria_fin_df$Group,levels = c(sort(unique(Bacteria_fin_df$Group[Bacteria_fin_df$Group != "Others"])), "Others"))
Bacteria_fin_df$Group <- factor(Bacteria_fin_df$Group,levels = sort(unique(Bacteria_fin_df$Group)))
Bacteria_max_df <- aggregate(Bacteria_fin_df['Value'], by = list(Group = Bacteria_fin_df$Group), FUN = max)
Bacteria_mean_df <- aggregate(Bacteria_fin_df['Value'], by = list(Group = Bacteria_fin_df$Group), FUN = mean)

n <- nrow(Bacteria_mean_df)
Bacteria_dunn_res <- dunn_test(Bacteria_fin_df, Value ~ Group, p.adjust.method = 'fdr')
Bacteria_dunn_res_df <- data.frame(Bacteria_dunn_res[-1])

Bacteria_pvalue_df <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
  for(j in (i + 1):n){ 
    k <- k + 1
    Bacteria_pvalue_df[i,j] <- Bacteria_dunn_res_df$p.adj[k]
    Bacteria_pvalue_df[j,i] <- Bacteria_dunn_res_df$p.adj[k]
  }
}

Bacteria_letter_df <- orderPvalue(Bacteria_mean_df$Group, Bacteria_mean_df$Value, 0.05, Bacteria_pvalue_df, console = TRUE)
Bacteria_fin_df$Group <- as.factor(Bacteria_fin_df$Group)
Bacteria_letter_df <- Bacteria_letter_df[levels(Bacteria_fin_df$Group),]
Bacteria_letter_vector <- Bacteria_letter_df$groups

color_mapping <- color_manual2

Bacteria_fin_df <- Bacteria_fin_df %>%
  mutate(
    point_color = ifelse(
      Roles == "Peripherals",
      "grey80",
      color_mapping[as.character(Group)]
    )
  )

p4 <- ggplot(Bacteria_fin_df, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color, group = Group),
    position = position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8, size = 1
  ) +
  scale_color_identity() + 
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
  geom_text(data = Bacteria_max_df, aes(y = Value * 1.2, label = Bacteria_letter_vector), 
            position = position_dodge(0.9), size = 7 / 2.835) + 
  labs(
    title = 'Superpathway',
    x = '',
    y = 'Bacteria',
  ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )

Fungi_fin_df <- result_list[["Fungi"]] %>%
  mutate(
    level2 = case_when(
      level2 == "Cellular community - prokaryotes" ~ "Cellular community",
      TRUE ~ level2
    ),
    level2 = replace_na(level2, "Others")
  ) %>%
  rename(Group = level2, Value = Degree, Roles = roles)
Fungi_fin_df <- Fungi_fin_df %>% filter(Group != "Others")
names(Fungi_fin_df) <- c('Group', 'Value' ,'Roles')
#Fungi_fin_df$Group <- factor(Fungi_fin_df$Group,levels = c(sort(unique(Fungi_fin_df$Group[Fungi_fin_df$Group != "Others"])), "Others"))
Fungi_fin_df$Group <- factor(Fungi_fin_df$Group,levels = sort(unique(Fungi_fin_df$Group)))
Fungi_max_df <- aggregate(Fungi_fin_df['Value'], by = list(Group = Fungi_fin_df$Group), FUN = max)
Fungi_mean_df <- aggregate(Fungi_fin_df['Value'], by = list(Group = Fungi_fin_df$Group), FUN = mean)

n <- nrow(Fungi_mean_df)
Fungi_dunn_res <- dunn_test(Fungi_fin_df, Value ~ Group, p.adjust.method = 'fdr')
Fungi_dunn_res_df <- data.frame(Fungi_dunn_res[-1])

Fungi_pvalue_df <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
  for(j in (i + 1):n){ 
    k <- k + 1
    Fungi_pvalue_df[i,j] <- Fungi_dunn_res_df$p.adj[k]
    Fungi_pvalue_df[j,i] <- Fungi_dunn_res_df$p.adj[k]
  }
}

Fungi_letter_df <- orderPvalue(Fungi_mean_df$Group, Fungi_mean_df$Value, 0.05, Fungi_pvalue_df, console = TRUE)
Fungi_fin_df$Group <- as.factor(Fungi_fin_df$Group)
Fungi_letter_df <- Fungi_letter_df[levels(Fungi_fin_df$Group),]
Fungi_letter_vector <- Fungi_letter_df$groups

color_mapping <- color_manual2

Fungi_fin_df <- Fungi_fin_df %>%
  mutate(
    point_color = ifelse(
      Roles == "Peripherals",
      "grey80",
      color_mapping[as.character(Group)]
    )
  )

p5 <- ggplot(Fungi_fin_df, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color, group = Group),
    position = position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8, size = 1
  ) +
  scale_color_identity() + 
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
  geom_text(data = Fungi_max_df, aes(y = Value * 1.2, label = Fungi_letter_vector), 
            position = position_dodge(0.9), size = 7 / 2.835) + 
  labs(
    title = '',
    x = '',
    y = 'Fungi',
  ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )


Protist_fin_df <- result_list[["Protist"]] %>%
  mutate(
    level2 = case_when(
      level2 == "Cellular community - prokaryotes" ~ "Cellular community",
      TRUE ~ level2
    ),
    level2 = replace_na(level2, "Others")
  )

bacteria_level2 <- result_list[["Bacteria"]] %>%
  mutate(
    level2 = case_when(
      level2 == "Cellular community - prokaryotes" ~ "Cellular community",
      TRUE ~ level2
    ),
    level2 = replace_na(level2, "Others")
  ) %>%
  distinct(level2) %>%
  pull(level2)

protist_level2 <- Protist_fin_df %>% 
  distinct(level2) %>% 
  pull(level2)

missing_level2 <- setdiff(bacteria_level2, protist_level2)
supplement_df <- data.frame(
  level2 = missing_level2,
  Degree = 0,
  roles = "Peripherals"
)

Protist_fin_df <- Protist_fin_df %>%
  bind_rows(supplement_df) %>%
  rename(Group = level2, Value = Degree, Roles = roles)
Protist_fin_df <- Protist_fin_df %>% filter(Group != "Others")

names(Protist_fin_df) <- c('Group', 'Value' ,'Roles')
#Protist_fin_df$Group <- factor(Protist_fin_df$Group,levels = c(sort(unique(Protist_fin_df$Group[Protist_fin_df$Group != "Others"])), "Others"))
Protist_fin_df$Group <- factor(Protist_fin_df$Group,levels = sort(unique(Protist_fin_df$Group)))
Protist_max_df <- aggregate(Protist_fin_df['Value'], by = list(Group = Protist_fin_df$Group), FUN = max)
Protist_mean_df <- aggregate(Protist_fin_df['Value'], by = list(Group = Protist_fin_df$Group), FUN = mean)

n <- nrow(Protist_mean_df)
Protist_dunn_res <- dunn_test(Protist_fin_df, Value ~ Group, p.adjust.method = 'fdr')
Protist_dunn_res_df <- data.frame(Protist_dunn_res[-1])

Protist_pvalue_df <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
  for(j in (i + 1):n){ 
    k <- k + 1
    Protist_pvalue_df[i,j] <- Protist_dunn_res_df$p.adj[k]
    Protist_pvalue_df[j,i] <- Protist_dunn_res_df$p.adj[k]
  }
}

Protist_letter_df <- orderPvalue(Protist_mean_df$Group, Protist_mean_df$Value, 0.05, Protist_pvalue_df, console = TRUE)
Protist_fin_df$Group <- as.factor(Protist_fin_df$Group)
Protist_letter_df <- Protist_letter_df[levels(Protist_fin_df$Group),]
Protist_letter_vector <- Protist_letter_df$groups

color_mapping <- color_manual2

Protist_fin_df <- Protist_fin_df %>%
  mutate(
    point_color = ifelse(
      Roles == "Peripherals",
      "grey80",
      color_mapping[as.character(Group)]
    )
  )
Protist_fin_df <- Protist_fin_df %>% filter(Group != "Others")
p6 <- ggplot(Protist_fin_df, aes(x = Group, y = Value)) + 
  geom_point(
    aes(color = point_color, group = Group),
    position = position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8, size = 1
  ) +
  scale_color_identity() + 
  geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
  geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
  geom_text(data = Protist_max_df, aes(y = Value * 1.2, label = Protist_letter_vector), 
            position = position_dodge(width=0.9), size = 7 / 2.835) + 
  labs(
    title = '',
    x = '',
    y = 'Protist',
  ) + 
    theme_bw() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right"
    )
final_plot <- wrap_plots(p1,p4,p2,p5,p3,p6 ,ncol = 2) + plot_annotation(tag_levels = "a")


width <- 17.5
height <- 20
name <- "./FigureS14"
ggsave(paste0(name, ".jpg"), final_plot, width = width, height = height, units = "cm",dpi = 900,bg = "white")
ggsave(paste0(name, ".pdf"), final_plot, width = width, height = height, units = "cm",dpi = 900,bg = "white")
ggsave(paste0(name, ".tiff"), final_plot, width = width, height = height, units = "cm",dpi = 900,bg = "white")