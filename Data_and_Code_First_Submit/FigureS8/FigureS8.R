### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
amplicon <- c("16S", "ITS", "Protist")
color_manual <- c("#7fc97f", "#beaed4", "#fdc086")

# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
tree_order_color <- read.csv("data/Tree_top10_order_color.csv")
metadata_rs <- read.csv("data/rhizosphere_metadata_merge_info.csv")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
all_df <- data.frame()
for (amp in amplicon) {
    raw_df <- read.csv(paste0("data/", amp, "_sample_absolute_abundance.csv"))
    names(raw_df) <- c("Sample", "Value")
    
    raw_df$Group <- ifelse(amp == "Protist", "18S", amp)
    all_df <- bind_rows(all_df, raw_df)
}

name <- paste0(dir_name, "/FigureS8")
write.csv(all_df, paste0(name, ".csv"), quote = F, row.names = F)

# 1. all
data_df <- all_df
data_df$Value <- log10(data_df$Value + 1)
data_df$Group <- factor(data_df$Group, levels = c("16S", "ITS", "18S"))
summarise_df <- data_df %>%
    group_by(Group) %>%
    summarise(
        n = n(),
        mean = mean(Value, na.rm = TRUE),
        min = min(Value, na.rm = TRUE),
        max = max(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE),
        se = sd / sqrt(n),
        ci = qt(0.975, df = n - 1) * se
    ) %>%
    mutate(allmax = max(max))

# Non-parametric test (Kruskal-Wallis rank-sum test with Dunn's post hoc test)
kruskal_test <- kruskal.test(Value ~ Group, data = data_df)
sig_value <- kruskal_test$p.value
dunn_test <- rstatix::dunn_test(data_df, Value ~ Group, p.adjust.method = p_adjust_method)
pvalue_df <- data.frame(KruskalWallis = sig_value, dunn_test[c("group1", "group2", "p.adj")])
pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                        ifelse(pvalue_df$p.adj < 0.01, "**", 
                                ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
n <- nrow(summarise_df)
pvalue_matrix <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
    for(j in (i + 1):n) { 
        k <- k + 1
        pvalue_matrix[i, j] <- pvalue_df$p.adj[k]
        pvalue_matrix[j, i] <- pvalue_df$p.adj[k]
    }
}
letter_df <- agricolae::orderPvalue(summarise_df$Group, summarise_df$mean, 0.05, pvalue_matrix, console = TRUE)
letter_df <- letter_df[levels(data_df$Group),]
summarise_df$label <- letter_df$groups

p_all <- ggplot(data = data_df, mapping = aes(x = Group, y = Value)) +
    geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
    geom_text(
        data = summarise_df,
        mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
        position = position_dodge(0.9),
        size = 7 / 2.835
    ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Absolute abundance\n(gene copies / g soil)"
    ) + 
    scale_color_manual(values = color_manual) +
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
write.csv(data_df, paste0(name, "_data.csv"), row.names = F)
write.csv(summarise_df, paste0(name, "_summarise.csv"), row.names = F)
write.csv(pvalue_df, paste0(name, "_pvalue.csv"), row.names = F)

# 2. 16S
amp <- "16S"
data_df <- subset(all_df, Group == amp)[1:2]
data_df <- merge(data_df, metadata_rs[c("FileID", "Order")], by.x = "Sample", by.y = "FileID")
names(data_df)[3] <- "Group"
data_df$Value <- log10(data_df$Value + 1)
data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, data_df$Group, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

summarise_df <- data_df %>%
    group_by(Group) %>%
    summarise(
        n = n(),
        mean = mean(Value, na.rm = TRUE),
        min = min(Value, na.rm = TRUE),
        max = max(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE),
        se = sd / sqrt(n),
        ci = qt(0.975, df = n - 1) * se
    ) %>%
    mutate(allmax = max(max))

# Non-parametric test (Kruskal-Wallis rank-sum test with Dunn's post hoc test)
kruskal_test <- kruskal.test(Value ~ Group, data = data_df)
sig_value <- kruskal_test$p.value
dunn_test <- rstatix::dunn_test(data_df, Value ~ Group, p.adjust.method = p_adjust_method)
pvalue_df <- data.frame(KruskalWallis = sig_value, dunn_test[c("group1", "group2", "p.adj")])
pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                          ifelse(pvalue_df$p.adj < 0.01, "**", 
                                 ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
n <- nrow(summarise_df)
pvalue_matrix <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
    for(j in (i + 1):n) { 
        k <- k + 1
        pvalue_matrix[i, j] <- pvalue_df$p.adj[k]
        pvalue_matrix[j, i] <- pvalue_df$p.adj[k]
    }
}
letter_df <- agricolae::orderPvalue(summarise_df$Group, summarise_df$mean, 0.05, pvalue_matrix, console = TRUE)
letter_df <- letter_df[levels(data_df$Group),]
summarise_df$label <- letter_df$groups

p_16s <- ggplot(data = data_df, mapping = aes(x = Group, y = Value)) +
    geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
    # geom_text(
    #     data = summarise_df,
    #     mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
    #     position = position_dodge(0.9),
    #     size = 7 / 2.835
    # ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = paste0("Absolute abundance (", amp, ")\n(gene copies / g soil)")
    ) + 
    scale_color_manual(values = tree_order_color$Color) +
    # scale_y_continuous(limits = c(8, 12)) +
    annotate("text", 
             x = -Inf, y = Inf,
             label = "n.s.",
             hjust = -0.5,
             vjust = 1.5,
             size = 7 / 2.835,
             color = "black") + 
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

# 3. its
amp <- "ITS"
data_df <- subset(all_df, Group == amp)[1:2]
data_df <- merge(data_df, metadata_rs[c("FileID", "Order")], by.x = "Sample", by.y = "FileID")
names(data_df)[3] <- "Group"
data_df$Value <- log10(data_df$Value + 1)
data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, data_df$Group, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

summarise_df <- data_df %>%
    group_by(Group) %>%
    summarise(
        n = n(),
        mean = mean(Value, na.rm = TRUE),
        min = min(Value, na.rm = TRUE),
        max = max(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE),
        se = sd / sqrt(n),
        ci = qt(0.975, df = n - 1) * se
    ) %>%
    mutate(allmax = max(max))

# Non-parametric test (Kruskal-Wallis rank-sum test with Dunn's post hoc test)
kruskal_test <- kruskal.test(Value ~ Group, data = data_df)
sig_value <- kruskal_test$p.value
dunn_test <- rstatix::dunn_test(data_df, Value ~ Group, p.adjust.method = p_adjust_method)
pvalue_df <- data.frame(KruskalWallis = sig_value, dunn_test[c("group1", "group2", "p.adj")])
pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                          ifelse(pvalue_df$p.adj < 0.01, "**", 
                                 ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
n <- nrow(summarise_df)
pvalue_matrix <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
    for(j in (i + 1):n) { 
        k <- k + 1
        pvalue_matrix[i, j] <- pvalue_df$p.adj[k]
        pvalue_matrix[j, i] <- pvalue_df$p.adj[k]
    }
}
letter_df <- agricolae::orderPvalue(summarise_df$Group, summarise_df$mean, 0.05, pvalue_matrix, console = TRUE)
letter_df <- letter_df[levels(data_df$Group),]
summarise_df$label <- letter_df$groups

p_its <- ggplot(data = data_df, mapping = aes(x = Group, y = Value)) +
    geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
    # geom_text(
    #     data = summarise_df,
    #     mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
    #     position = position_dodge(0.9),
    #     size = 7 / 2.835
    # ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = paste0("Absolute abundance (", amp, ")\n(gene copies / g soil)")
    ) + 
    scale_color_manual(values = tree_order_color$Color) +
    # scale_y_continuous(limits = c(6, 11)) +
    annotate("text", 
             x = -Inf, y = Inf,
             label = "n.s.",
             hjust = -0.5,
             vjust = 1.5,
             size = 7 / 2.835,
             color = "black") + 
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

# 4. 18S
amp <- "18S"
data_df <- subset(all_df, Group == amp)[1:2]
data_df <- merge(data_df, metadata_rs[c("FileID", "Order")], by.x = "Sample", by.y = "FileID")
names(data_df)[3] <- "Group"
data_df$Value <- log10(data_df$Value + 1)
data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, data_df$Group, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

summarise_df <- data_df %>%
    group_by(Group) %>%
    summarise(
        n = n(),
        mean = mean(Value, na.rm = TRUE),
        min = min(Value, na.rm = TRUE),
        max = max(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE),
        se = sd / sqrt(n),
        ci = qt(0.975, df = n - 1) * se
    ) %>%
    mutate(allmax = max(max))

# Non-parametric test (Kruskal-Wallis rank-sum test with Dunn's post hoc test)
kruskal_test <- kruskal.test(Value ~ Group, data = data_df)
sig_value <- kruskal_test$p.value
dunn_test <- rstatix::dunn_test(data_df, Value ~ Group, p.adjust.method = p_adjust_method)
pvalue_df <- data.frame(KruskalWallis = sig_value, dunn_test[c("group1", "group2", "p.adj")])
pvalue_df$label <- ifelse(pvalue_df$p.adj < 0.001, "***",
                          ifelse(pvalue_df$p.adj < 0.01, "**", 
                                 ifelse(pvalue_df$p.adj < 0.05, "*", "n.s.")))
n <- nrow(summarise_df)
pvalue_matrix <- matrix(1, ncol = n, nrow = n)
k <- 0
for(i in 1:(n - 1)) { 
    for(j in (i + 1):n) { 
        k <- k + 1
        pvalue_matrix[i, j] <- pvalue_df$p.adj[k]
        pvalue_matrix[j, i] <- pvalue_df$p.adj[k]
    }
}
letter_df <- agricolae::orderPvalue(summarise_df$Group, summarise_df$mean, 0.05, pvalue_matrix, console = TRUE)
letter_df <- letter_df[levels(data_df$Group),]
summarise_df$label <- letter_df$groups

p_18s <- ggplot(data = data_df, mapping = aes(x = Group, y = Value)) +
    geom_jitter(mapping = aes(color = Group), width = 0.2, alpha = 0.5, size = 2, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +
    # geom_text(
    #     data = summarise_df,
    #     mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
    #     position = position_dodge(0.9),
    #     size = 7 / 2.835
    # ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = paste0("Absolute abundance (", amp, ")\n(gene copies / g soil)")
    ) + 
    scale_color_manual(values = tree_order_color$Color) +
    # scale_y_continuous(limits = c(6, 11)) +
    annotate("text", 
             x = -Inf, y = Inf,
             label = "n.s.",
             hjust = -0.5,
             vjust = 1.5,
             size = 7 / 2.835,
             color = "black") + 
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

p <- cowplot::plot_grid(
    p_all,
    p_16s + theme(plot.margin = margin(0, 0, 0, 0.5, "cm")),
    p_its, 
    p_18s + theme(plot.margin = margin(0, 0, 0, 0.5, "cm")), 
    labels = c("a", "b", "c", "d"), label_size = 10, label_fontface = "bold",
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 2, nrow = 2
)

width <- 16
height <- 12
name <- paste0(dir_name, "/FigureS8")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------