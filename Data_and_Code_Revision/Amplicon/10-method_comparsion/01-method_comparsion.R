### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "01-method_comparsion"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}


# Import package
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(Hmisc)
library(corrplot)

### Define variable -----------------------------------------------------------
amplicon<- c("16S", "ITS", "Protist")
r_threshold <- 0.3
color_manual <- c(rep(c("#7fc97f", "#beaed4", "#fdc086"), each = 10))
# -----------------------------------------------------------------------------


### Get results ---------------------------------------------------------------
all_top_order_df <- all_relative_df <- all_absolute_df <- data.frame()
for (amp in amplicon) {
    top_order_df <- read.csv(paste0("../06-community_composition/02-community_composition_Order/", amp, "_top10Order_color_meanValue.csv"), header = T)
    top_order_df <- top_order_df[1:10,]
    
    absolute_df <- read.csv(paste0("../01-sort_data/04-barplot_data/", amp, "_barplot_data_absolute_sorted.csv"), header = T)
    absolute_df <- absolute_df %>%
        filter(Level == "Order", Taxonomy %in%  top_order_df$Taxonomy) %>%
        tibble::column_to_rownames(var = "Taxonomy") %>%
        select(-Mean, -Level)
    
    relative_df <- read.csv(paste0("../01-sort_data/04-barplot_data/", amp, "_barplot_data_relative_sorted.csv"), header = T)
    relative_df <- relative_df %>%
        filter(Level == "Order", Taxonomy %in%  top_order_df$Taxonomy) %>%
        tibble::column_to_rownames(var = "Taxonomy") %>%
        select(-Mean, -Level)
    top_order_df$Clade <- ifelse(amp == "16S", "Bacteria", ifelse(amp == "ITS", "Fungi", "Protist"))
    all_top_order_df <- bind_rows(all_top_order_df,top_order_df )
    all_relative_df <- bind_rows(all_relative_df, relative_df)
    all_absolute_df <- bind_rows(all_absolute_df, absolute_df)
}

sub_relative_df <- data.frame(t(all_relative_df[,!is.na(colSums(all_relative_df))]))
sub_relative_df <- sub_relative_df[all_top_order_df$Taxonomy]

sub_absolute_df <- data.frame(t(all_absolute_df[,!is.na(colSums(all_absolute_df))]))
sub_absolute_df <- sub_absolute_df[all_top_order_df$Taxonomy]

cor_relative <- rcorr(as.matrix(sub_relative_df), type = 'spearman')
cor_absolute <- rcorr(as.matrix(sub_absolute_df), type = 'spearman')

relative_r <- cor_relative$r
relative_r[upper.tri(relative_r)] <- NA
relative_r_df <- melt(relative_r, na.rm = T, value.name = 'R')

relative_p <- cor_relative$P
relative_p[upper.tri(relative_p)] <- NA
diag(relative_p) <- 0
relative_p_df <- melt(relative_p, na.rm = T, value.name = 'P')

relative_fin_df <- merge(relative_r_df, relative_p_df, by = c('Var1', 'Var2'))
relative_fin_df$Method <- ifelse(relative_fin_df$Var1 == relative_fin_df$Var2, 'Diag', 'Relative')

absolute_r <- cor_absolute$r
absolute_r[lower.tri(absolute_r)] <- NA
diag(absolute_r) <- NA
absolute_r_df <- melt(absolute_r, na.rm = T, value.name = 'R')

absolute_p <- cor_absolute$P
absolute_p[lower.tri(absolute_p)] <- NA
absolute_p_df <- melt(absolute_p, na.rm = T, value.name = 'P')

absolute_fin_df <- merge(absolute_r_df, absolute_p_df, by = c('Var1', 'Var2'))
absolute_fin_df$Method <- 'Absolute'


fin_df <- rbind(relative_fin_df, absolute_fin_df)
fin_df$R <- ifelse(abs(fin_df$R) < r_threshold, 0, fin_df$R)
fin_df$R <- ifelse(fin_df$P > 0.001, 0, fin_df$R)
write_df <- subset(fin_df, abs(R) > r_threshold & P < 0.001 & Method != 'Diag')
write.csv(write_df, paste0(dir_name, '/correlation.csv'), row.names = F)

fin_matrix <- dcast(fin_df, Var1 ~ Var2, value.var = 'R')
row.names(fin_matrix) <- fin_matrix$Var1
fin_matrix <- as.matrix(fin_matrix[-1])

width <- 16
height <- 16

pdf(paste0(dir_name, '/corrplot.pdf'), width = width / 2.54, height = height / 2.54)
corrplot(fin_matrix, diag = F, tl.col = 'black', cl.cex = 0.5, tl.cex = 0.5, cl.pos = 'b', cl.length = 5, cl.ratio = 0.2, tl.srt = 45)
dev.off()


# Relative
tmp_relative_data <- sub_relative_df %>%
    pivot_longer(
        cols = everything(),
        names_to = "Group",
        values_to = "Value"
    )
tmp_relative_data$Group <- factor(tmp_relative_data$Group, levels = rev(all_top_order_df$Taxonomy))
tmp_relative_data$Value <- log10(tmp_relative_data$Value + 1E-12)

summarise_df <- tmp_relative_data %>%
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

p <- ggplot(data = tmp_relative_data, mapping = aes(x = Group, y = Value, color = Group)) +
    
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.5, stroke = 0) +
    geom_boxplot(width = 0.68, na.rm = TRUE, outlier.alpha = 0) +
    
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Relative abundance"
    ) + 
    scale_color_manual(values = rev(color_manual)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(-4, 0))+
    theme_minimal() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
    )


width <- 12.7
height <- 3
name <- paste0(dir_name, "/abuncance_relative")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


# Absolute
tmp_absolute_data <- sub_absolute_df %>%
    pivot_longer(
        cols = everything(),
        names_to = "Group",
        values_to = "Value"
    )
tmp_absolute_data$Group <- factor(tmp_absolute_data$Group, levels = all_top_order_df$Taxonomy)

tmp_absolute_data$Value <- log10(tmp_absolute_data$Value + 1)
summarise_df <- tmp_absolute_data %>%
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

p <- ggplot(data = tmp_absolute_data, mapping = aes(x = Group, y = Value, color = Group)) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.5, stroke = 0) +
    geom_boxplot(width = 0.68, na.rm = TRUE, outlier.alpha = 0) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Absolute abundance"
    ) + 
    scale_color_manual(values = color_manual) +
    scale_y_continuous(expand = c(0, 0), breaks = c(5, 7, 9, 11)) +
    coord_cartesian(ylim = c(5, 11))+
    theme_minimal() + theme(
        text = element_text(color = "black", size = 6),
        plot.title = element_text(size = 7, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
        panel.spacing = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.key.size = unit(0.25, "cm"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
    )

width <- 12.7
height <- 3
name <- paste0(dir_name, "/abuncance_absolute")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")



# Sort Data
names(relative_fin_df)[1:4] <- c("Taxa1", "Taxa2", "RPM_R", "RMP_P")
names(absolute_fin_df)[1:4] <- c("Taxa2", "Taxa1", "QPM_R", "QMP_P")
sorted_df <- merge(relative_fin_df[1:4], absolute_fin_df[1:4])
sorted_df <- merge(sorted_df, all_top_order_df[c("Taxonomy", "Clade")], by.x = "Taxa1", by.y = "Taxonomy")
names(sorted_df)[ncol(sorted_df)] <- "Taxa1_clade"
sorted_df <- merge(sorted_df, all_top_order_df[c("Taxonomy", "Clade")], by.x = "Taxa2", by.y = "Taxonomy")
names(sorted_df)[ncol(sorted_df)] <- "Taxa2_clade"

sorted_df <- sorted_df[c(2, 1, 7:8, 3:6)]
name <- paste0(dir_name, "/detail_correlation")

write.csv(sorted_df, paste0(name, ".csv"), quote = F, row.names = F)
# -----------------------------------------------------------------------------
