### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "03-qpcr_vs_rpm"
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
qpcr_data <- read.csv("../00-rawdata/qpcr/qpcr_results.csv", header = T, row.names = 1)


### Get results ---------------------------------------------------------------
all_top_order_df <- all_relative_df <- all_qpcr_df <- data.frame()
for (amp in amplicon) {
    top_order_df <- read.csv(paste0("../06-community_composition/02-community_composition_Order/", amp, "_top10Order_color_meanValue.csv"), header = T)
    top_order_df <- top_order_df[1:10,]
    
    relative_df <- read.csv(paste0("../01-sort_data/04-barplot_data/", amp, "_barplot_data_relative_sorted.csv"), header = T)
    relative_df <- relative_df %>%
        filter(Level == "Order", Taxonomy %in%  top_order_df$Taxonomy) %>%
        tibble::column_to_rownames(var = "Taxonomy") %>%
        select(-Mean, -Level)
    relative_df <- relative_df[row.names(qpcr_data)]
    
    kingdom_label <- ifelse(amp == "16S", "Bacteria", ifelse(amp == "ITS", "Fungi", "Protists"))
    qpcr_df <- relative_df
    for (i in names(qpcr_df)) {
        qpcr_df[i] <- qpcr_df[i] * qpcr_data[i, kingdom_label] * 200 * 100
    }
    
    top_order_df$Clade <- ifelse(amp == "16S", "Bacteria", ifelse(amp == "ITS", "Fungi", "Protist"))
    all_top_order_df <- bind_rows(all_top_order_df, top_order_df)
    all_relative_df <- bind_rows(all_relative_df, relative_df)
    all_qpcr_df <- bind_rows(all_qpcr_df, qpcr_df)
}

sub_relative_df <- data.frame(t(all_relative_df[,!is.na(colSums(all_relative_df))]))
sub_relative_df <- sub_relative_df[all_top_order_df$Taxonomy]

sub_qpcr_df <- data.frame(t(all_qpcr_df[,!is.na(colSums(all_qpcr_df))]))
sub_qpcr_df <- sub_qpcr_df[all_top_order_df$Taxonomy]

cor_relative <- rcorr(as.matrix(sub_relative_df), type = 'spearman')
cor_qpcr <- rcorr(as.matrix(sub_qpcr_df), type = 'spearman')

relative_r <- cor_relative$r
relative_r[upper.tri(relative_r)] <- NA
relative_r_df <- melt(relative_r, na.rm = T, value.name = 'R')

relative_p <- cor_relative$P
relative_p[upper.tri(relative_p)] <- NA
diag(relative_p) <- 0
relative_p_df <- melt(relative_p, na.rm = T, value.name = 'P')

relative_fin_df <- merge(relative_r_df, relative_p_df, by = c('Var1', 'Var2'))
relative_fin_df$Method <- ifelse(relative_fin_df$Var1 == relative_fin_df$Var2, 'Diag', 'Relative')

qpcr_r <- cor_qpcr$r
qpcr_r[lower.tri(qpcr_r)] <- NA
diag(qpcr_r) <- NA
qpcr_r_df <- melt(qpcr_r, na.rm = T, value.name = 'R')

qpcr_p <- cor_qpcr$P
qpcr_p[lower.tri(qpcr_p)] <- NA
qpcr_p_df <- melt(qpcr_p, na.rm = T, value.name = 'P')

qpcr_fin_df <- merge(qpcr_r_df, qpcr_p_df, by = c('Var1', 'Var2'))
qpcr_fin_df$Method <- 'qpcr'


fin_df <- rbind(relative_fin_df, qpcr_fin_df)
fin_df$R <- ifelse(abs(fin_df$R) < r_threshold, 0, fin_df$R)
fin_df$R <- ifelse(fin_df$P > 0.05, 0, fin_df$R)
write_df <- subset(fin_df, abs(R) > r_threshold & P < 0.05 & Method != 'Diag')
write.csv(write_df, paste0(dir_name, '/correlation.csv'), row.names = F)

fin_matrix <- dcast(fin_df, Var1 ~ Var2, value.var = 'R')
row.names(fin_matrix) <- fin_matrix$Var1
fin_matrix <- as.matrix(fin_matrix[-1])

width <- 16
height <- 16

pdf(paste0(dir_name, '/corrplot.pdf'), width = width / 2.54, height = height / 2.54)
corrplot(fin_matrix, diag = F, tl.col = 'black', cl.cex = 0.5, tl.cex = 0.5, cl.pos = 'b', cl.length = 5, cl.ratio = 0.2, tl.srt = 45)
dev.off()


# Sort Data
names(relative_fin_df)[1:4] <- c("Taxa1", "Taxa2", "RPM_R", "RMP_P")
names(qpcr_fin_df)[1:4] <- c("Taxa2", "Taxa1", "qPCR_R", "QMP_P")
sorted_df <- merge(relative_fin_df[1:4], qpcr_fin_df[1:4])
sorted_df <- merge(sorted_df, all_top_order_df[c("Taxonomy", "Clade")], by.x = "Taxa1", by.y = "Taxonomy")
names(sorted_df)[ncol(sorted_df)] <- "Taxa1_clade"
sorted_df <- merge(sorted_df, all_top_order_df[c("Taxonomy", "Clade")], by.x = "Taxa2", by.y = "Taxonomy")
names(sorted_df)[ncol(sorted_df)] <- "Taxa2_clade"

sorted_df <- sorted_df[c(2, 1, 7:8, 3:6)]
name <- paste0(dir_name, "/detail_correlation")

write.csv(sorted_df, paste0(name, ".csv"), quote = F, row.names = F)
# -----------------------------------------------------------------------------


