### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(1994)

# Import Packages
library(ggplot2)
library(vegan)
library(RColorBrewer)

# Create Directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------

### Import Data ----------------------------------------------------------------
metadata_rs <- read.csv('data/rhizosphere_metadata_merge_info.csv')
metadata_tree <- read.csv('data/tree_metadata_merge_info.csv')
tree_order_color <- read.csv('data/Tree_top10_order_color.csv')
# ------------------------------------------------------------------------------


### Get Results ----------------------------------------------------------------
# 1 Amplicon (16S)
set.seed(1994)
matrix_df <- read.table(paste0("data/diversity/16S/distance-matrix.tsv"), row.names = 1)

group_df <- metadata_rs[metadata_rs$FileID %in% rownames(matrix_df),]
rownames(group_df) <- group_df$FileID

sub_group <- group_df["Order"]

names(sub_group) <- "Group"
tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
tmp_dist <- as.dist(tmp_matrix)

pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
variance <- pcoa_sig['Model', 'R2']
p.val <- pcoa_sig['Model', 'Pr(>F)']

cmdscale_total <- cmdscale(tmp_dist, eig = T)
eig <- cmdscale_total$eig

points <- data.frame(cmdscale_total$points)
data_df <- cbind(points, Group = tmp_group[rownames(points),])

data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, tree_order_color$Order, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

p_text <- ifelse(p.val < 0.001, "p < 0.001", 
                 ifelse(p.val < 0.01, "p < 0.01", 
                        ifelse(p.val < 0.05, "p < 0.05", paste0("p = ", format(p.val, digits = 2)))))
title_text <- "Amplicon 16S"

p_16s <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')


# 2 Amplicon (ITS)
set.seed(1994)
matrix_df <- read.table(paste0("data/diversity/ITS/distance-matrix.tsv"), row.names = 1)

group_df <- metadata_rs[metadata_rs$FileID %in% rownames(matrix_df),]
rownames(group_df) <- group_df$FileID

sub_group <- group_df["Order"]

names(sub_group) <- "Group"
tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
tmp_dist <- as.dist(tmp_matrix)

pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
variance <- pcoa_sig['Model', 'R2']
p.val <- pcoa_sig['Model', 'Pr(>F)']

cmdscale_total <- cmdscale(tmp_dist, eig = T)
eig <- cmdscale_total$eig

points <- data.frame(cmdscale_total$points)
data_df <- cbind(points, Group = tmp_group[rownames(points),])

data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, tree_order_color$Order, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

p_text <- ifelse(p.val < 0.001, "p < 0.001", 
                 ifelse(p.val < 0.01, "p < 0.01", 
                        ifelse(p.val < 0.05, "p < 0.05", paste0("p = ", format(p.val, digits = 2)))))
title_text <- "Amplicon ITS"

p_its <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

# 3 Amplicon (Protist)
set.seed(1994)
matrix_df <- read.table(paste0("data/diversity/Protist/distance-matrix.tsv"), row.names = 1)

group_df <- metadata_rs[metadata_rs$FileID %in% rownames(matrix_df),]
rownames(group_df) <- group_df$FileID

sub_group <- group_df["Order"]

names(sub_group) <- "Group"
tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
tmp_dist <- as.dist(tmp_matrix)

pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
variance <- pcoa_sig['Model', 'R2']
p.val <- pcoa_sig['Model', 'Pr(>F)']

cmdscale_total <- cmdscale(tmp_dist, eig = T)
eig <- cmdscale_total$eig

points <- data.frame(cmdscale_total$points)
data_df <- cbind(points, Group = tmp_group[rownames(points),])

data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, tree_order_color$Order, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

p_text <- ifelse(p.val < 0.001, "p < 0.001", 
                 ifelse(p.val < 0.01, "p < 0.01", 
                        ifelse(p.val < 0.05, "p < 0.05", paste0("p = ", format(p.val, digits = 2)))))
title_text <- "Amplicon 18S"

p_pro <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')


p_amplicon <- cowplot::plot_grid(
    p_16s,
    p_its, 
    p_pro, 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 1, nrow = 3
)


# 4 Metagenome (Bacteria)
set.seed(1994)
matrix_df <- read.csv(paste0("data/diversity/Meta/Bacteria_KO_bray.csv"), row.names = 1)

group_df <- metadata_tree[metadata_tree$TreeID %in% rownames(matrix_df),]
rownames(group_df) <- group_df$TreeID

sub_group <- group_df["Order"]

names(sub_group) <- "Group"
tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
tmp_dist <- as.dist(tmp_matrix)

pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
variance <- pcoa_sig['Model', 'R2']
p.val <- pcoa_sig['Model', 'Pr(>F)']

cmdscale_total <- cmdscale(tmp_dist, eig = T)
eig <- cmdscale_total$eig

points <- data.frame(cmdscale_total$points)
data_df <- cbind(points, Group = tmp_group[rownames(points),])

data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, tree_order_color$Order, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

p_text <- ifelse(p.val < 0.001, "p < 0.001", 
                 ifelse(p.val < 0.01, "p < 0.01", 
                        ifelse(p.val < 0.05, "p < 0.05", paste0("p = ", format(p.val, digits = 2)))))
title_text <- "Metagenome bacteria"

p_meta_bacteria <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

# 5 Metagenome (Fungi)
set.seed(1994)
matrix_df <- read.csv(paste0("data/diversity/Meta/Fungi_KO_bray.csv"), row.names = 1)

group_df <- metadata_tree[metadata_tree$TreeID %in% rownames(matrix_df),]
rownames(group_df) <- group_df$TreeID

sub_group <- group_df["Order"]

names(sub_group) <- "Group"
tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
tmp_dist <- as.dist(tmp_matrix)

pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
variance <- pcoa_sig['Model', 'R2']
p.val <- pcoa_sig['Model', 'Pr(>F)']

cmdscale_total <- cmdscale(tmp_dist, eig = T)
eig <- cmdscale_total$eig

points <- data.frame(cmdscale_total$points)
data_df <- cbind(points, Group = tmp_group[rownames(points),])

data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, tree_order_color$Order, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

p_text <- ifelse(p.val < 0.001, "p < 0.001", 
                 ifelse(p.val < 0.01, "p < 0.01", 
                        ifelse(p.val < 0.05, "p < 0.05", paste0("p = ", format(p.val, digits = 2)))))
title_text <- "Metagenome fungi"

p_meta_fungi <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

# 6 Metagenome (Protist)
set.seed(1994)
matrix_df <- read.csv(paste0("data/diversity/Meta/Protist_KO_bray.csv"), row.names = 1)

group_df <- metadata_tree[metadata_tree$TreeID %in% rownames(matrix_df),]
rownames(group_df) <- group_df$TreeID

sub_group <- group_df["Order"]

names(sub_group) <- "Group"
tmp_group <- sub_group[sub_group$Group != "" & !(is.na(sub_group$Group)), , drop = F]
tmp_matrix <- matrix_df[rownames(tmp_group), rownames(tmp_group),]
tmp_dist <- as.dist(tmp_matrix)

pcoa_sig <- adonis2(tmp_dist ~ Group, data = tmp_group, permutations = 1000)
variance <- pcoa_sig['Model', 'R2']
p.val <- pcoa_sig['Model', 'Pr(>F)']

cmdscale_total <- cmdscale(tmp_dist, eig = T)
eig <- cmdscale_total$eig

points <- data.frame(cmdscale_total$points)
data_df <- cbind(points, Group = tmp_group[rownames(points),])

data_df$Group <- ifelse(data_df$Group %in% tree_order_color$Order, tree_order_color$Order, "Others")
data_df$Group <- factor(data_df$Group, levels = tree_order_color$Order)

p_text <- ifelse(p.val < 0.001, "p < 0.001", 
                 ifelse(p.val < 0.01, "p < 0.01", 
                        ifelse(p.val < 0.05, "p < 0.05", paste0("p = ", format(p.val, digits = 2)))))
title_text <- "Metagenome protist"

p_meta_protist <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.25, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

p_meta <- cowplot::plot_grid(
    p_meta_bacteria,
    p_meta_fungi, 
    p_meta_protist, 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 1, nrow = 3
)

p_legend <- ggplot(data_df, aes(x = X1, y = X2, color = Group)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    labs(x = paste('PCo 1 (', format(100*eig[1]/sum(eig), digits=4), '%)', sep = ''),
         y = paste('PCo 2 (', format(100*eig[2]/sum(eig), digits=4), '%)', sep = '')) + 
    ggtitle(paste(format(100*variance, digits = 3), '% of variance; ', 
                  p_text, ' (', title_text, ')', sep = '')) +
    scale_color_manual(values = tree_order_color$Color, name = "Order") + 
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.key.size = unit(0.5, 'cm'),
          legend.box.spacing = unit(0.1,"cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'right')

legend_p <- cowplot::get_legend(p_legend)

p <- cowplot::plot_grid(
    p_amplicon,
    p_meta, 
    legend_p,
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 3, rel_widths = c(3,3,1)
)
width <- 17
height <- 15
name <- paste0(dir_name, "/FigureS18")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------