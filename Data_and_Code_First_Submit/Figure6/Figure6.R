### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(rstatix)
library(agricolae)
library(RColorBrewer)
library(ggpmisc)

### Import data ----------------------------------------------------------------
metadata_rs <- read.csv("data/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("data/tree_metadata_merge_info.csv")
tree_color_order <- read.csv("data/Tree_top10_order_color.csv")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
tax <- "Order"

## 1. Amplicon
kingdom <- c("16S", "ITS", "Protist")
type <- "Amplicon"
if (type == "Amplicon") {metadata <- metadata_rs} else {metadata <- metadata_tree}
all_data_df <- all_max_df<- all_sig_df <- data.frame()
factor_level <- c()

for (kin in kingdom) {
    raw_df <- read.table(paste0("data/diversity/", kin, "/alpha-diversity.tsv"), row.names = 1)
    names(raw_df) <- "Value"
    
    merge_df <- merge(raw_df, metadata, by.x = "row.names", by.y = "FileID")
    names(merge_df)[1] <- "Sample"
    
    res_df <- merge_df[c(tax, "Value", "Sample")]
    names(res_df)[1] <- "Group"
    
    fin_df <- res_df[res_df$Group != "" & !(is.na(res_df$Group)),]
    if(tax == "Class") {
        tax_color_df <- tree_color_class
        tax_color_df[,1][nrow(tax_color_df)] <- "Other class"
    } else {
        tax_color_df <- tree_color_order
        tax_color_df[,1][nrow(tax_color_df)] <- "Other order"
    }
    
    fin_df$Group <- ifelse(fin_df$Group %in% tax_color_df[,1], fin_df$Group, tax_color_df[,1][nrow(tax_color_df)])
    fin_df$Group <- factor(fin_df$Group, levels = tax_color_df[,1])

    max_df <- fin_df %>%
        group_by(Group) %>%
        summarise(Value = max(Value)) %>%
        ungroup() %>%
        mutate(GroupMax = max(Value))  # 所有组的全局最大值（每行重复）
    
    mean_df <- fin_df %>%
        group_by(Group) %>%
        summarise(Value = mean(Value))  # 每组的平均值
    
    pttest_sig <- pairwise.t.test(fin_df$Value, fin_df$Group, p.adjust.method = "none", paired = FALSE)
    pttest_sig_df <- na.omit(reshape2::melt(pttest_sig$p.value)[c(2, 1, 3)])
    names(pttest_sig_df) <- c("group1", "group2", "p.adj")
    
    n <- length(unique(fin_df$Group))
    pvalue_df <- matrix(1, ncol = n, nrow = n)
    k <- 0
    for(i in 1:(n - 1)) { 
        for(j in (i + 1):n) { 
            k <- k + 1
            pvalue_df[i, j] <- pttest_sig_df$p.adj[k]
            pvalue_df[j, i] <- pttest_sig_df$p.adj[k]
        }
    }
    label_df <- agricolae::orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = TRUE)
    label_df <- label_df[levels(fin_df$Group),]
    max_df$Label <- label_df$groups
    
    kin_text <- ifelse(kin == "Protist", paste0(type, " (18S)"), paste0(type, " (", kin, ")"))
    
    fin_df$Kingdom <- kin_text
    max_df$Kingdom <- kin_text
    pttest_sig_df$Kingdom <- kin_text
    factor_level <- c(factor_level, kin_text)
    
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, max_df)
    all_sig_df <- rbind(all_sig_df, pttest_sig_df)
}        

name <- paste0(dir_name, "/alpha_diversity_", tax, "_", type)
write.csv(all_sig_df, paste0(name, '_sig.csv'), row.names = F)   
write.csv(all_data_df, paste0(name, '_data.csv'), row.names = F)
write.csv(all_max_df, paste0(name, '_mean.csv'), row.names = F)


all_data_df$Kingdom <- factor(all_data_df$Kingdom, levels = factor_level)
all_max_df$Kingdom <- factor(all_max_df$Kingdom, levels = factor_level)

p1 <- ggplot(all_data_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = all_max_df, aes(x = Group, y = Value + GroupMax * 0.1, label = Label),
              position = position_dodge(0.9), size = 7 / 2.835) + 
    labs(
        x = NULL,
        y = paste0("Shannon index"),
    ) + 
    scale_color_manual(values = tax_color_df$Color) +

    facet_grid(
        Kingdom ~ .,
        scales = "free"
    ) +
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = "black", hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5), 
          axis.title = element_text(size = 7, color = "black"), 
          axis.text = element_text(size = 6, color = "black"), 
          legend.title = element_text(size = 7, color = "black"), 
          legend.text = element_text(size = 6,  color = "black"), 
          strip.text = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = "none")


## 2. Metagenome
kingdom <- c("Bacteria", "Fungi", "Protist")
type <- "Metagenome"
if (type == "Amplicon") {metadata <- metadata_rs} else {metadata <- metadata_tree}
all_data_df <- all_max_df<- all_sig_df <- data.frame()
factor_level <- c()

for (kin in kingdom) {
    raw_df <- read.table(paste0("data/diversity/Meta/", kin, "_KO_alpha_diversity.tsv"), row.names = 1)
    names(raw_df) <- "Value"
    
    merge_df <- merge(raw_df, metadata, by.x = "row.names", by.y = "TreeID")
    names(merge_df)[1] <- "Tree"
    
    res_df <- merge_df[c(tax, "Value", "Tree")]
    names(res_df)[1] <- "Group"
    
    fin_df <- res_df[res_df$Group != "" & !(is.na(res_df$Group)),]
    if(tax == "Class") {
        tax_color_df <- tree_color_class
        tax_color_df[,1][nrow(tax_color_df)] <- "Other class"
    } else {
        tax_color_df <- tree_color_order
        tax_color_df[,1][nrow(tax_color_df)] <- "Other order"
    }
    
    fin_df$Group <- ifelse(fin_df$Group %in% tax_color_df[,1], fin_df$Group, tax_color_df[,1][nrow(tax_color_df)])
    fin_df$Group <- factor(fin_df$Group, levels = tax_color_df[,1])
    
    max_df <- fin_df %>%
        group_by(Group) %>%
        summarise(Value = max(Value)) %>%
        ungroup() %>%
        mutate(GroupMax = max(Value))  # 所有组的全局最大值（每行重复）
    
    mean_df <- fin_df %>%
        group_by(Group) %>%
        summarise(Value = mean(Value))  # 每组的平均值
    
    pttest_sig <- pairwise.t.test(fin_df$Value, fin_df$Group, p.adjust.method = "none", paired = FALSE)
    pttest_sig_df <- na.omit(reshape2::melt(pttest_sig$p.value)[c(2, 1, 3)])
    names(pttest_sig_df) <- c("group1", "group2", "p.adj")
    
    n <- length(unique(fin_df$Group))
    pvalue_df <- matrix(1, ncol = n, nrow = n)
    k <- 0
    for(i in 1:(n - 1)) { 
        for(j in (i + 1):n) { 
            k <- k + 1
            pvalue_df[i, j] <- pttest_sig_df$p.adj[k]
            pvalue_df[j, i] <- pttest_sig_df$p.adj[k]
        }
    }
    label_df <- agricolae::orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = TRUE)
    label_df <- label_df[levels(fin_df$Group),]
    max_df$Label <- label_df$groups
    
    kin_text <- paste0(type, " (", kin, ")")
    
    fin_df$Kingdom <- kin_text
    max_df$Kingdom <- kin_text
    pttest_sig_df$Kingdom <- kin_text
    factor_level <- c(factor_level, kin_text)
    
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, max_df)
    all_sig_df <- rbind(all_sig_df, pttest_sig_df)
}        

name <- paste0(dir_name, "/alpha_diversity_", tax, "_", type)
write.csv(all_sig_df, paste0(name, '_sig.csv'), row.names = F)   
write.csv(all_data_df, paste0(name, '_data.csv'), row.names = F)
write.csv(all_max_df, paste0(name, '_mean.csv'), row.names = F)


all_data_df$Kingdom <- factor(all_data_df$Kingdom, levels = factor_level)
all_max_df$Kingdom <- factor(all_max_df$Kingdom, levels = factor_level)

y_lable_height <- c(rep(0.1, 11), rep(0.5, 22))
p2 <- ggplot(all_data_df, aes(x = Group, y = Value)) + 
    geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 3, stroke = 0) +
    geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = all_max_df, aes(x = Group, y = Value + y_lable_height, label = Label),
              position = position_dodge(0.9), size = 7 / 2.835) + 
    labs(
        x = NULL,
        y = NULL,
    ) + 
    scale_color_manual(values = tax_color_df$Color) +
    
    facet_grid(
        Kingdom ~ .,
        scales = "free"
    ) +
    
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = "black", hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5), 
          axis.title = element_text(size = 7, color = "black"), 
          axis.text = element_text(size = 6, color = "black"), 
          legend.title = element_text(size = 7, color = "black"), 
          legend.text = element_text(size = 6,  color = "black"), 
          strip.text = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = "none")


# 1. Amplicon
raw_df <- read.csv("results/alpha_diversity_Order_Amplicon_data.csv")
raw_df$Type <- "Amplicon"
raw_df$Kingdom <- ifelse(grepl("16S", raw_df$Kingdom), "Bacteria", 
                         ifelse(grepl("ITS", raw_df$Kingdom), "Fungi", "Protist"))
res_df <- merge(raw_df, metadata_rs[c("FileID", "TreeID")], by.x = "Sample", by.y = "FileID")

aggregated_df <- res_df %>%
    group_by(Group, Kingdom, Type, TreeID) %>%
    summarise(
        Value = mean(Value),
        .groups = "drop"
    ) %>%
    rename(Tree = TreeID)

# 2. Metagenome
meta_df <- read.csv("results/alpha_diversity_Order_Metagenome_data.csv")
meta_df$Type <- "Metagenome"
meta_df$Kingdom <- ifelse(grepl("Bacteria", meta_df $Kingdom), "Bacteria", 
                          ifelse(grepl("Fungi", meta_df $Kingdom), "Fungi", "Protist"))

# 3. Merge data
tree_list <- intersect(aggregated_df$Tree, meta_df$Tree)

amplicon_df <- aggregated_df %>%
    filter(Tree %in% tree_list) %>%
    select(Group, Kingdom, Tree, Value) %>%
    rename(Amplicon = Value)

metagenome_df <- meta_df %>%
    filter(Tree %in% tree_list) %>%
    select(Group, Kingdom, Tree, Value) %>%
    rename(Metagenome = Value)

data_df <- merge(amplicon_df, metagenome_df, by = c("Group", "Kingdom", "Tree"))
data_df$Group <- factor(data_df$Group, levels = c(tree_color_order$Order, "Other order"))

p3 <- ggplot(data_df, aes(x = Metagenome, y = Amplicon)) +
    geom_point(aes(color = Group), size = 1) +  # 按 Group 着色点
    geom_smooth(method = "lm", formula = y ~ x, 
                color = "#F8766D", linewidth = 0.7) +
    labs(
        x = "Shannon index of metagenomic KO functions",
        y = "Shannon index of microbial species"
    ) +
    facet_wrap(~ Kingdom, ncol = 1, scales = "free", strip.position = "right")+
    scale_color_manual(
        values = tree_color_order$Color, name = "Plant",
        guide = guide_legend(nrow = 2)
    ) +
    guides(color=guide_legend(nrow= 3 , byrow= TRUE )) +
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = "black", hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = "black", hjust = 0.5), 
          axis.title = element_text(size = 7, color = "black"), 
          axis.text = element_text(size = 6, color = "black"), 
          strip.text = element_text(size = 6, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
          legend.margin = margin(-0.1, 0, 0, 0, "cm"),
          legend.title = element_text(size = 7, color = "black"), 
          legend.text = element_text(size = 6, color = "black"), 
          legend.key.size = unit(0.25, "cm"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"
    ) +
    stat_poly_eq(
        aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~")),
        parse = TRUE,
        label.x = 0.95,
        label.y = 0.95,
        size = 7 / 2.835,
        p.digits = 3,
        small.p = TRUE
    )


p <- cowplot::plot_grid(
    p1,
    p2, 
    p3,
    ncol = 3, rel_widths = c(1, 1, 1) 
)
width <- 17.5
height <- 15
name <- paste0(dir_name, "/Figure6")

ggsave(paste0(dir_name, "/Figure6.png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(dir_name, "/Figure6.pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------


