### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set Seeds
set.seed(1994)

# Import Packages
library(rstatix)
library(agricolae)
library(ggplot2)
library(RColorBrewer)

# Create Directory
dir_name <- "results"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}
# ------------------------------------------------------------------------------


### Define Variables ------------------------------------------------------------
kingdom <- c('Bacteria', 'Fungi', 'Protist')
feature <- c('Degree', 'Closeness', 'Betweenness')
num <- 10
# ------------------------------------------------------------------------------


### Import Data ----------------------------------------------------------------
tax_all <- read.csv('data/All_core_ASV_taxonomy.csv', row.names = 1)
top_order <- read.csv('data/top_order.csv')

# ------------------------------------------------------------------------------

### Get Results ----------------------------------------------------------------
# 1. within
hub_info <- read.csv("data/Within_Microbiome_network_hub_info.csv", row.names = 1)
raw_df <- merge(hub_info, tax_all, by = 'row.names')
names(raw_df)[2:4] <- c('Degree', 'Closeness', 'Betweenness')
all_data_df <- all_max_df <- all_sig_df <- data.frame()
for (kin in kingdom) {
    f <- "Degree"
    
    res_df <- raw_df[raw_df$Clade == kin, c('Order', f, "roles")]
    names(res_df) <- c('Group', 'Value', 'Role')
    fin_df <- res_df[res_df$Group %in% top_order[,kin],]
    
    ### 删除小于3个点的数据
    order_number <- table(fin_df$Group)
    delete_list <- names(order_number)[order_number < 3] 
    fin_df <- fin_df[!(fin_df$Group %in% delete_list),]
    
    color_df <- data.frame(
        Taxonomy = top_order[, kin],
        Color = colorRampPalette(brewer.pal(9, 'Set1'))(num+1)[1:num]
    )
    color_df <- color_df[color_df$Taxonomy %in% fin_df$Group,]
    
    if (kin == "Protist") {
        fin_df$Group[fin_df$Group == "ATCC50593-Flamella-WIM80-lineage"] <- "ATCC50593-Flamella-WIM80"
        color_df$Taxonomy[color_df$Taxonomy == "ATCC50593-Flamella-WIM80-lineage"] <- "ATCC50593-Flamella-WIM80"
    }
    
    
    fin_df$Group <- factor(fin_df$Group, levels = color_df$Taxonomy)
    
    max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
    max_df$GroupMax <- max(max_df$Value)
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
    letter_df <- letter_df[levels(fin_df$Group),]
    max_df$Label <- letter_df$groups
    
    fin_df$Clade <- kin
    max_df$Clade <- kin
    dunn_res_df$Clade <- kin
    
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, max_df)
    all_sig_df <- rbind(all_sig_df, dunn_res_df)
}

all_hub_data_df <- all_data_df[all_data_df$Role != "Peripherals",]
all_nonhub_data_df <- all_data_df[all_data_df$Role == "Peripherals",]

p1 <- ggplot() + 
    geom_point(data = all_nonhub_data_df, aes(x = Group, y = Value, fill = "Group"), 
               position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 2, stroke = 0, color = "#999999") +
    geom_point(data = all_hub_data_df, aes(x = Group, y = Value, fill = "Group"), 
               position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.8, size = 2, stroke = 0, color = "red") +
    
    geom_boxplot(data = all_data_df, aes(x = Group, y = Value), width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(data = all_data_df, aes(x = Group, y = Value), width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = all_max_df, aes(x = Group, y = Value + GroupMax * 0.2, label = Label), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(
        x = '',
        y = paste0(f, ' centrality of within-kingdom networks'),
    ) + 
    facet_wrap(vars(Clade), scales = "free", strip.position = "right") +
    theme_bw() + 
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
name <- paste0(dir_name, "/FigureS11ABC")
write.csv(all_sig_df, paste0(name, '_sig_fino.csv'), row.names = F)   
write.csv(all_data_df, paste0(name, '_data.csv'), row.names = F)
write.csv(all_max_df, paste0(name, '_data_mean.csv'), row.names = F)

# 2. cross
hub_info <- read.csv("data/Cross_Microbiome_network_hub_info.csv", row.names = 1)
raw_df <- merge(hub_info, tax_all, by = 'row.names')
names(raw_df)[2:4] <- c('Degree', 'Closeness', 'Betweenness')
all_data_df <- all_max_df <- all_sig_df <- data.frame()
for (kin in kingdom) {
    f <- "Degree"
    
    res_df <- raw_df[raw_df$Clade == kin, c('Order', f, "roles")]
    names(res_df) <- c('Group', 'Value', 'Role')
    fin_df <- res_df[res_df$Group %in% top_order[,kin],]
    
    ### 删除小于3个点的数据
    order_number <- table(fin_df$Group)
    delete_list <- names(order_number)[order_number < 3] 
    fin_df <- fin_df[!(fin_df$Group %in% delete_list),]
    
    color_df <- data.frame(
        Taxonomy = top_order[, kin],
        Color = colorRampPalette(brewer.pal(9, 'Set1'))(num+1)[1:num]
    )
    color_df <- color_df[color_df$Taxonomy %in% fin_df$Group,]
    
    if (kin == "Protist") {
        fin_df$Group[fin_df$Group == "ATCC50593-Flamella-WIM80-lineage"] <- "ATCC50593-Flamella-WIM80"
        color_df$Taxonomy[color_df$Taxonomy == "ATCC50593-Flamella-WIM80-lineage"] <- "ATCC50593-Flamella-WIM80"
    }
    

    fin_df$Group <- factor(fin_df$Group, levels = color_df$Taxonomy)
    
    max_df <- aggregate(fin_df['Value'], by = list(Group = fin_df$Group), FUN = max)
    max_df$GroupMax <- max(max_df$Value)
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
    letter_df <- letter_df[levels(fin_df$Group),]
    max_df$Label <- letter_df$groups
    
    fin_df$Clade <- kin
    max_df$Clade <- kin
    dunn_res_df$Clade <- kin
    
    all_data_df <- rbind(all_data_df, fin_df)
    all_max_df <- rbind(all_max_df, max_df)
    all_sig_df <- rbind(all_sig_df, dunn_res_df)
}

all_hub_data_df <- all_data_df[all_data_df$Role != "Peripherals",]
all_nonhub_data_df <- all_data_df[all_data_df$Role == "Peripherals",]

p2 <- ggplot() + 
    geom_point(data = all_nonhub_data_df, aes(x = Group, y = Value, fill = "Group"), 
               position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 2, stroke = 0, color = "#999999") +
    geom_point(data = all_hub_data_df, aes(x = Group, y = Value, fill = "Group"), 
               position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.8, size = 2, stroke = 0, color = "red") +
    
    geom_boxplot(data = all_data_df, aes(x = Group, y = Value), width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(data = all_data_df, aes(x = Group, y = Value), width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(data = all_max_df, aes(x = Group, y = Value + GroupMax * 0.2, label = Label), 
              position = position_dodge(0.9), size = 2.5) + 
    labs(
        x = '',
        y = paste0(f, ' centrality of cross-kingdom networks'),
    ) + 
    facet_wrap(vars(Clade), scales = "free", strip.position = "right") +
    theme_bw() + 
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
name <- paste0(dir_name, "/FigureS11DEF")
write.csv(all_sig_df, paste0(name, '_sig_fino.csv'), row.names = F)   
write.csv(all_data_df, paste0(name, '_data.csv'), row.names = F)
write.csv(all_max_df, paste0(name, '_data_mean.csv'), row.names = F)

p <- cowplot::plot_grid(
    p1 + theme(plot.margin = margin(0, 0, 0, 0, "cm")),
    p2 + theme(plot.margin = margin(0, 0, 0, 0, "cm")), 
    align = "hv", axis = "tblr", hjust = 0, vjust = 0,
    ncol = 1, nrow = 2
)

width <- 17
height <- 16
name <- paste0(dir_name, "/FigureS11")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


# ------------------------------------------------------------------------------



