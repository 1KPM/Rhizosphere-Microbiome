### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "11-closeness_and_betweenness"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(patchwork)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
feature <- c('Betweenness', 'Closeness')
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
p_list <- NULL
for (fea in feature) {
    all_data_df <- all_summarise_df<- all_pvalue_df <- data.frame()
    
    
    for (typ in type) {
        first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
        second_prefix <- ifelse(typ == "all", "02", "01")
        
        dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
        
        hub_info <- read.csv(paste0(dir_path, '/', typ, '_network_hub_info.csv'),row.names = 1)
        tax_table <- read.csv(paste0(dir_path, '/', typ, '_tax_table.csv'), row.names = 1)
        
        tmp_df <- merge(hub_info, tax_table, by = 'row.names')
        names(tmp_df)[2:4] <- c('Degree', 'Closeness', 'Betweenness')
        tmp_df <- tmp_df %>%
            mutate(Clade = case_when(Clade == "Protist" ~ "Protists", TRUE ~ Clade),
                   Clade = factor(Clade, levels = c('Bacteria', 'Fungi', 'Protists')))
        
    
        data_df <- tmp_df[c('Clade', fea)]
        names(data_df) <- c('Group', 'Value')
        
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
        
        type_text <- ifelse(typ == "all", "Whole network", ifelse(typ == "inter", "Interkingdom network", "Intrakingdom network"))
        data_df$Type <- paste0(type_text, " (ASVs)")
        summarise_df$Type <- paste0(type_text, " (ASVs)")
        pvalue_df$Type <- paste0(type_text, " (ASVs)")
        
        all_data_df <- bind_rows(all_data_df, data_df)
        all_summarise_df <- bind_rows(all_summarise_df, summarise_df)
        all_pvalue_df <- bind_rows(all_pvalue_df, pvalue_df)
        
        name <- paste0(dir_name, "/", fea)
        write.csv(all_data_df, paste0(name, "_data.csv"), row.names = F)
        write.csv(all_summarise_df, paste0(name, "_summarise.csv"), row.names = F)
        write.csv(all_pvalue_df, paste0(name, "_pvalue.csv"), row.names = F)
    }
    
    all_data_df$Type <- factor(all_data_df$Type, 
                               levels = c("Whole network (ASVs)", "Interkingdom network (ASVs)", "Intrakingdom network (ASVs)"))
    all_summarise_df$Type <- factor(all_summarise_df$Type, 
                                    levels = c("Whole network (ASVs)", "Interkingdom network (ASVs)", "Intrakingdom network (ASVs)"))
    p <- ggplot(all_data_df, aes(x = Group, y = Value)) + 
        geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
                   alpha = 0.4, size = 3, stroke = 0) +
        geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
        geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
        geom_text(
            data = all_summarise_df,
            mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
            position = position_dodge(0.9),
            size = 7 / 2.835
        ) +
        labs(
            x = '',
            y = paste0(fea, ' centrality'),
        ) + 
        facet_grid(
            Type ~ .,
            scales = "free_y",
            switch = NULL,
            space = "fixed",
        ) +
        theme_bw() + 
        scale_color_manual(values = color_manual) + 
        theme(
            text = element_text(color = "black", size = 6),
            plot.title = element_text(size = 7, hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.title = element_text(size = 7),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6, color = "black"),
            # axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
            panel.spacing = unit(0.1, "cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.position = 'none')
    p_list[[fea]] <- p
    width <- 7
    height <- 15
    name <- paste0(dir_name, "/", fea)
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}

p <- p_list[["Betweenness"]] | p_list[["Closeness"]]

width <- 14
height <- 15
name <- paste0(dir_name, "/Betweenness_and_Closeness")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------
