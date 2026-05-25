### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "10-centrality"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

typ <- "inter"
centrality <- c('Degree', 'Closeness', 'Betweenness')
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
all_data_df <- all_summarise_df<- all_pvalue_df <- data.frame()
for (ord in names(top10order_list)) {
    dir_path <- paste0("03-get_top10order_network_property")
    
    hub_info <- read.csv(paste0(dir_path, '/', ord, "_", typ, '_network_hub_info.csv'), row.names = 1)
    tax_table <- read.csv(paste0(dir_path, '/', ord, "_",  typ, '_tax_table.csv'), row.names = 1)
    
    tmp_df <- merge(hub_info, tax_table, by = 'row.names')
    names(tmp_df)[2:4] <- centrality
    tmp_df <- tmp_df %>%
        mutate(Clade = case_when(Clade == "Protist" ~ "Protists", TRUE ~ Clade),
               Clade = factor(Clade, levels = c('Bacteria', 'Fungi', 'Protists')))
    
    for (cen in centrality) {
        data_df <- tmp_df[c('Clade', cen)]
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
        
        summarise_df <- summarise_df %>% mutate(Order = ord, Centrality = cen)
        data_df <- data_df %>% mutate(Order = ord, Centrality = cen)
        pvalue_df <- pvalue_df %>% mutate(Order = ord, Centrality = cen)
        
        all_data_df <- bind_rows(all_data_df, data_df)
        all_summarise_df <- bind_rows(all_summarise_df, summarise_df)
        all_pvalue_df <- bind_rows(all_pvalue_df, pvalue_df)
    }
}

name <- paste0(dir_name, "/", typ, "_top10order_centrality")
write.csv(all_data_df, paste0(name, "_data.csv"), row.names = F)
write.csv(all_summarise_df, paste0(name, "_summarise.csv"), row.names = F)
write.csv(all_pvalue_df, paste0(name, "_pvalue.csv"), row.names = F)



for (cen in centrality) {
    data_df <- all_data_df %>% filter(Centrality == cen) %>%
        mutate(Order = factor(Order, levels = names(top10order_list)))
    summarise_df <- all_summarise_df %>% filter(Centrality == cen) %>%
        mutate(Order = factor(Order, levels = names(top10order_list)))
    
    p <- ggplot(data_df, aes(x = Group, y = Value)) + 
        geom_point(aes(color = Group), position = position_jitterdodge(dodge.width = 0.6), 
                   alpha = 0.4, size = 3, stroke = 0) +
        geom_boxplot(width = 0.3, alpha = 0.2, na.rm = TRUE) +   
        geom_violin(width = 0.5, alpha = 0.2, na.rm = TRUE) +  
        geom_text(
            data = summarise_df,
            mapping = aes(x = Group, y = max + allmax * 0.3, label = label),
            position = position_dodge(0.9),
            size = 7 / 2.835
        ) +
        labs(
            x = '',
            y = paste0(cen, ' centrality'),
        ) + 
        facet_wrap(~ Order, nrow = 1, ncol = 10) +
        theme_bw() + 
        scale_color_manual(values = color_manual) + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.text = element_text(size = 6, color = 'black'), 
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
              panel.spacing = unit(0.1, "cm"),
              legend.box.spacing = unit(0.1,"cm"),
              legend.key.size = unit(0.25, "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              legend.position = 'none')
    width <- 17.5
    height <- 5
    name <- paste0(dir_name, "/", typ, "_top10order_", cen)
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
    ggsave(paste0(name, ".tiff"), p, width = width, height = height, dpi = 600, units = "cm", compression = "lzw")
}
# ------------------------------------------------------------------------------