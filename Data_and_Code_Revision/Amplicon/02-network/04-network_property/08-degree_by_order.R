### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "08-degree_by_order"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(ggh4x)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
kingdom <- c('Bacteria', 'Fungi', 'Protist')
fea <- "Degree"
num <- 10
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
core_taxonomy <- read.csv("../../01-sort_data/02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
top_order <- read.csv("../../01-sort_data/07-top_order/top_order.csv", header = T)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
all_data_df <- all_summarise_df<- all_pvalue_df <- data.frame()

for (typ in type) {
    first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
    second_prefix <- ifelse(typ == "all", "02", "01")
    
    dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
    hub_info <- read.csv(paste0(dir_path, '/', typ, '_network_hub_info.csv'),row.names = 1)
    
    raw_df <- merge(hub_info, core_taxonomy, by = 'row.names')
    names(raw_df)[2:4] <- c('Degree', 'Closeness', 'Betweenness')
    
    for (kin in kingdom) {
        tmp_df <- raw_df %>%
            filter(Clade == kin) %>%
            select(Group = Order, Value = all_of(fea), Role = roles) %>%
            group_by(Group) %>%
            filter(n() >= 3) %>%
            ungroup()
        
        color_df <- tibble(
            Taxonomy = top_order[, kin],
            Color = colorRampPalette(brewer.pal(9, 'Set1'))(num + 1)[1:num]
        )
        
        data_df <- tmp_df %>%
            filter(Group %in% color_df$Taxonomy) %>%
            mutate(Group = factor(Group, levels = color_df$Taxonomy))
        
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
        
        network <- ifelse(typ == "all", "Whole network", paste0(str_to_title(typ), "kingdom network"))
        clade <- ifelse(kin == "Protist", "Protists", kin)
        
        data_df$Network <- network
        summarise_df$Network <- network
        pvalue_df$Network <- network
        
        data_df$Kingdom <- clade
        summarise_df$Kingdom <- clade
        pvalue_df$Kingdom <- clade
        
        all_data_df <- bind_rows(all_data_df, data_df)
        all_summarise_df <- bind_rows(all_summarise_df, summarise_df)
        all_pvalue_df <- bind_rows(all_pvalue_df, pvalue_df)
        
    }
    
}

name <- paste0(dir_name, "/degree_by_order")
write.csv(all_data_df, paste0(name, "_data.csv"), row.names = F)
write.csv(all_summarise_df, paste0(name, "_summarise.csv"), row.names = F)
write.csv(all_pvalue_df, paste0(name, "_pvalue.csv"), row.names = F)


network_levels <- c("Whole network", "Interkingdom network", "Intrakingdom network")
kingdom_levels <- c("Bacteria", "Fungi", "Protists")

all_data_df <- all_data_df %>%
    mutate(
        Network = factor(Network, levels = network_levels),
        Kingdom = factor(Kingdom, levels = kingdom_levels)
    )

all_summarise_df <- all_summarise_df %>%
    mutate(
        Network = factor(Network, levels = network_levels),
        Kingdom = factor(Kingdom, levels = kingdom_levels)
    )

all_hub_data_df <- all_data_df %>%
    filter(Role != "Peripherals")
all_nonhub_data_df <- all_data_df %>%
    filter(Role == "Peripherals")

p <- ggplot() + 
    geom_point(data = all_nonhub_data_df, aes(x = Group, y = Value, fill = "Group"), 
               position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.4, size = 2, stroke = 0, color = "#999999") +
    geom_point(data = all_hub_data_df, aes(x = Group, y = Value, fill = "Group"), 
               position = position_jitterdodge(dodge.width = 0.6), 
               alpha = 0.8, size = 2, stroke = 0, color = "red") +
    
    geom_boxplot(data = all_data_df, aes(x = Group, y = Value), width = 0.3, alpha = 0.2, na.rm = TRUE) +   
    geom_violin(data = all_data_df, aes(x = Group, y = Value), width = 0.5, alpha = 0.2, na.rm = TRUE) +  
    geom_text(
        data = all_summarise_df,
        mapping = aes(x = Group, y = max + allmax * 0.1, label = label),
        position = position_dodge(0.9),
        size = 7 / 2.835
    ) +
    labs(
        title = NULL,
        subtitle = NULL,
        x = NULL,
        y = "Degree centrality"
    ) + 
    facet_grid2(
        rows = vars(Network), 
        cols = vars(Kingdom), 
        scales = "free",
        independent = "y",
        space = "free_x"
    ) +
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

width <- 16
height <- 16
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
# ------------------------------------------------------------------------------