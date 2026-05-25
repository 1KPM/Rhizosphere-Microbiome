### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "03-network_plot"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(igraph)
library(ggClusterNet)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

treatment <- c("B", "BF", "BFP", "All")
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
network_color <- read.csv("../01-sort_data/05-network_color/top_5percent_taxa.csv", header = T)

# ------------------------------------------------------------------------------

### Sort data ----------------------------------------------------------------
seed_df <- data.frame(
    Treatment = treatment, 
    Seed = c(1995, 1996, 1995, 1997)
)

# 定义一个简易清理函数：去除向量中的 NA 和 空字符串 ("")
clean_taxa <- function(x) {
    x[!is.na(x) & x != ""]
}

# 批量清理并直接释放为全局变量！
# lapply(network_color, clean_taxa) 会返回一个干净的列表。
# list2env() 会把这个列表中的元素，直接按照列名（如 top_16s_phylum）变成你环境里的独立变量。
list2env(lapply(network_color, clean_taxa), envir = .GlobalEnv)
# ------------------------------------------------------------------------------



### Get results ----------------------------------------------------------------
for (tre in treatment) {
    tmp_edge <- read.csv(paste0("01-get_network_data/", tre, "_edge_property_raw.csv"))
    tmp_node <- read.csv(paste0("01-get_network_data/", tre, "_node_property_raw.csv"))
    
    otu_table <- read.csv(paste0("02-get_network_property/",tre, "_otu_table.csv"), check.names = F, row.names = 1)
    r_table <- read.csv(paste0("02-get_network_property/",tre, "_r_table.csv"), row.names = 1)
    tax_table <- read.csv(paste0("02-get_network_property/",tre, "_tax_table.csv"), row.names = 1)
    
    igraph <- graph_from_data_frame(tmp_edge, directed = F, vertices = tmp_node)
    node_property <- node_properties(igraph)
    
    # 按照maptree内聚算法改进离散点排布
    set.seed(seed_df$Seed[seed_df$Treatment == tre])
    r_matrix <- as.matrix(r_table)
    maptree <- model_maptree2(cor = r_matrix, method = 'cluster_fast_greedy')
    maptree_node <- maptree[[1]]
    maptree_edge <- edgeBuild(cor = r_matrix, node = maptree_node)
    names(maptree_edge)[names(maptree_edge) == 'cor'] <- 'Correlation'
    
    # node节点注释
    maptree_node <- nodeadd(plotcord = maptree_node, otu_table = otu_table, tax_table = tax_table)
    node_property <- data.frame(node_property)
    names(node_property) <- c('Degree', 'Closeness', 'Betweenness', 'CenDegree')
    maptree_node <- merge(maptree_node, data.frame(node_property), by = 'row.names')
    
    maptree_node <- maptree_node %>%
        mutate(
            # 1. 整理 Phylum 因子：如果不在 top 列表中，则归为 Other
            Phylum = case_when(
                !(Phylum %in% top_16s_phylum) ~ 'Other bacteria',
                TRUE ~ Phylum  # 其余情况（即在 top 列表中的）保持原名
            ),
            
            # 2. 整理 Order 因子：如果不在 top 列表中，则归为 Other
            Order = case_when(
                !(Order %in% top_16s_order) ~ 'Other bacteria',
                TRUE ~ Order   # 其余情况保持原名
            ),
            
            # 3. 转换为因子并严格指定 levels 顺序
            Phylum = factor(Phylum, levels = c(
                top_16s_phylum, 'Other bacteria'
            )),
            
            Order = factor(Order, levels = c(
                top_16s_order, 'Other bacteria'
            ))
        )
    
    # 按Order着色
    color_order <- colorRampPalette(brewer.pal(12, 'Paired'))(7)

    p <- ggplot() + 
        geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Correlation), 
                     alpha = 0.01, linewidth = 0.1, data = maptree_edge) + 
        geom_point(aes(X1, X2, fill = Order, size = Degree ** 2), pch = 21, data = maptree_node) +
        scale_fill_manual(values = color_order) +
        scale_x_continuous(breaks = NULL) + 
        scale_y_continuous(breaks = NULL) +
        theme_classic() + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text = element_text(size = 6, color = 'black'), 
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              legend.position = 'right')

    width <- 24
    height <- 18
    name <- paste0(dir_name, '/', tre, '_network_maptree_Phylum')
    
    ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
    ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
}
# ------------------------------------------------------------------------------

