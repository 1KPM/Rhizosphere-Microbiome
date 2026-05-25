### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "08-network_color"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
ratio_threshold <- 0.05
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
core_table <- read.csv("05-core_table/All_core_feature_table_absolute.csv", row.names = 1)
core_taxonomy <- read.csv("02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. 计算所有Core Feature的Sum并与tax合并
sum_df <- data.frame(Sum = rowSums(core_table)) %>%
    rownames_to_column(var = "ASVID")

top_tax <- core_taxonomy  %>%
    inner_join(sum_df, by = "ASVID")

# 2. 定义通用函数：提取指定 Clade 和 分类层级 (Phylum/Order) 的 Top taxa
get_top_taxa <- function(data, clade_name, tax_level) {
    data %>%
        filter(Clade == clade_name) %>%
        group_by(.data[[tax_level]]) %>%
        summarise(Sum = sum(Sum, na.rm = TRUE), .groups = 'drop') %>%
        mutate(Ratio = Sum / sum(Sum, na.rm = TRUE)) %>%
        # 过滤掉空名称，并保留 Ratio >= 5 的类群
        filter(.data[[tax_level]] != "", Ratio >= ratio_threshold) %>%
        arrange(desc(Sum)) %>%
        pull(.data[[tax_level]]) # 提取为向量
}

# 3. 批量获取结果，生成列表
top_5percent_taxa <- list(
    top_16s_phylum = get_top_taxa(top_tax, "Bacteria", "Phylum"),
    top_16s_order  = get_top_taxa(top_tax, "Bacteria", "Order"),
    top_its_phylum = get_top_taxa(top_tax, "Fungi", "Phylum"),
    top_its_order  = get_top_taxa(top_tax, "Fungi", "Order"),
    top_pro_phylum = get_top_taxa(top_tax, "Protist", "Phylum"),
    top_pro_order  = get_top_taxa(top_tax, "Protist", "Order")
)

# 4. 统一列表长度，转换为 Data Frame 并输出
max_len <- max(lengths(top_5percent_taxa))
top_5percent_df <- as.data.frame(lapply(top_5percent_taxa, function(x) {
    length(x) <- max_len # 长度不足的自动用 NA 补齐
    return(x)
}))

write.csv(top_5percent_df, paste0(dir_name, '/top_5percent_taxa.csv'), row.names = F, quote = F, na = '')
# ------------------------------------------------------------------------------

