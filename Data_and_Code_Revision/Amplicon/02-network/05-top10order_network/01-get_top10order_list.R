### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "01-get_top10order_list"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
metadata_rs <- read.csv('../../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv')
absolute_abundance <- read.csv('../../01-sort_data/03-absolute_abundance/all_sample_copies.csv')

# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------
metadata_core <- metadata_rs %>%
    filter(FileID %in% absolute_abundance$FileID)

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. 获取 Top 10 Plant Order的统计信息
top_summary <- metadata_core %>%
    count(Order, sort = TRUE) %>%
    slice_head(n = 10) 

target_orders <- top_summary$Order

data_df <- metadata_core %>%
    filter(Order %in% target_orders) %>%
    group_by(Order) %>%
    mutate(row_id = row_number()) %>%
    ungroup() %>%
    select(row_id, Order, FileID) %>%
    pivot_wider(names_from = Order, values_from = FileID) %>%
    select(all_of(target_orders))


write.csv(data_df, paste0(dir_name, '/', 'top10order_list.csv'), row.names = F)
# ------------------------------------------------------------------------------


