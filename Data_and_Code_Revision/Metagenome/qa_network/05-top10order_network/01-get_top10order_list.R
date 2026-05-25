# ******************************************************************************
# @File: 01-get_top10order_list.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-03 11:13:24
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


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
metadata_rs <- read.csv('../01-all_network/00-data/metadata.csv')
absolute_abundance <- read.csv('../01-all_network/00-data/All_core0.2_quantitative.csv',check.names = F,row.names = 1)

# ------------------------------------------------------------------------------
metadata_core <- metadata_rs %>% filter(TreeID %in% colnames(absolute_abundance))

### Get results ----------------------------------------------------------------
# 1. 获取 Top 10 Plant Order的统计信息
target_orders <- c("Fabales","Rosales","Malpighiales","Gentianales","Sapindales","Lamiales","Malvales","Myrtales","Arecales","Asparagales")

data_df <- metadata_core %>%
  filter(Order %in% target_orders) %>%
  group_by(Order) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  select(row_id, Order, TreeID) %>%
  pivot_wider(names_from = Order, values_from = TreeID) %>%
  select(all_of(target_orders))


write.csv(data_df, paste0(dir_name, '/', 'top10order_list.csv'), row.names = F)
# ------------------------------------------------------------------------------


