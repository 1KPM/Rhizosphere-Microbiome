# ******************************************************************************
# @File: 03-robustness.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-03 10:04:59
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
dir_name <- "03-robustness"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

source('get_network_robustness.R')

### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

type <- c("all", "inter", "intra")
target_clades <- c('Bacteria', 'Fungi', 'Protists')
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)
# ------------------------------------------------------------------------------

### Get results ----------------------------------------------------------------
for (typ in type) {
    first_prefix <- ifelse(typ == "all", "01", ifelse(typ == "inter", "02", "03"))
    second_prefix <- ifelse(typ == "all", "02", "01")
    
    dir_path <- paste0("../", first_prefix, "-", typ, "_network/", second_prefix, "-get_", typ, "_network_property")
    
    hub_info <- read.csv(paste0(dir_path, '/', typ, '_network_hub_info.csv'),row.names = 1)
    r_table <- read.csv(paste0(dir_path, '/', typ, '_r_table.csv'), row.names = 1)
    ko_table <- read.csv(paste0(dir_path, '/', typ, '_ko_table.csv'), row.names = 1)
    
    hub_info <- hub_info %>%
      mutate(Clade = case_when(
        substr(.[[1]], 1, 1) == "b" ~ "Bacteria",
        substr(.[[1]], 1, 1) == "f" ~ "Fungi",
        substr(.[[1]], 1, 1) == "p" ~ "Protist",
        TRUE ~ NA_character_
      ))
    
    
    hub_info <- hub_info %>%
      mutate(Clade = case_when(Clade == "Protist" ~ "Protists", TRUE ~ Clade),
             Clade = factor(Clade, levels = c('Bacteria', 'Fungi', 'Protists')))
    
    tmp_df <- hub_info %>%
        filter(roles != 'Peripherals') %>%
        arrange(desc(degree))
    
    # 1. By cumulative KO
    data_df <- data.frame(row.names = 1:nrow(tmp_df))
    for (i in 1:nrow(tmp_df)) {
        target_remove_robustness <-
            get_target_remove_robustness(r_table, ko_table, rm_feature = tmp_df$name[1:i], table_type = 'absolute')
        random_remove_robustness_df <-
            get_random_remove_robustness(r_table, ko_table, rm_percent = i/nrow(r_table), table_type = 'absolute')
        random_remove_robustness <- mean(random_remove_robustness_df$remain_percen)
        
        data_df[i, c('target', 'random')] <- c(target_remove_robustness, random_remove_robustness)
    }
    
    name <- paste0(dir_name, "/", typ, "_robustness_by_cumulative_ko")
    write.csv(data_df, paste0(name, ".csv"), quote = F)
    
    # 2. By Kingdom
    data_df <- data.frame(row.names = target_clades)
    for (cla in target_clades) {
        rm_feature <- tmp_df[tmp_df$Clade == cla, 'name']
        n <- length(rm_feature)
        target_remove_robustness <-
            get_target_remove_robustness(r_table, ko_table, rm_feature = rm_feature, table_type = 'absolute')
        random_remove_robustness_df <-
            get_random_remove_robustness(r_table, ko_table, rm_percent = n/nrow(r_table), table_type = 'absolute')
        random_remove_robustness <- mean(random_remove_robustness_df$remain_percen)
        data_df[cla, c('number', 'target', 'random')] <- c(n, target_remove_robustness, random_remove_robustness)
    }
    
    name <- paste0(dir_name, "/", typ, "_robustness_by_kingdom")
    write.csv(data_df, paste0(name, ".csv"), quote = F)
}

# ------------------------------------------------------------------------------