# ******************************************************************************
# @File: 02-get_absolute_beta_diversity.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-01-17 15:30:04
# @License: Copyright (C) 2026 Mingxing Wang. All rights reserved.
# @Reference: Mingxing Wang
# @Description: 
# ******************************************************************************


### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- "./"
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "02-get_absolute_beta_diversity"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(ape)
library(vegan)
library(GUniFrac)


args = commandArgs(trailingOnly = TRUE)
kin = args[1]

### Get results ----------------------------------------------------------------
absolute_copies <- read.delim(
    paste0("../../Rawdata/feature_table/", kin, "/fin_feature_table_absolute.tsv"), row.names = 1, check.names = FALSE)
rooted_tree <- read.tree(paste0("../../Rawdata/phylogeny/", kin, "/tree.nwk"))

name <- paste0(dir_name, "/", kin, "_beta_diversity_absolute_")
absolute_df <- as.data.frame(t(absolute_copies))

absolute_bd <- vegdist(absolute_df, method = 'bray')
absolute_bd <- as.data.frame(as.matrix(absolute_bd))
write.csv(absolute_bd, paste0(name, '_bray_curtis.csv'))

# 这等同于 QIIME 2 中的 Jaccard
absolute_jd <- vegdist(absolute_df, method = 'jaccard', binary = TRUE)
absolute_jd <- as.data.frame(as.matrix(absolute_jd))
write.csv(absolute_jd, paste0(name, '_jaccard.csv'))

absolute_GUniFrac <- GUniFrac(absolute_df, rooted_tree, alpha=c(0, 0.5, 1))
absolute_weighted <-  absolute_GUniFrac$unifracs[, , 'd_1']
absolute_unweighted <-  absolute_GUniFrac$unifracs[, , 'd_UW']
absolute_generalized <-  absolute_GUniFrac$unifracs[, , 'd_0.5']


write.csv(absolute_weighted, paste0(name, '_weighted.csv'))
write.csv(absolute_unweighted, paste0(name, '_unweighted.csv'))
write.csv(absolute_generalized, paste0(name, '_generalized.csv'))

# ------------------------------------------------------------------------------
