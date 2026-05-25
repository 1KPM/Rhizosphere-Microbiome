# ******************************************************************************
# @File: 01-get_absolute_alpha_diversity.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-01-17 15:18:37
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
dir_name <- "01-get_absolute_alpha_diversity"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(ape)
library(vegan)
library(picante)

# Define function
alpha_diversity <- function(df, tree = NULL, base = exp(1)) {
    est <- estimateR(df)
    ASV <- est[1,]
    Chao <- est[2,]
    ACE <- est[4,]
    Shannon <- diversity(df, index = 'shannon', base = base)
    Simpson <- diversity(df, index = 'simpson')
    Pielou <- Shannon / log(ASV, base)
    Coverage <- 1 - rowSums(df == 1) / rowSums(df)
    res_df <- data.frame(ASV, Chao, ACE, Shannon, Simpson, Pielou, Coverage)
    
    if (!is.null(tree)) {
        PD_whole_tree <- pd(df, tree, include.root = TRUE)[1]
        names(PD_whole_tree) <- 'PD'
        res_df <- cbind(PD_whole_tree, res_df)
    }
    
    return(res_df)
}
# ------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
kin = args[1]


### Get results ----------------------------------------------------------------
copies_df <- read.delim(
    paste0("../../Rawdata/feature_table/", kin, "/fin_feature_table_absolute.tsv"), row.names = 1, check.names = FALSE)
rooted_tree <- read.tree(paste0("../../Rawdata/phylogeny/", kin, "/tree.nwk"))

absolute_copies_df <- as.data.frame(t(copies_df))
absolute_div_df <- alpha_diversity(absolute_copies_df, tree = rooted_tree)
absolute_div_df <- absolute_div_df[c('PD', 'ASV', 'Chao', 'Shannon', 'Simpson', 'Pielou')]
name <- paste0(dir_name, "/", kin, "_alpha_diversity_absolute.csv")
write.csv(absolute_div_df, name)
# ------------------------------------------------------------------------------


