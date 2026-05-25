### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "04-barplot_data"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)


### Define variable -----------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
method <- c("relative", "absolute")
tax_level7 <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_level8 <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
for (met in method) {
    for (amp in amplicon) {
        file_name <- paste0("../00-rawdata/feature_table/", amp, '/fin_feature_table_', met, '.tsv')
        data_df <- read.csv(file_name, sep = '\t', row.names = 1)
        tax_df <- read.csv(paste0("../01-sort_data/02-taxonomy/", amp, "_ASV_taxonomy.csv"), row.names = 1)
        
        if (amp == "Protist") {
            tax_level <- tax_level8
        } else {
            tax_level <- tax_level7
        }
        value_cols <- names(data_df)[!names(data_df) %in% tax_level]
        
        
        # 1. 将Feature Table按照注释信息分类汇总
        raw_df <- merge(data_df, tax_df[tax_level], by = "row.names")
        raw_df <- raw_df[-1]
        
        res_df <- raw_df %>%
            group_by(across(Kingdom:Species)) %>%
            summarise(
                across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop"
            )
        
        name <- paste0(dir_name, "/", amp, "_barplot_data_", met, ".csv")
        write.csv(res_df, name, row.names = F)
        
        # 2. 将Feature Table按照每个分类层级的注释信息分别分类汇总
        all_df <- data.frame()
        for (level in tax_level) {
            sum_df <- res_df %>%
                group_by(across(all_of(level))) %>%
                summarise(across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
            
            level_summary <- sum_df %>%
                mutate(
                    Mean = rowMeans(select(., all_of(value_cols)), na.rm = TRUE),
                    Level = level
                ) %>%
                rename(Taxonomy = all_of(level)) %>%
                select(Taxonomy, Level, Mean, everything())
            
            all_df <- bind_rows(all_df, level_summary)
        }
        
        name <- paste0(dir_name, "/", amp, "_barplot_data_", met, "_sorted.csv")
        write.csv(all_df, name, row.names = F)
    }
}
# ------------------------------------------------------------------------------
