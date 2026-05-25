### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "05-core_table"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
amplicon <- c("16S", "ITS", "Protist")
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. absolute
core_table_df <- data.frame()
for (amp in amplicon) {
    tmp_df <- read.delim(paste0("../00-rawdata/feature_table/", amp, "/core_feature_table_absolute.tsv"))
    core_table_df <- bind_rows(core_table_df, tmp_df)
}

tax_df <- read.csv("02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)

row.names(core_table_df) <- tax_df[match(core_table_df$FeatureID, tax_df$FeatureID), 'ASVID']

data_df <- core_table_df[-1]

name <- paste0(dir_name, "/All_core_feature_table_absolute")
write.csv(data_df, paste0(name, ".csv"), quote = F, row.names = T)


# 2. relative
## 2.1 将所有文件读取到一个列表中，而不是逐个合并
df_list <- lapply(amplicon, function(amp) {
    read.delim(paste0("../00-rawdata/feature_table/", amp, "/core_feature_table_relative.tsv"))
})

## 2.2 获取所有数据框共有的列名 (names)
common_cols <- Reduce(intersect, lapply(df_list, names))

## 2.3 仅保留每个数据框中属于 common_cols 的列，最后一次性合并
core_table_df <- bind_rows(lapply(df_list, function(df) df[, common_cols, drop = FALSE]))

tax_df <- read.csv("02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)

row.names(core_table_df) <- tax_df[match(core_table_df$FeatureID, tax_df$FeatureID), 'ASVID']

data_df <- core_table_df[-1]

name <- paste0(dir_name, "/All_core_feature_table_relative")
write.csv(data_df, paste0(name, ".csv"), quote = F, row.names = T)
# ------------------------------------------------------------------------------