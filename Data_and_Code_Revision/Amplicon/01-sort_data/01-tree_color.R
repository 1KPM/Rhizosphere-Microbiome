### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "01-tree_color"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Import data ----------------------------------------------------------------
metadata_rs <- read.csv('../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv')
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
# 1. Top4 Class
fre_class <- data.frame(table(metadata_rs$Class))
fre_class <- fre_class[order(-fre_class$Freq), ]
top_class <- c(as.character(fre_class$Var1[1:4]), "Others")
class_color_manual <- c(colorRampPalette(brewer.pal(8, 'Paired'))(length(top_class)-1), "#999999")
class_color_df <- data.frame(Class = top_class, Color = class_color_manual)
name <- paste0(dir_name, "/Tree_top4_class_color.csv")
write.csv(class_color_df, name, row.names = F)

# 2. Top10 Order
fre_order <- data.frame(table(metadata_rs$Order))
fre_order <- fre_order[order(-fre_order$Freq), ]
top_order <- c(as.character(fre_order$Var1[1:10]), "Others")
order_color_manual <- colorRampPalette(brewer.pal(9, 'Set1'))(length(top_order))
order_color_df <- data.frame(Order = top_order, Color = order_color_manual)
name <- paste0(dir_name, "/Tree_top10_order_color.csv")
write.csv(order_color_df, name, row.names = F)
# ------------------------------------------------------------------------------
