### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "07-top_order"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
kingdom <- c('Bacteria', 'Fungi', 'Protist')
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
tax_all <- read.csv('02-taxonomy/All_core_ASV_taxonomy.csv', row.names = 1)
# ------------------------------------------------------------------------------
num = 10
top_order <- data.frame(row.names = 1:num)
for (k in kingdom) {
    raw_df <- tax_all[tax_all$Clade == k,]
    order_df <- data.frame(t(table(raw_df$Order)))
    order_df <- order_df[order(-order_df$Freq),]
    order_df <- order_df[!(order_df$Var2 %in% c('')),]
    top_order[k] <- order_df$Var2[1:num]
}
write.csv(top_order, paste0(dir_name, '/top_order.csv'), row.names = F)

