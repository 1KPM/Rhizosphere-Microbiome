### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "09-enrichment_database"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
core_taxonomy <- read.csv("02-taxonomy/All_core_ASV_taxonomy.csv", row.names = 1)

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
order2asv <- core_taxonomy %>%
    select(Order, ASVID)
write.csv(order2asv, paste0(dir_name, "/order2asv.csv"), quote = F, row.names = F)

order2kingdom <- core_taxonomy %>%
    select(Order, Clade) %>%
    rename(Kingdom = "Clade") %>%
    distinct() %>%
    mutate(Kingdom = case_when(Kingdom == "Protist" ~ "Protists", TRUE ~ Kingdom))
write.csv(order2kingdom, paste0(dir_name, "/order2kingdom.csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------