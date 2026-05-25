# ******************************************************************************
# @File: 09-ko_distribution.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-12 19:22:43
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
dir_name <- "08-ko_distribution"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"
net <- "inter"

# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
color_manual <- colorRampPalette(brewer.pal(3, 'Accent'))(3)

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
all_hub_df <- NULL
all_ko_df <- NULL

for (ord in names(top10order_list)) {
  file_path <- paste0("03-get_top10order_network_property/", ord, "_", net, "_network_hub_info.csv")
  hub_info <- read.csv(file_path, row.names = 1)
  
  hub_df <- data.frame(Order = ord, KO = rownames(hub_info)[hub_info$roles != "Peripherals"])
  ko_df <- data.frame(Order = ord, KO = rownames(hub_info))
  
  all_hub_df <- bind_rows(all_hub_df, hub_df)
  all_ko_df <- bind_rows(all_ko_df, ko_df)
}

all_hub_df <- all_hub_df %>%
  mutate(Kingdom = case_when(
    grepl("^b", KO) ~ "Bacteria", 
    grepl("^f", KO) ~ "Fungi",
    TRUE ~ "Protists"
  ))

name <- paste0(dir_name, "/top10order_hub_ko")
write.csv(all_hub_df, paste0(name, "_data.csv"), quote = F, row.names = F)

all_ko_df <- all_ko_df %>%
  mutate(Kingdom = case_when(
    grepl("^b", KO) ~ "Bacteria", 
    grepl("^f", KO) ~ "Fungi",
    TRUE ~ "Protists"
  ))

name <- paste0(dir_name, "/top10order_all_ko")
write.csv(all_hub_df, paste0(name, "_data.csv"), quote = F, row.names = F)

all_ko_df$Type <- "All KOs"
all_hub_df$Type <- "Hub KOs"

result_df <- bind_rows(all_hub_df, all_ko_df)

data_df <- result_df %>%
  count(Order, Kingdom, Type, name = "Count") %>%
  mutate(Order = factor(Order, levels = names(top10order_list)))

label_df <- result_df %>%
  count(Order, Type, name = "Count")

p <- ggplot(data_df, aes(x = Order, y = Count)) +
  geom_bar(aes(fill = Kingdom), stat = "identity", position = position_fill(reverse = F), width = 0.68) +
  geom_text(
    data = label_df,
    mapping = aes(x = Order, y = 1.05, label = Count),
    position = position_dodge(0.9),
    size = 5 / 2.835
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Percentage of KO numbers"
  ) + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_manual, name = "Kingdom") +
  facet_grid(cols = vars(Type), scales = "free", space = "free_x") +
  
  theme_bw() + theme(
    text = element_text(color = "black", size = 5),
    plot.title = element_text(size = 7, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 7),
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.box.spacing = unit(0.01,"cm"),
    strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
    panel.spacing = unit(0, "cm"),
    plot.margin = margin(5.5, 1, 5.5, 1, "pt"),
    legend.position = "none"
  )

name <- paste0(dir_name, "/top10order_ko_distribution")
width <- 11.5
height <- 5.5
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")
write.csv(data_df , paste0(name, "_data_df .csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------
