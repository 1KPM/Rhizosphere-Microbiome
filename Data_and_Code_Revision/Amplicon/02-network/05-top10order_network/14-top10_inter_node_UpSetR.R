### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "14-top10_inter_node_UpSe"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(clusterProfiler)
library(patchwork)
library(ggh4x)
library(UpSetR)
library(vegan)
library(data.table)
library(zoo)
library(seriation)


### Get results ----------------------------------------------------------------

# 1. Figure5A
net <- "inter"


All_list <- c("Fabales","Arecales","Malpighiales","Sapindales","Lamiales",
              "Rosales","Asparagales","Myrtales","Malvales","Gentianales")

tmp_list <- NULL
for (ord in All_list) {
  file_path <- paste0("03-get_top10order_network_property/", ord, "_", net, "_network_hub_info.csv")
  hub_info <- read.csv(file_path, header = T, row.names = 1)
  hub_ko <- row.names(hub_info[hub_info$roles != 'Peripherals', ])
  tmp_list[[ord]] <- hub_ko
}

all_elements <- unique(unlist(tmp_list))
tmp_data <- fromList(tmp_list)
rownames(tmp_data) <- all_elements

extract_df <- tmp_data %>%
  mutate(Element_ID = rownames(.)) %>%
  rowwise() %>%
  mutate(
    Intersections = paste(names(tmp_data)[c_across(-Element_ID) == 1], collapse = "; ")
  ) %>%
  ungroup()

summary_table <- extract_df %>%
  group_by(Intersections) %>%
  summarise(
    Count = n(),
    Members = paste(Element_ID, collapse = "; ") 
  ) %>%
  arrange(desc(Count))

# 修改为9个分号（10个orders之间会有9个分号）
common_hub_ko <- summary_table$Members[str_count(summary_table$Intersections, ";") == 9]
common_hub_ko <- unlist(strsplit(common_hub_ko, split = "; "))

name <- paste0(dir_name, "/", net, "_top10order_node")
write.csv(summary_table, paste0(name, "_upsetR_results.csv"), quote = F, row.names = F)
write.csv(common_hub_ko, paste0(name, "_common_hub_ko.csv"), quote = F, row.names = F)

# 修改为至少7个分号（至少存在于8个orders中）
combinations <- summary_table[str_count(summary_table$Intersections, ";") >= 7, "Intersections"] 
combinations <- strsplit(combinations$Intersections, "; ")

combinations3 <- strsplit(summary_table$Intersections[1:43], "; ")
combinations3 <- c(combinations3, combinations[-1])
target_indices <- c(42, 44:50)
queries_list <- lapply(target_indices, function(idx) {
    list(
        query = intersects,
        params = list(combinations3[[idx]]),
        color = ifelse(idx == 42, "#E64B35", "#3c5488"),  # 这里固定颜色即可
        active = TRUE
    )
})

p1 <- upset(tmp_data,
                      sets = All_list,
                      keep.order = TRUE,
                      order.by = "freq",
                      mb.ratio = c(0.45, 0.55),
                      main.bar.color = "#3C5488",
                      sets.bar.color = "#7E6148",
                      matrix.color = "#2F2F2F",
                      shade.color = "#E0E0E0",
                      shade.alpha = 0.25,
                      intersections = combinations3,
                      text.scale = c(0.6, 0.8, 0.6, 0.8, 1.0, 0.8),
                      point.size = 1.0,
                      line.size = 0.7,
                      mainbar.y.label = NULL,
                      sets.x.label = NULL,
                      queries = queries_list,
                      query.legend = "none"
)

width <- 17.5
height <- 6.5
name <- paste0(dir_name, "/", net, "_node")
png(paste0(name, "_UpSetR.png"), width = width, height = height, units = "cm", res = 600); print(p1); dev.off()
pdf(paste0(name, "_UpSetR.pdf"), width = width / 2.54, height = height / 2.54); print(p1); dev.off()
