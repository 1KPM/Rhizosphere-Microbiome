
### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "10-Upset_edge"
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


net <- "inter"


All_list <- c("Fabales","Arecales","Malpighiales","Sapindales","Lamiales",
              "Rosales","Asparagales","Myrtales","Malvales","Gentianales")


#-------------------------------------------------------------------------------
tmp_list <- NULL
for (ord in All_list) {
  file_path <- paste0("03-get_top10order_network_property/", ord, "_", net, "_edge_property_raw.csv")
  edge_mt2 <- read.csv(file_path, header = T)
  edge <- edge_mt2 %>% select(node1,node2,label)
  tmp_list[[ord]] <- edge
}

all_pairs <- do.call(rbind, lapply(tmp_list, function(df) df[, c("node1", "node2")]))
all_pairs <- as.data.frame(t(apply(all_pairs, 1, sort)))
colnames(all_pairs) <- c("node1", "node2")
unique_edges <- unique(all_pairs)

edge_names <- names(tmp_list)

edge_maps <- map(tmp_list, function(df) {
  df %>%
    rowwise() %>%
    mutate(edge_key = paste(sort(c(node1, node2)), collapse = "_")) %>%
    ungroup() %>%
    select(edge_key, label) %>%
    deframe()
})

unique_edges <- unique_edges %>%
  rowwise() %>%
  mutate(edge_key = paste(sort(c(node1, node2)), collapse = "_")) %>%
  ungroup()

setDT(unique_edges)
unique_edges[, edge_key := paste(pmin(node1, node2), pmax(node1, node2), sep = "_")]

edge_names <- names(tmp_list)

for (nm in edge_names) {
  dt <- as.data.table(tmp_list[[nm]])[, .(node1, node2, label)]
  dt[, edge_key := paste(pmin(node1, node2), pmax(node1, node2), sep = "_")]
  dt[, label := fcase(label == "+", 1L,
                      label == "-", -1L,
                      default = 0L)]
  
  setkey(unique_edges, edge_key)
  setkey(dt, edge_key)
  unique_edges <- merge(unique_edges, dt[, .(edge_key, label)], 
                        by = "edge_key", all.x = TRUE, sort = FALSE)
  setnames(unique_edges, "label", nm)
  unique_edges[is.na(get(nm)), (nm) := 0L]
}

unique_edges[, edge_key := NULL]

write.csv(unique_edges,paste0(dir_name,'/unique_edges.csv'),row.names = F)

all_edge_result <- read.csv(paste0(dir_name,'/unique_edges.csv'))
tmp_data <- all_edge_result %>% 
  mutate(rowname = paste(node1, node2, sep = "-")) %>%
  column_to_rownames("rowname") %>%
  select(3:12) %>%  # 改为3:12，因为只有10个orders
  mutate(across(everything(), ~ifelse(.x == -1, 0, .x)))

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

summary_table <- summary_table[summary_table$Intersections != "", ]

# 修改为9个分号（10个orders之间会有9个分号）
common_edge <- summary_table$Members[str_count(summary_table$Intersections, ";") == 9]
common_edge <- unlist(strsplit(common_edge, split = "; "))

name <- paste0(dir_name, "/", net, "_top10order_edge")
write.csv(summary_table, paste0(name, "_upsetR_results.csv"), quote = F, row.names = F)
write.csv(summary_table[,1:2], paste0(name, "_upsetR_results_simple.csv"), quote = F, row.names = F)
write.csv(common_edge, paste0(name, "_common_edge.csv"), quote = F, row.names = F)

# 修改为至少8个分号（至少存在于9个orders中）
combinations <- strsplit(summary_table$Intersections[1:50], "; ")


p2 <- upset(tmp_data,
                     sets = All_list,
                     keep.order = TRUE,
                     order.by = "freq",
                     mb.ratio = c(0.45, 0.55),
                     main.bar.color = "#3C5488",
                     sets.bar.color = "#7E6148",
                     matrix.color = "#2F2F2F",
                     shade.color = "#E0E0E0",
                     shade.alpha = 0.25,
                     number.angles = 90,
                     # scale.intersections = "log2",
                     intersections = combinations,
                     text.scale = c(0.6, 0.8, 0.6, 0.8, 1.0, 0.7),
                     point.size = 1.0,
                     line.size = 0.7,
                     mainbar.y.label = NULL,
                     # mainbar.y.label = "",   # 修改为字符串
                     sets.x.label = NULL,
                     queries = list(
                       list(
                         query = intersects,
                         params = list(All_list),
                         color = "#E64B35",
                         active = TRUE
                       )
                     ),
                     query.legend = "none"
)

width <- 18
height <- 6.5
name <- paste0(dir_name, "/", net, "_edge")
png(paste0(name, "_UpSetR.png"), width = width, height = height, units = "cm", res = 600); print(p2); dev.off()
pdf(paste0(name, "_UpSetR.pdf"), width = width / 2.54, height = height / 2.54); print(p2); dev.off()
