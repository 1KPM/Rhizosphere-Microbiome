# ******************************************************************************
# @File: 05-get_stable_hub_ko.R
# @Author: Mingxing Wang
# @Email: xing592798030@163.com
# @Date: 2026-03-10 16:56:05
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
dir_name <- "04-get_stable_hub_ko"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)
library(UpSetR)
library(ggvenn)

### Define variable -----------------------------------------------------------
p_adjust_method <- "fdr"

# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
top10order_list <- read.csv("01-get_top10order_list/top10order_list.csv")
abundance <- read.csv("../01-all_network/00-data/All_core0.2_quantitative.csv",row.names = 1,check.names = F)
# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
net <- "inter"

# 1. UpSetR
tmp_list <- NULL
for (ord in names(top10order_list)) {
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
        # 获取该元素所属的所有集合名称，用 "; " 连接
        Intersections = paste(names(tmp_data)[c_across(-Element_ID) == 1], collapse = "; ")
    ) %>%
    ungroup()

# 生成汇总表：每一类交集包含哪些元素
summary_table <- extract_df %>%
    group_by(Intersections) %>%
    summarise(
        Count = n(),
        # 将该组下的所有 ID 用分号连接，方便在 Excel 查看
        Members = paste(Element_ID, collapse = "; ") 
    ) %>%
    arrange(desc(Count))

# 在10个目中都为core ko的定义为common hub ko
common_hub_ko <- summary_table$Members[str_count(summary_table$Intersections, ";") == 9]
common_hub_ko <- unlist(strsplit(common_hub_ko, split = "; "))

name <- paste0(dir_name, "/", net, "_top10order")
write.csv(summary_table, paste0(name, "_upsetR_results.csv"), quote = F, row.names = F)
write.csv(common_hub_ko, paste0(name, "_common_hub_ko.csv"), quote = F, row.names = F)

combinations <- summary_table[str_count(summary_table$Intersections, ";") >= 7, "Intersections"] 
combinations <- strsplit(combinations$Intersections, "; ")

p <- upset(tmp_data, 
           nsets = 10,
           order.by = "freq",
           mb.ratio = c(0.55, 0.45),
           main.bar.color = "steelblue",
           sets.bar.color = "darkorange",
           intersections = c(as.list(names(tmp_data)[1:10]), combinations),
           text.scale = c(7/8, 6/7, 7/8, 6/7, 6/7, 1),
           point.size = 1,
           line.size = 0.5)
width <- 16
height <- 8
png(paste0(name, "_UpSetR.png"), width = width, height = height, units = "cm", res = 600); print(p); dev.off()
pdf(paste0(name, "_UpSetR.pdf"), width = width / 2.54, height = height / 2.54); print(p); dev.off()


# 2. Pie Plot
file_path <- switch(
    net,
    "all"   = "../01-all_network/02-get_all_network_property/all_network_hub_info.csv",
    "inter" = "../02-inter_network/01-get_inter_network_property/inter_network_hub_info.csv",
    "../03-intra_network/01-get_intra_network_property/intra_network_hub_info.csv"
)

raw_hub <- read.csv(file_path, row.names = 1)
raw_hub_ko <- row.names(raw_hub[raw_hub$roles != "Peripherals",])


# 明确背景全集
background_ko <- rownames(abundance)

# 严谨计算 2x2 列联表的四个象限 (Quadrant)
# [1] 既是 raw_hub 又是 common_hub (A ∩ B)
n_both <- length(intersect(raw_hub_ko, common_hub_ko))

# [2] 是 raw_hub，但不是 common_hub (A - B)
n_raw_not_common <- length(setdiff(raw_hub_ko, common_hub_ko))

# [3] 是 common_hub，但不是 raw_hub (B - A)
n_common_not_raw <- length(setdiff(common_hub_ko, raw_hub_ko))

# [4] 既不是 raw_hub，也不是 common_hub (全集 - (A ∪ B))
union_ko <- union(raw_hub_ko, common_hub_ko)
n_neither <- length(setdiff(background_ko, union_ko))

# 3. 构建列联表 (Contingency Table)
tmp_df <- data.frame(
    'Common_hub' = c(n_both, n_common_not_raw),
    'Non_common_hub' = c(n_raw_not_common, n_neither),
    row.names = c('Raw_hub', 'Non_raw_hub')
)

# 4. 运行卡方检验
sig_res <- chisq.test(tmp_df)

# （顺带修复下方绘图代码中未定义的变量名）
data_list <- list()
data_list[['Raw hub']] <- raw_hub_ko 
data_list[['Common hub']] <- common_hub_ko

p <- ggvenn(data = data_list,
            show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
            label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
            show_percentage = T,      # 显示每一组的百分比
            digits = 1,               # 百分比的小数点位数
            fill_color = c("#E41A1C", "#1E90FF"), # 填充颜色
            fill_alpha = 0.5,         # 填充透明度
            stroke_color = "white",   # 边缘颜色
            stroke_alpha = 0.5,       # 边缘透明度
            stroke_size = 0.5,        # 边缘粗细
            stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
            set_name_color = "black", # 组名颜色
            set_name_size = 7 / 2.835,        # 组名大小
            text_color = "black",     # 交集个数颜色
            text_size = 6 / 2.835             # 交集个数文字大小
) +
    labs(title = ifelse(sig_res$p.value < 0.001, 'p < 0.001\n', 
                       ifelse(sig_res$p.value < 0.01, 'p < 0.01\n', 
                              ifelse(sig_res$p.value < 0.05, 'p < 0.05\n', 'p > 0.05\n'))),
        subtitle = NULL
    ) + 
    theme(
        plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 7, color = 'black', hjust = 0.5)
    )

name <- paste0(dir_name, "/raw_vs_common_venn")
width <- 6
height <- 6
ggsave(paste0(name, '.png'), p, width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
ggsave(paste0(name, '.pdf'), p, width = width, height = height)

stable_hub_ko <- intersect(raw_hub_ko, common_hub_ko)
write.csv(stable_hub_ko, paste0(dir_name, "/", net, "_stable_hub_ko.csv"), quote = F, row.names = F)
# ------------------------------------------------------------------------------


