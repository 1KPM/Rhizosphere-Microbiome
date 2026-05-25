### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(rstatix)
library(ggpubr)
library(agricolae)

# 读取Excel文件
data <- read.csv("Protistan_Fabales_vs_others.csv")  #整个代码将Bacterial替换为Fungal或Protistan即可
plant_colors <- read.csv("./data/Tree_top10_order_color.csv")

# 1. 统计每个比较组中鉴定为1Fabales和其他植物目的biomarker taxa数目

# 提取比较组列名（排除第一列Taxa）
comparison_cols <- colnames(data)[-1]

# 创建结果数据框
results <- data.frame(
  Comparison = character(),
  Fabales_count = integer(),
  Other_count = integer(),
  stringsAsFactors = FALSE
)

# 对每个比较组进行统计
for(col in comparison_cols) {
  # 统计该列中为"1Fabales"的个数
  fabales_count <- sum(data[[col]] == "1Fabales", na.rm = TRUE)
  
  # 统计该列中为其他植物目（非"1Fabales"且非空）的个数
  other_count <- sum(!is.na(data[[col]]) & data[[col]] != "1Fabales" & data[[col]] != "" & data[[col]] != "-")
  
  results <- rbind(results, data.frame(
    Comparison = col,
    Fabales_count = fabales_count,
    Other_count = other_count
  ))
}

# 查看统计结果
print(results)

# 2. 将数据转换为长格式用于绘图
results_long <- results %>%
  pivot_longer(
    cols = c(Fabales_count, Other_count),
    names_to = "Enrichment_Type",
    values_to = "Count"
  ) %>%
  mutate(Enrichment_Type = factor(Enrichment_Type, 
                                  levels = c("Fabales_count", "Other_count"),
                                  labels = c("Increased", "Decreased")))

results_long$Order <- gsub("vs", "", results_long$Comparison)
results_long$Order <- factor(results_long$Order, levels = plant_colors$Order[-1])


# 3. 进行配对t检验
# 提取Fabales和Other的计数
fabales_counts <- results$Fabales_count
other_counts <- results$Other_count

# 进行配对t检验
t_test_result <- t.test(fabales_counts, other_counts, paired = TRUE)

# 4. 绘制箱线图+散点图
p1 <- ggplot(results_long, aes(x = Enrichment_Type, y = Count)) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.3,
    outlier.shape = NA,
    fill = "lightblue"
  ) +
  geom_point(
    aes(color = Order),  
    position = position_jitter(width = 0.1, height = 0),
    size = 3,
    alpha = 1
  ) +
  # 手动设置点的颜色
  scale_color_manual(
    values = plant_colors$Color[-1], 
  ) +
  # 添加p值标注
  stat_compare_means(
    comparisons = list(c("Increased", "Decreased")),
    method = "t.test",
    paired = TRUE,
    label = "p.format",
    tip.length = 0,
    vjust = -0.5,
    size = 3
  ) +
  labs(
    title = "LEfSe Biomarkers",
    x = "",
    y = "Number of Taxa",
    color = "Fabales vs *"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),
    axis.title = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.height = unit(0.3, "cm"),  # 减小图例项的高度
    legend.key.width = unit(0.3, "cm"),   # 减小图例项的宽度
    legend.spacing.y = unit(0.1, "cm"),   # 减小图例项之间的垂直间距
  ) +
  # 设置y轴从-100开始
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    limits = c(-100, NA)  
  )

print(p1)

# 保存图形
ggsave("./plots/LEfSe_Protistan_taxa_biomarker.png", p1, 
       width = 7, height = 7, dpi = 300, bg = "white", units = "cm", limitsize = FALSE)
ggsave("./plots/LEfSe_Protistan_taxa_biomarker.pdf", p1, 
       width = 7, height = 7, bg = "white", units = "cm", limitsize = FALSE)


# 5. 统计前10个最频繁鉴定为1Fabales的biomarker taxa
# 创建一个数据框来记录每个taxa在所有比较组中被鉴定为1Fabales的次数
taxa_fabales <- data.frame(
  Taxa = data[[1]],  # 第一列是Taxa名称
  Fabales_freq = apply(data[, -1], 1, function(x) sum(x == "1Fabales", na.rm = TRUE)),
  stringsAsFactors = FALSE
)

# 获取前10个
top10_fabales <- taxa_fabales %>%
  arrange(desc(Fabales_freq)) %>%
  head(10)

# 4. 统计前10个最频繁鉴定为其他植物目的biomarker taxa
# 创建一个数据框来记录每个taxa在所有比较组中被鉴定为其他植物目的次数
taxa_other <- data.frame(
  Taxa = data[[1]],  # 第一列是Taxa名称
  Other_freq = apply(data[, -1], 1, function(x) {
    sum(x != "1Fabales" & !is.na(x) & x != "" & x != "-", na.rm = TRUE)
  }),
  stringsAsFactors = FALSE
)

# 获取前10个
top10_other <- taxa_other %>%
  arrange(desc(Other_freq)) %>%
  head(10)

# 6. 绘制柱形图
# color_manual <- c(rep(c("#7fc97f", "#beaed4", "#fdc086"), each = 10)) # 三界的颜色

# 为Fabales taxa绘图
p2 <- ggplot(top10_fabales, aes(x = reorder(Taxa, Fabales_freq), y = Fabales_freq)) +
  geom_bar(stat = "identity", fill = "#fdc086", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 10 Most Frequent Biomarker",
       x = "",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))

print(p2)

# 保存图形
ggsave("./plots/Top_10_Most_Frequent_Fabales_Biomarker_Protistan_Taxa.png", p2,
       width = 21, height = 7, dpi = 300, bg = "white", units = "cm", limitsize = FALSE)
ggsave("./plots/Top_10_Most_Frequent_Fabales_Biomarker_Protistan_Taxa.pdf", p2, 
       width = 21, height = 7, bg = "white", units = "cm", limitsize = FALSE)


# 为其他植物目taxa绘图
p3 <- ggplot(top10_other, aes(x = reorder(Taxa, Other_freq), y = Other_freq)) +
  geom_bar(stat = "identity", fill = "darkred", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 10 Most Frequent Biomarker",
       x = "",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))

print(p3)

# 保存图形
ggsave("./plots/Top_10_Most_Frequent_Fabales_Decreased_Protistan_Taxa.png", p3,
       width = 21, height = 7, dpi = 300, bg = "white", units = "cm", limitsize = FALSE)
ggsave("./plots/Top_10_Most_Frequent_Fabales_Decreased_Protistan_Taxa.pdf", p3, 
       width = 21, height = 7, bg = "white", units = "cm", limitsize = FALSE)





