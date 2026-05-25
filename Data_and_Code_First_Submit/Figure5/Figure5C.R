### 加载必要包 ----------------------------------------------------------------
library(dplyr)
library(tibble)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(RColorBrewer)
library(ggnewscale)
library(reshape2)
library(scales)

### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

### 数据准备 ------------------------------------------------------------------
# 导入数据
core_taxonomy <- read.csv('./data/All_core_ASV_taxonomy.csv', row.names = 1)

lamda_file <- read.csv("./data/absolute_format_core_ASV_phylogenetic_signal_protist.csv", row.names = 1) %>% 
  rownames_to_column("FeatureID")

tree <- treeio::read.newick("./data/pro_core_feature_tree.nwk", node.label = "support")

rf_file <-  read.csv("./data/pro_qa_log_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
row.names(rf_file) <- gsub("\\|.*", "", row.names(rf_file))
rf_file$ASVID <- row.names(rf_file) 
table(rf_file$P.adj < 0.05)

# 颜色配置
class_colors <- read.csv("./data/Protist_top10Class_color_meanValue.csv")
plant_colors <- read.csv("./data/Tree_top10Order_color.csv")

# 数据合并
group_file <- core_taxonomy %>%
  left_join(lamda_file, by = "FeatureID") %>%
  left_join(rf_file, by = "ASVID") %>%
  filter(FeatureID %in% tree@phylo$tip.label) %>%
  column_to_rownames("FeatureID")

# 分类处理
group_file <- group_file %>%
  mutate(
    # 处理分类合并
    Class = case_when(
      Class %in% class_colors$Taxa[1:10] ~ Class,       # 保留前10个指定分类
      Class == "" ~ "Unassigned",                       # 处理空字符串为Unassigned
      TRUE ~ "Others"                                   # 其他情况归为Others
    ),
    
    # 创建因子水平
    Class = factor(Class,levels = class_colors$Taxa),
    
    # 创建系统发育信号标签
    Group = if_else(
      FDR < 0.05 & P.adj < 0.05, Group, #lambda显著，且为rf鉴定biomarker 
      ""
    )
  )

table(group_file$Group)
table(group_file$FDR < 0.05)
table(group_file$P.adj < 0.05)

### 可视化设置 ----------------------------------------------------------------

### 数据整合与分组 -----------------------------------------------------------
# 创建分组列表（关键步骤）
groupInfo <- split(
  rownames(group_file), 
  group_file$Class  # 使用微生物Class作为分组依据
)

# 使用groupOTU整合分组信息
tree_grouped <- groupOTU(
  tree,
  groupInfo,
  group_name = "Taxonomy"  # 自定义分组变量名
)

### 可视化构建 --------------------------------------------------------------

### 显示lamda结果，颜色深浅代表系统发育信号强弱
lambda_plot <- group_file %>%
  rownames_to_column("label") %>%
  select(label, lambda) 

### 显示randomforest&lambda结果, 是某个植物的就绘制一个圆点
rf_plot <- group_file %>%
  rownames_to_column("label") %>%
  select(label, Group) %>%
  # 过滤无效值，并确保Group为因子
  filter(!is.na(Group) & Group != "") %>%
  mutate(Group = factor(Group))




### 完整可视化流程 ----------------------------------------------------------
# 基础树结构（包含分组信息）
base_tree <- ggtree(tree_grouped, layout = 'fan', right = TRUE, linewidth = 0.15, open.angle = 15) +  
  geom_aline(aes(color = Taxonomy), linetype = 'solid', size = 0.2, linewidth = 0.15, key_glyph = 'rect', show.legend = TRUE) +
  scale_color_manual(values = class_colors$Color, breaks = class_colors$Taxa) + 
  new_scale_color()

base_tree

# 最终绘图
final_plot <- base_tree +
  
  # 系统发育信号层
  geom_fruit(
    data = lambda_plot,
    geom = geom_tile,
    mapping = aes(y=label, x = "Lambda", fill= lambda),
    color = NA,  # 去除边框
    width = 0.12,
    offset = 0.05, #调整offset
    axis.params = list(
      axis="x",
      line.color="white",
      text.size = 0.6,
      nbreak = 2, 
      text.angle = 85, 
      vjust = 0, 
      hjust = 1)
  ) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "red"))(1000), 
                       name = "Lambda") +
  new_scale_fill()  +
  
  geom_fruit(
    data = rf_plot,
    geom = geom_point, 
    mapping = aes(y = label, x = "Coevolved", color = Group), 
    size = 0.01,
    pwidth = 0.1,
    offset = 0.05, #调整offset
    axis.params = list(
      axis="x",
      line.color="white",
      text.size = 0.6,
      nbreak = 2, 
      text.angle = 85, 
      vjust = 0, 
      hjust = 1),
    grid.params = list(vline = TRUE)
  ) +
  scale_color_manual(
    values = plant_colors$Color2, 
    breaks = plant_colors$Order,
    name = "Coevolved"
  ) + 
  new_scale_color() 


#########################tree all core ASV
ggsave(final_plot,
       filename = "./result/Figure5C_protistan_core_ASVs_lambda_rf.pdf", width = 12, height = 12,  dpi = 900, bg = "white")

###
final_plot2 <- final_plot + theme(legend.position = "none") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggsave(final_plot2,
       filename = "./result/Figure5C_protistan_core_ASVs_lambda_rf_nolegend.pdf",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot2,
       filename = "./result/Figure5C_protistan_core_ASVs_lambda_rf_nolegend.tiff",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

ggsave(final_plot2,
       filename = "./result/Figure5C_protistan_core_ASVs_lambda_rf_nolegend.jpg",
       width = 115/2, height = 115/2, units = "mm", dpi = 900, bg = "white")

