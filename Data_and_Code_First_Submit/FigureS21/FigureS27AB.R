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
library(Vennerable)

### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)


core_pathway <- read.csv('../rawdata/KO/kegg_pathway_for_1KPM.csv', row.names = 1)
names(core_pathway)[1] <- "KOID"

core_annotation <- read.csv('../rawdata/KO/K_gene_name.csv', row.names = 2)
names(core_annotation)
core_annotation <- core_annotation[,-1]
core_annotation$KOLabel <- paste(row.names(core_annotation), core_annotation$gene, sep = "|")
core_annotation <- core_annotation[,c(3,1,2)]
core_annotation <- core_annotation %>% rownames_to_column("KOID")

bac_lamda_file <- read.csv("../phylogenetic_signal/metagenomic/Bacteria_KO_phylogenetic_signal_qa(log).csv", row.names = 1) %>% 
  rownames_to_column("KOID")


####统计韦恩图
class <- read.csv("../metagenomics_randomforest/bac_ko_qa_log_Class_tree_rf_res_diff_KO.csv", sep = ",", header = T, row.names = 1)
row.names(class) <- gsub("\\|.*", "", row.names(class))
class$KOID <- row.names(class)

order <- read.csv("../metagenomics_randomforest/bac_ko_qa_log_tree_rf_res_diff_KO.csv", sep = ",", header = T, row.names = 1)
row.names(order) <- gsub("\\|.*", "", row.names(order))
order$KOID <- row.names(order)

family <- read.csv("../metagenomics_randomforest/bac_ko_qa_log_Family_tree_rf_res_diff_KO.csv", sep = ",", header = T, row.names = 1)
row.names(family) <- gsub("\\|.*", "", row.names(family))
family$KOID <- row.names(family)



bac_lamda_rf_file <- bac_lamda_file %>%
  left_join(core_annotation, by = "KOID") %>%
  left_join(core_pathway, by = "KOID") %>%
  left_join(class[class$P.adj < 0.05, c("KOID", "Group")], by = "KOID") %>%
  left_join(order[order$P.adj < 0.05, c("KOID", "Group")], by = "KOID") %>%
  left_join(family[family$P.adj < 0.05, c("KOID", "Group")], by = "KOID") 

names(bac_lamda_rf_file)[17:19] <- c("Tree_Class", "Tree_Order", "Tree_Family") 

# write.csv(bac_lamda_rf_file, file = "bac_KO_qa_log_lamda_indicators.csv")

bac_lamda_ko_list <- unique(bac_lamda_rf_file[bac_lamda_rf_file$FDR < 0.05,]$KOID)


class_venn <- Venn(list('Phylogenetically conserved KOs' = bac_lamda_ko_list, 
                        'Class indicators' = subset(class, P.adj < 0.05)$KOID))

win.metafile("bac_KO_class_venn.emf")
plot(class_venn, doWeights = TRUE, type="circles",  
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

order_venn <- Venn(list('Phylogenetically conserved KOs' = bac_lamda_ko_list, 
                        'Order indicators' = subset(order, P.adj < 0.05)$KOID))

win.metafile("bac_KO_order_venn.emf")
plot(order_venn, doWeights = TRUE, type="circles",  
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

family_venn <- Venn(list('Phylogenetically conserved KOs' = bac_lamda_ko_list, 
                         'Family indicators' = subset(family, P.adj < 0.05)$KOID))

win.metafile("bac_KO_family_venn.emf")
plot(family_venn, doWeights = TRUE, type="circles", 
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

#####################################################################
####统计co-evolved KO 的superpathway贡献barplot

bac_lamda_rf_file <- read.csv(file = "bac_KO_qa_log_lamda_indicators.csv")
bac_lamda_rf_file$Taxa <- bac_lamda_rf_file$level2
table(bac_lamda_rf_file$level2) # superpathway level

ko_superpathway_color <- read.csv('../rawdata/KO/superpathway_color.csv')


raw_df <- bac_lamda_rf_file[, c("KOID","lambda", "P", "FDR", "Taxa", "Tree_Class", "Tree_Order", "Tree_Family")] ###原文件是在KOID上合并的，在superpathway上有重复，如同一个KO多次在一个pathway上出现
raw_df <- unique.data.frame(raw_df) #去重复
raw_df <- raw_df[!is.na(raw_df$Taxa),] # 7920个

all_df <- data.frame(t(table(raw_df$Taxa)))
bac_lamda_df <- data.frame(t(table(raw_df[raw_df$FDR < 0.05, 'Taxa'])))
bac_rf_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]), 'Taxa'])))
co_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]) & raw_df$FDR < 0.05, 'Taxa'])))
bac_rf_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]), 'Taxa'])))
co_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]) & raw_df$FDR < 0.05, 'Taxa'])))
bac_rf_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]), 'Taxa'])))
co_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]) & raw_df$FDR < 0.05, 'Taxa'])))

res_df <- all_df %>%
  left_join(bac_lamda_df, by = "Var2") %>%
  left_join(bac_rf_df_class, by = "Var2") %>%
  left_join(co_df_class, by = "Var2") %>%
  left_join(bac_rf_df_order, by = "Var2") %>%
  left_join(co_df_order, by = "Var2") %>%
  left_join(bac_rf_df_family, by = "Var2") %>%
  left_join(co_df_family, by = "Var2") 

res_df <- res_df[,c(2,3,5,7,9,11,13,15,17)]
res_df[is.na(res_df)] <- 0
names(res_df) <- c("Superpathway", "All", "Phylogenetically conserved", 
                   "Class indicators", "Class co-evolved",
                   "Order indicators", "Order co-evolved",
                   "Family indicators", "Family co-evolved")

fin_df <- melt(res_df, id.vars = 'Superpathway')

color_df <- data.frame(
  Superpathway = ko_superpathway_color$superpathway,
  Color = ko_superpathway_color$fill_color2
)
color_df <- color_df[color_df$Superpathway %in% fin_df$Superpathway,]
fin_df$Superpathway <- factor(fin_df$Superpathway, levels = color_df$Superpathway)
color_manual <- color_df$Color

theme_used <-
  theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 7, color = 'black'),
        axis.line = element_blank(),
        axis.text = element_text(size = 6, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 7, color = 'black'),
        legend.text = element_text(size = 6,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')

p <- ggplot(fin_df, aes(x = variable, y = value, fill = Superpathway)) +
  geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
  labs(title = '',
       subtitle = '',
       x = '',
       y = '') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_manual, name = "Superpathway") +
  theme_bw() +  # 显示外边框
  theme_used +
  guides(fill = guide_legend( ncol = 1, byrow = TRUE)) 

p
width = 5.5*2; height = 8.5
name = './Superpathway distribution bacterial co-evolved KOs'
ggsave(paste0(name, '.png'), p, width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
ggsave(paste0(name, '.pdf'), p, width = width, height = height, units = 'cm')

#############################################################
table(raw_df["Tree_Class"])
table(raw_df["Tree_Order"])
table(raw_df["Tree_Family"])

class_used <- read.csv("../rawdata/amplicon/Tree_top_class_color.csv")
order_used <- read.csv("../rawdata/amplicon/Tree_top_order_color.csv")
family_used <- read.csv("../rawdata/amplicon/Tree_top_family_color.csv")

##############################################不同 水平
# 定义处理函数
process_and_plot <- function(tax_level, level_order, file_suffix, facet_width_factor) {
  # 统一获取分类级别名称（带P值过滤）
  tax_levels <- names(table(
    raw_df[!is.na(raw_df[, tax_level]) & raw_df$FDR < 0.05, tax_level]
  ))
  
  # 初始化结果数据框
  names(bac_lamda_df)[2] <- "Superpathway"
  res_df <- bac_lamda_df[-1]
  
  # 批量处理数据
  for (current_level in tax_levels) {
    # 创建指标数据（该分类下的所有taxa，无P值过滤）
    bac_rf_data <- data.frame(
      Superpathway = names(table(raw_df[raw_df[, tax_level] == current_level, 'Taxa'])),
      Freq_rf = as.numeric(table(raw_df[raw_df[, tax_level] == current_level, 'Taxa'])),
      stringsAsFactors = FALSE
    )
    names(bac_rf_data)[2] <- paste(current_level, "indicators")
    
    # 创建共进化数据（该分类下且P值显著的taxa）
    co_data <- data.frame(
      Superpathway = names(table(
        raw_df[raw_df[, tax_level] == current_level & raw_df$FDR < 0.05, 'Taxa']
      )),
      Freq_co = as.numeric(table(
        raw_df[raw_df[, tax_level] == current_level & raw_df$FDR < 0.05, 'Taxa']
      )),
      stringsAsFactors = FALSE
    )
    names(co_data)[2] <- paste(current_level, "co-evolved")
    
    # 合并数据
    res_df <- res_df %>% 
      left_join(bac_rf_data, by = "Superpathway") %>% 
      left_join(co_data, by = "Superpathway")
  }
  
  # 处理NA值
  res_df[is.na(res_df)] <- 0
  
  # 转换长格式
  fin_df <- melt(res_df[-2], id.vars = 'Superpathway')
  
  # 提取绘图变量
  fin_df <- fin_df %>% mutate(
    Tree = factor(gsub(" .*", "", variable), levels = level_order),
    Biomarker = factor(gsub(".* ", "", variable), levels = c("indicators", "co-evolved"))
  )
  
  # 创建颜色映射
  color_df <- data.frame(
    Superpathway = ko_superpathway_color$superpathway,
    Color = ko_superpathway_color$fill_color2
  ) %>% filter(Superpathway %in% fin_df$Superpathway)
  
  fin_df$Superpathway <- factor(fin_df$Superpathway, levels = color_df$Superpathway)
  
  # 根据分类层级调整字体大小
  strip_text_size <- if(tax_level %in% c("Tree_Family")) 3 else 6
  
  # 生成图形
  p <- ggplot(fin_df, aes(x = Biomarker, y = value, fill = Superpathway)) +
    geom_bar(stat = 'identity', width = 0.68) +
    labs(x = '', y = 'Number of KOs', fill = "Superpathway") +
    scale_fill_manual(values = color_df$Color) +
    facet_grid(~ Tree, scales = "free_x", space = "free", switch = "x") +
    theme_bw() + 
    theme_used +
    theme(strip.text.x = element_text(size = strip_text_size, color = "black"))+
    guides(fill = guide_legend( ncol = 1, byrow = TRUE)) +
    theme(legend.position = "none")
  
  # 保存图形
  width <- 5.5 * facet_width_factor
  height <- 8.5
  base_name <- paste0('./Superpathway count bacterial co-evolved KOs-', file_suffix)
  
  ggsave(paste0(base_name, '.png'), p, 
         width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
  ggsave(paste0(base_name, '.pdf'), p, 
         width = width, height = height, units = 'cm')
}

# 批量处理各分类水平
process_and_plot("Tree_Class", class_used$Class, "Class", 2)
process_and_plot("Tree_Order", order_used$Order, "Order", 4) 
process_and_plot("Tree_Family", family_used$Family, "Family", 4)
