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

###
core_taxonomy <- read.csv('../rawdata/amplicon/taxonomy/All_core_ASV_taxonomy.csv', row.names = 1)

lamda_file <- read.csv("../phylogenetic_signal/amplicon/core_asv/absolute_format_core_ASV_phylogenetic_signal_Fungi.csv", row.names = 1) %>% 
  rownames_to_column("FeatureID")


####统计韦恩图
order <- read.csv("../amplicons_randomforest/fun_qa_log_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
row.names(order) <- gsub("\\|.*", "", row.names(order))
order$ASVID <- row.names(order)
class <- read.csv("../amplicons_randomforest/fun_qa_log_Class_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
row.names(class) <- gsub("\\|.*", "", row.names(class))
class$ASVID <- row.names(class)
family <- read.csv("../amplicons_randomforest/fun_qa_log_Family_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
row.names(family) <- gsub("\\|.*", "", row.names(family))
family$ASVID <- row.names(family)
genus <- read.csv("../amplicons_randomforest/fun_qa_log_Genus_rf_res_diff_ASV.csv", sep = ",", header = T, row.names = 1)
row.names(genus) <- gsub("\\|.*", "", row.names(genus))
genus$ASVID <- row.names(genus)

lambda_rf_file <- lamda_file %>%
  left_join(core_taxonomy[,c("FeatureID", "ASVID", "ASVLabel", "Class", "Order", "Family", "Genus")], by = "FeatureID") %>%
  left_join(class[,c("ASVID", "Group")], by = "ASVID") %>%
  left_join(order[,c("ASVID", "Group")], by = "ASVID") %>%
  left_join(family[,c("ASVID", "Group")], by = "ASVID") %>%
  left_join(genus[,c("ASVID", "Group")], by = "ASVID") %>%
  column_to_rownames("ASVID")

names(lambda_rf_file)[13:16] <- c("Tree_Class", "Tree_Order", "Tree_Family", "Tree_Genus") 

# write.csv(lambda_rf_file, file = "fun_core_asv_qa_log_lambda_indicators.csv")


########
lamda_list <- row.names(lambda_rf_file[lambda_rf_file$FDR < 0.05,])

class_venn <- Venn(list('Phylogenetically conserved ASVs' = lamda_list, 
                        'Class indicators' = class$ASVID))

win.metafile("Fungal_class_venn.emf")
plot(class_venn, doWeights = TRUE, type="circles",  
     show = list(Faces = T, DarkMatter = FALSE))

dev.off()

order_venn <- Venn(list('Phylogenetically conserved ASVs' = lamda_list, 
                        'Order indicators' = order$ASVID))

win.metafile("Fungal_order_venn.emf")
plot(order_venn, doWeights = TRUE, type="circles",  
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()


family_venn <- Venn(list('Phylogenetically conserved ASVs' = lamda_list, 
                         'Family indicators' = family$ASVID))

win.metafile("Fungal_family_venn.emf")
plot(family_venn, doWeights = TRUE, type="circles", 
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

#############################################################################
####统计coevolved ASV贡献barplot
lambda_rf_file <- read.csv(file = "fun_core_asv_qa_log_lambda_indicators.csv")
table(subset(lambda_rf_file, FDR < 0.05 & Tree_Order != "NA")$Tree_Order)

lambda_rf_file$Taxa <- lambda_rf_file$Class

top_class_fun <- read.csv('../rawdata/amplicon/ITS_top10Class_color_meanValue.csv')

raw_df <- lambda_rf_file
raw_df$Taxa <- ifelse(raw_df$Taxa == "", "Unassigned", raw_df$Taxa)
raw_df$Taxa <- ifelse(raw_df$Taxa %in% top_class_fun$Taxa, raw_df$Taxa, 'Others')

all_df <- data.frame(t(table(raw_df$Taxa)))
lambda_df <- data.frame(t(table(raw_df[raw_df$FDR< 0.05, 'Taxa'])))
rf_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]), 'Taxa'])))
co_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]) & raw_df$FDR< 0.05, 'Taxa'])))
rf_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]), 'Taxa'])))
co_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]) & raw_df$FDR< 0.05, 'Taxa'])))
rf_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]), 'Taxa'])))
co_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]) & raw_df$FDR< 0.05, 'Taxa'])))

res_df <- all_df %>%
  left_join(lambda_df, by = "Var2") %>%
  left_join(rf_df_class, by = "Var2") %>%
  left_join(co_df_class, by = "Var2") %>%
  left_join(rf_df_order, by = "Var2") %>%
  left_join(co_df_order, by = "Var2") %>%
  left_join(rf_df_family, by = "Var2") %>%
  left_join(co_df_family, by = "Var2") 

res_df <- res_df[,c(2,3,5,7,9,11,13,15,17)]
res_df[is.na(res_df)] <- 0
names(res_df) <- c("Taxonomy", "All", "Phylogenetically conserved", 
                   "Class indicators", "Class co-evolved",
                   "Order indicators", "Order co-evolved",
                   "Family indicators", "Family co-evolved")

fin_df <- melt(res_df, id.vars = 'Taxonomy')

color_df <- data.frame(
  Taxonomy = top_class_fun$Taxa,
  Color = top_class_fun$Color
)
color_df <- color_df[color_df$Taxonomy %in% fin_df$Taxonomy,]
fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)
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
        # legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')

p <- ggplot(fin_df, aes(x = variable, y = value, fill = Taxonomy)) +
  geom_bar(stat = 'identity', position = position_fill(reverse = F), width = 0.68) +
  labs(title = '',
       subtitle = '',
       x = '',
       y = '') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_manual, name = "Taxa") +
  theme_bw() +  # 显示外边框
  theme_used
p
width = 5.5*2; height = 6.5
name = './Class distribution co-evolved Fungi ASVs'
ggsave(paste0(name, '.png'), p, width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
ggsave(paste0(name, '.pdf'), p, width = width, height = height, units = 'cm')

#############################################################
table(raw_df["Tree_Class"])
table(raw_df["Tree_Order"])
table(raw_df["Tree_Family"])
table(raw_df["Tree_Genus"])

class_used <- read.csv("../rawdata/amplicon/Tree_top_class_color.csv")
order_used <- read.csv("../rawdata/amplicon/Tree_top_order_color.csv")
family_used <- read.csv("../rawdata/amplicon/Tree_top_family_color.csv")
genus_used <- read.csv("../rawdata/amplicon/Tree_top_genus_color.csv")

##############################################不同 水平
# 定义处理函数
process_and_plot <- function(tax_level, level_order, file_suffix, facet_width_factor) {
  # 统一获取分类级别名称（带P值过滤）
  tax_levels <- names(table(
    raw_df[!is.na(raw_df[, tax_level]) & raw_df$FDR < 0.05, tax_level]
  ))
  
  # 初始化结果数据框
  names(lambda_df)[2] <- "Taxonomy"
  res_df <- lambda_df[-1]
  
  # 批量处理数据
  for (current_level in tax_levels) {
    # 创建指标数据（该分类下的所有taxa，无P值过滤）
    rf_data <- data.frame(
      Taxonomy = names(table(raw_df[raw_df[, tax_level] == current_level, 'Taxa'])),
      Freq_rf = as.numeric(table(raw_df[raw_df[, tax_level] == current_level, 'Taxa'])),
      stringsAsFactors = FALSE
    )
    names(rf_data)[2] <- paste(current_level, "indicators")
    
    # 创建共进化数据（该分类下且P值显著的taxa）
    co_data <- data.frame(
      Taxonomy = names(table(
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
      left_join(rf_data, by = "Taxonomy") %>% 
      left_join(co_data, by = "Taxonomy")
  }
  
  # 处理NA值
  res_df[is.na(res_df)] <- 0
  
  # 转换长格式
  fin_df <- melt(res_df[-2], id.vars = 'Taxonomy')
  
  # 提取绘图变量
  fin_df <- fin_df %>% mutate(
    Tree = factor(gsub(" .*", "", variable), levels = level_order),
    Biomarker = factor(gsub(".* ", "", variable), levels = c("indicators", "co-evolved"))
  )
  
  # 创建颜色映射
  color_df <- data.frame(
    Taxonomy = top_class_fun$Taxa,
    Color = top_class_fun$Color
  ) %>% filter(Taxonomy %in% fin_df$Taxonomy)
  
  fin_df$Taxonomy <- factor(fin_df$Taxonomy, levels = color_df$Taxonomy)
  
  # 根据分类层级调整字体大小
  strip_text_size <- if(tax_level %in% c("Tree_Family", "Tree_Genus")) 3 else 6
  
  # 生成图形
  p <- ggplot(fin_df, aes(x = Biomarker, y = value, fill = Taxonomy)) +
    geom_bar(stat = 'identity', width = 0.68) +
    labs(x = '', y = 'Number of ASVs', fill = "Taxa") +
    scale_fill_manual(values = color_df$Color) +
    facet_grid(~ Tree, scales = "free_x", space = "free", switch = "x") +
    theme_bw() + 
    theme_used +
    theme(strip.text.x = element_text(size = strip_text_size, color = "black")) +
    theme(legend.position = "none")
  
  # 保存图形
  width <- 5.5 * facet_width_factor
  height <- 6.5
  base_name <- paste0('./Class count co-evolved Fungal ASVs-', file_suffix)
  
  ggsave(paste0(base_name, '.png'), p, 
         width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
  ggsave(paste0(base_name, '.pdf'), p, 
         width = width, height = height, units = 'cm')
}

# 批量处理各分类水平
process_and_plot("Tree_Class", class_used$Class, "Class", 2)
process_and_plot("Tree_Order", order_used$Order, "Order", 4) 
process_and_plot("Tree_Family", family_used$Family, "Family", 4)



