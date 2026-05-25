### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

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
library("clusterProfiler")

####数据准备####
plant_colors <- read.csv("../../../../kegg_abundance/Tree_top_order_color.csv")
core_annotation <- read.csv('../../../../kegg_abundance/K_gene_name.csv', row.names = 2)
names(core_annotation)
core_annotation <- core_annotation[,-1]
core_annotation$KOLabel <- paste(row.names(core_annotation), core_annotation$gene, sep = "|")
core_annotation <- core_annotation[,c(3,1,2)]
core_annotation <- core_annotation %>% rownames_to_column("KOID")

core_pathway <- read.csv('../../../../kegg_abundance/kegg_pathway_for_1KPM.csv', row.names = 1)
names(core_pathway)[1] <- "KOID"

########################################################################
metadata_tree <- read.csv(file = "../../../../metadata/tree_metadata_merge_info.csv", row.names = 1)
metadata_tree$Order <- ifelse(metadata_tree$Order %in% plant_colors$Order, metadata_tree$Order, "Others")

pro_ko_ra <- read.csv(file = "../../../../kegg_abundance/relative/Protists_core0.2_tpm.csv", row.names = 1, header = T)
names(pro_ko_ra) <- gsub("RS.", "", names(pro_ko_ra))
pro_ko_qa <- read.csv(file = "../../../../kegg_abundance/absolute/Protists_core0.2_quantitative.csv", row.names = 1, header = T)
names(pro_ko_qa) <- gsub("RS.", "", names(pro_ko_qa))


pro_lamda_file <- read.csv("../../../../12.phylogenetic_signal_λ/metagenomic/core_KO/Protists_core0.2_qa_log_phylogenetic_signal.csv", row.names = 1) %>% 
  rownames_to_column("KOID")
table(pro_lamda_file$padj < 0.05) ###160个系统发图显著的KO

pro_rf_file <-  read.csv("../../../../13.random_forest/metagenome/pro_ko_core0.2_qa_log_tree_rf_res_diff_KO.csv", sep = ",", header = T, row.names = 1)
# 该文件只包含显著性矫正P值小于0.2的KO
row.names(pro_rf_file) <- gsub("\\|.*", "", row.names(pro_rf_file))
pro_rf_file$KOID <- row.names(pro_rf_file) 
table(pro_rf_file$P.adj < 0.05) ###198个随机森林显著的KO 

# 数据合并
pro_group_file <- pro_lamda_file %>%
  left_join(pro_rf_file, by = "KOID") %>%
  left_join(core_annotation, by = "KOID") 

table(pro_group_file$Group)
table(pro_group_file$padj < 0.05)  ###160 ko
table(pro_group_file$P < 0.05)
table(pro_group_file$P.adj < 0.05) ###198 ko

######################################################################
######################################################################

#####################################################################################################
####统计韦恩图
pro_ko_list <- pro_group_file$KOID
pro_lamda_ko_list <- subset(pro_group_file, padj < 0.05)$KOID
pro_rf_ko_list <- subset(pro_group_file, P.adj < 0.05)$KOID
order_co_ko_list <- subset(pro_group_file, padj < 0.05 & P.adj < 0.05)$KOID #100


order_venn <- Venn(list('Phylogenetically conserved KOs' = pro_lamda_ko_list, 
                        'Order indicators' = pro_rf_ko_list))

win.metafile("./plots/pro_KO_order_venn.emf")
plot(order_venn, doWeights = TRUE, type="circles", 
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

###Class和Family水平
class <- read.csv("../../../../13.random_forest/metagenome/pro_ko_qa_log_Class_tree_rf_res_diff_KO.csv", sep = ",", header = T, row.names = 1)
row.names(class) <- gsub("\\|.*", "", row.names(class))
class$KOID <- row.names(class)
family <- read.csv("../../../../13.random_forest/metagenome/pro_ko_qa_log_Family_tree_rf_res_diff_KO.csv", sep = ",", header = T, row.names = 1)
row.names(family) <- gsub("\\|.*", "", row.names(family))
family$KOID <- row.names(family)

class_venn <- Venn(list('Phylogenetically conserved KOs' = pro_lamda_ko_list, 
                        'Class indicators' = subset(class, P.adj < 0.05)$KOID))

win.metafile("./plots/pro_KO_class_venn.emf")
plot(class_venn, doWeights = TRUE, type="circles",  
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

family_venn <- Venn(list('Phylogenetically conserved KOs' = pro_lamda_ko_list, 
                         'Family indicators' = subset(family, P.adj < 0.05)$KOID))

win.metafile("./plots/pro_KO_family_venn.emf")
plot(family_venn, doWeights = TRUE, type="circles", 
     show = list(Faces = T, DarkMatter = FALSE))
dev.off()

#####################################################################
###superpathway 贡献barplot
#############################################################################
####统计co-evoluted KO贡献barplot

pro_lamda_rf_file <- pro_lamda_file %>%
  left_join(core_annotation, by = "KOID") %>%
  left_join(core_pathway[,-8], by = "KOID") %>%
  left_join(class[class$P.adj < 0.05, c("KOID", "Group","P.adj")], by = "KOID") %>%
  left_join(pro_rf_file[pro_rf_file$P.adj < 0.05, c("KOID", "Group", "P.adj")], by = "KOID") %>%  #pro_rf_file就是order水平
  left_join(family[family$P.adj < 0.05, c("KOID", "Group","P.adj")], by = "KOID") 

names(pro_lamda_rf_file)[17:22] <- c("Tree_Class", "rf_FDR_Class", "Tree_Order", "rf_FDR_Order", "Tree_Family", "rf_FDR_Family") 

write.csv(pro_lamda_rf_file, file = "./stats/pro_KO_qa_log_lamda_indicators.csv")
# write.csv(unique(pro_lamda_rf_file[,-c(11:16)]), file = "./stats/TableS7F-pro_KO_qa_log_lamda_indicators.csv")

table(pro_lamda_rf_file[!is.na(pro_lamda_rf_file$Tree_Order) & pro_lamda_rf_file$padj < 0.05,]$Tree_Order)
#########

pro_lamda_rf_file <- read.csv(file = "./stats/pro_KO_qa_log_lamda_indicators.csv")
pro_lamda_rf_file$Superpathway <- pro_lamda_rf_file$level2
table(pro_lamda_rf_file$level2) # superpathway level

ko_superpathway_color <- read.csv('../../../../metadata/superpathway_color.csv')


raw_df <- pro_lamda_rf_file[, c("KOID","lambda", "P", "padj", "Superpathway", "Tree_Class", "Tree_Order", "Tree_Family")] 
###原文件是在KOID上合并的，在superpathway上有重复，如同一个KO多次在一个pathway上出现
raw_df <- unique.data.frame(raw_df) #去重复
raw_df <- raw_df[!is.na(raw_df$Superpathway),] # 2318个

all_df <- data.frame(t(table(raw_df$Superpathway)))
pro_lamda_df <- data.frame(t(table(raw_df[raw_df$padj < 0.05, 'Superpathway'])))
pro_rf_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]), 'Superpathway'])))
co_df_class <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Class"]) & raw_df$padj < 0.05, 'Superpathway']))) # 没有显著的
co_df_class <- pro_rf_df_class; co_df_class$Freq <- 0 # 没有显著的, 用pro_rf_df_class替代，但是值为0
pro_rf_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]), 'Superpathway'])))
co_df_order <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Order"]) & raw_df$padj < 0.05, 'Superpathway'])))
pro_rf_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]), 'Superpathway'])))
co_df_family <- data.frame(t(table(raw_df[!is.na(raw_df["Tree_Family"]) & raw_df$padj < 0.05, 'Superpathway'])))

res_df <- all_df %>%
  left_join(pro_lamda_df, by = "Var2") %>%
  left_join(pro_rf_df_class, by = "Var2") %>%
  left_join(co_df_class, by = "Var2") %>%   
  left_join(pro_rf_df_order, by = "Var2") %>%
  left_join(co_df_order, by = "Var2") %>%
  left_join(pro_rf_df_family, by = "Var2") %>%
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
name = './plots/Superpathway_distribution_protistan_co-evolved KOs'
ggsave(paste0(name, '.png'), p, width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
ggsave(paste0(name, '.pdf'), p, width = width, height = height, units = 'cm')

#############################################################
table(raw_df["Tree_Class"])
table(raw_df["Tree_Order"])
table(raw_df["Tree_Family"])

class_used <- read.csv("../../../../metadata/Tree_top_class_color.csv")
order_used <- read.csv("../../../../metadata/Tree_top_order_color.csv")
family_used <- read.csv("../../../../metadata/Tree_top_family_color.csv")

##############################################不同 水平
# 定义处理函数
process_and_plot <- function(tax_level, level_order, file_suffix, facet_width_factor) {
  # 统一获取分类级别名称（带P值过滤）
  tax_levels <- names(table(
    raw_df[!is.na(raw_df[, tax_level]) & raw_df$padj < 0.05, tax_level]
  ))
  
  # 初始化结果数据框
  names(pro_lamda_df)[2] <- "Superpathway"
  res_df <- pro_lamda_df[-1]
  
  # 批量处理数据
  for (current_level in tax_levels) {
    # 创建指标数据（该分类下的所有Superpathway，无P值过滤）
    pro_rf_data <- data.frame(
      Superpathway = names(table(raw_df[raw_df[, tax_level] == current_level, 'Superpathway'])),
      Freq_rf = as.numeric(table(raw_df[raw_df[, tax_level] == current_level, 'Superpathway'])),
      stringsAsFactors = FALSE
    )
    names(pro_rf_data)[2] <- paste(current_level, "indicators")
    
    # 创建共进化数据（该分类下且P值显著的Superpathway）
    co_data <- data.frame(
      Superpathway = names(table(
        raw_df[raw_df[, tax_level] == current_level & raw_df$padj < 0.05, 'Superpathway']
      )),
      Freq_co = as.numeric(table(
        raw_df[raw_df[, tax_level] == current_level & raw_df$padj < 0.05, 'Superpathway']
      )),
      stringsAsFactors = FALSE
    )
    names(co_data)[2] <- paste(current_level, "co-evolved")
    
    # 合并数据
    res_df <- res_df %>% 
      left_join(pro_rf_data, by = "Superpathway") %>% 
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
    geom_bar(stat = 'identity', position = "stack", width = 0.68) +     #position = c("stack", "fill")
    labs(x = '', y = 'Number of KOs', fill = "Superpathway") +     # Number/Percentage
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
  base_name <- paste0('./plots/Superpathway_count_protistan_co-evolved-KOs-Number-', file_suffix)     # Number/Percentage
  
  ggsave(paste0(base_name, '.png'), p, 
         width = width, height = height, dpi = 600, type = 'cairo', units = 'cm')
  ggsave(paste0(base_name, '.pdf'), p, 
         width = width, height = height, units = 'cm')
}


# 批量处理各分类水平
process_and_plot("Tree_Class", class_used$Class, "Class", 2)
process_and_plot("Tree_Order", order_used$Order, "Order", 4) 
process_and_plot("Tree_Family", family_used$Family, "Family", 4)

##################################################################################################
#####################富集分析

# （1）所有具备系统发育信号的KO
# （2）所有rf鉴定为top植物纲目科biomarker的KO
# （3）具备系统发育信号的biomarker （co-evolved KO）

pro_lamda_rf_file <- read.csv(file = "./stats/pro_KO_qa_log_lamda_indicators.csv")
names(pro_lamda_rf_file)
pro_lamda_rf_enrich <- unique(pro_lamda_rf_file[,-c(1,11:17)]) ###2458行，2458个KO， 没错的

table(subset(pro_lamda_rf_enrich, padj < 0.05 & Tree_Order != "NA")$Tree_Order)

pro_core_ko_list <- pro_lamda_rf_enrich$KOID
pro_lamda_ko_list <- subset(pro_lamda_rf_enrich, padj < 0.05)$KOID
class_ko_list <- subset(pro_lamda_rf_enrich, Tree_Class != "NA")$KOID
class_co_ko_list <- subset(pro_lamda_rf_enrich, padj < 0.05 & Tree_Class != "NA")$KOID
order_ko_list <- subset(pro_lamda_rf_enrich, Tree_Order != "NA")$KOID
order_co_ko_list <- subset(pro_lamda_rf_enrich, padj < 0.05 & Tree_Order != "NA")$KOID
family_ko_list <- subset(pro_lamda_rf_enrich, Tree_Family != "NA")$KOID
family_co_ko_list <- subset(pro_lamda_rf_enrich, padj < 0.05 & Tree_Family != "NA")$KOID

table(pro_lamda_rf_enrich$Tree_Class)

#进行富集分析
pro_core_ko_result=enrichKEGG(pro_core_ko_list,
                          organism = "ko", #物种选择ko
                          pvalueCutoff = 0.05,  #p值cutoff
                          pAdjustMethod = "fdr", #fdr矫正p值
                          qvalueCutoff = 0.2   #q值cutoff
)

pro_lamda_ko_result=enrichKEGG(pro_lamda_ko_list,
                            organism = "ko", #物种选择ko
                            pvalueCutoff = 0.05,  #p值cutoff
                            pAdjustMethod = "fdr", #fdr矫正p值
                            qvalueCutoff = 0.2   #q值cutoff
)

class_ko_result=enrichKEGG(class_ko_list,
                           organism = "ko", #物种选择ko
                           pvalueCutoff = 0.05,  #p值cutoff
                           pAdjustMethod = "fdr", #fdr矫正p值
                           qvalueCutoff = 0.2   #q值cutoff
)
class_ko_result@result$qvalue <- 1

class_co_ko_result=enrichKEGG(class_co_ko_list,
                              organism = "ko", #物种选择ko
                              pvalueCutoff = 0.05,  #p值cutoff
                              pAdjustMethod = "fdr", #fdr矫正p值
                              qvalueCutoff = 0.2   #q值cutoff
)

# class_co_ko_list内容为空，没有富集结果，使用class_ko_result替代，但是pvalue改为1，count改为0
class_co_ko_result <- class_ko_result
class_co_ko_result@result$pvalue <- 1; class_co_ko_result@result$p.adjust <- 1; class_co_ko_result@result$Count <- 0

order_ko_result=enrichKEGG(order_ko_list,
                              organism = "ko", #物种选择ko
                              pvalueCutoff = 0.05,  #p值cutoff
                              pAdjustMethod = "fdr", #fdr矫正p值
                              qvalueCutoff = 0.2   #q值cutoff
)

order_co_ko_result=enrichKEGG(order_co_ko_list,
                              organism = "ko", #物种选择ko
                              pvalueCutoff = 0.05,  #p值cutoff
                              pAdjustMethod = "fdr", #fdr矫正p值
                              qvalueCutoff = 0.2   #q值cutoff
)

family_ko_result=enrichKEGG(family_ko_list,
                           organism = "ko", #物种选择ko
                           pvalueCutoff = 0.05,  #p值cutoff
                           pAdjustMethod = "fdr", #fdr矫正p值
                           qvalueCutoff = 0.2   #q值cutoff
)

family_co_ko_result=enrichKEGG(family_co_ko_list,
                              organism = "ko", #物种选择ko
                              pvalueCutoff = 0.05,  #p值cutoff
                              pAdjustMethod = "fdr", #fdr矫正p值
                              qvalueCutoff = 0.2   #q值cutoff
)

#保存富集分析结果
write.csv(file="./stats/KO_enrichment_all_pro_core_ko.csv", data.frame(pro_core_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_pro_lamda_ko.csv", data.frame(pro_lamda_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_plant_class_biomarker_ko.csv", data.frame(class_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_plant_class_coevolved_ko.csv", data.frame(class_co_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_plant_order_biomarker_ko.csv", data.frame(order_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_plant_order_coevolved_ko.csv", data.frame(order_co_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_plant_family_biomarker_ko.csv", data.frame(family_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_all_plant_family_coevolved_ko.csv", data.frame(family_co_ko_result), row.names=F)

# 绘制气泡图
barplot(pro_core_ko_result)
barplot(pro_lamda_ko_result)
barplot(order_ko_result)
barplot(order_co_ko_result)

####整合图
inter_enriched <- rbind.data.frame(data.frame(pro_core_ko_result),
                                   data.frame(pro_lamda_ko_result),
                                   data.frame(class_ko_result), 
                                   data.frame(class_co_ko_result),
                                   data.frame(order_ko_result),
                                   data.frame(order_co_ko_result),
                                   data.frame(family_ko_result), 
                                   data.frame(family_co_ko_result))
  
inter_enriched$Group <- c(rep("All", dim(data.frame(pro_core_ko_result))[1]), 
                          rep("Phylogenetically conserved", dim(data.frame(pro_lamda_ko_result))[1]),
                          rep("Class indicators", dim(data.frame(class_ko_result))[1]),
                          rep("Class co-evolved", dim(data.frame(class_co_ko_result))[1]),
                          rep("Order indicators", dim(data.frame(order_ko_result))[1]),
                          rep("Order co-evolved", dim(data.frame(order_co_ko_result))[1]),
                          rep("Family indicators", dim(data.frame(family_ko_result))[1]),
                          rep("Family co-evolved", dim(data.frame(family_co_ko_result))[1]))

inter_enriched$Group <- factor(inter_enriched$Group, levels = c("All", "Phylogenetically conserved", 
                                                                "Class indicators", "Class co-evolved",
                                                                "Order indicators", "Order co-evolved",
                                                                "Family indicators", "Family co-evolved"))

names(inter_enriched)
inter_enriched$qvalue

inter_enriched_g <- ggplot(inter_enriched, aes(x=-log2(qvalue), y=Description)) +
  geom_point(aes(size=log(Count,base=2), color=-log2(qvalue))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("white", "red"))(1000),
    values = rescale(c(min(-log2(inter_enriched$qvalue)),
                       median(-log2(inter_enriched$qvalue)),
                       max(-log2(inter_enriched$qvalue))))) +
  scale_size_continuous(name="log2(gene counts)") +
  labs(title="", x="-log2(qvalue)", y="") +
  facet_grid(~ Group) +
  theme_bw() +
  theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 7, color = 'black'),
        axis.line = element_blank(),
        axis.text = element_text(size = 6, color = 'black'),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 7, color = 'black'),
        legend.text = element_text(size = 6,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')

windows(12,12)
inter_enriched_g

ggsave(inter_enriched_g,
       filename = "./plots/KO_enrichment_all_plant_class_order_family_protistan_coevolved_ko.pdf",
       width = 20, height = 32, units = "cm", dpi = 900, bg = "white")
ggsave(inter_enriched_g,
       filename = "./plots/KO_enrichment_all_plant_class_order_family_protistan_coevolved_ko.png",
       width = 20, height = 32, units = "cm", dpi = 900, bg = "white")



####只富集"All"和 "Phylogenetically conserved"
inter_enriched2 <- subset(inter_enriched, Group %in% c("All","Phylogenetically conserved") & qvalue < 0.001)
table(table(inter_enriched2$Description) == 2)
table(inter_enriched2$Group)


inter_enriched_g2 <- ggplot(inter_enriched2, aes(x=-log2(qvalue), y=Description)) +
  geom_point(aes(size=log(Count,base=2), color=-log2(qvalue))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("white", "red"))(1000),
    values = rescale(c(min(-log2(inter_enriched$qvalue)),
                       median(-log2(inter_enriched$qvalue)),
                       max(-log2(inter_enriched$qvalue))))) +
  scale_size_continuous(name="log2(gene counts)") +
  labs(title="", x="-log2(qvalue)", y="") +
  facet_grid(~ Group) +
  theme_bw() +
  theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 7, color = 'black'),
        axis.line = element_blank(),
        axis.text = element_text(size = 6, color = 'black'),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 7, color = 'black'),
        legend.text = element_text(size = 6,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')

windows(8,12)
inter_enriched_g2

ggsave(inter_enriched_g2,
       filename = "./plots/KO_enrichment_pro_lamda_and_Phylogenetically_conserved_ko.pdf",
       width = 20, height = 32, units = "cm", dpi = 900, bg = "white")
ggsave(inter_enriched_g2,
       filename = "./plots/KO_enrichment_pro_lamda_and_Phylogenetically_conserved_ko.png",
       width = 20, height = 32, units = "cm", dpi = 900, bg = "white")




##################################################################################
###不同目水平的indicators 和co-evoled富集分析
table(pro_lamda_rf_enrich$Tree_Order)

# 豆目没有检测到共进化的真菌KO，使用Asparagales
Asparagales_ko_list <- subset(pro_lamda_rf_enrich, Tree_Order == "Asparagales")$KOID
Asparagales_co_ko_list <- subset(pro_lamda_rf_enrich, Tree_Order == "Asparagales" & padj < 0.05)$KOID
#进行富集分析
Asparagales_ko_result=enrichKEGG(Asparagales_ko_list,
                                organism = "ko", #物种选择ko
                                pvalueCutoff = 0.05,  #p值cutoff
                                pAdjustMethod = "fdr", #fdr矫正p值
                                qvalueCutoff = 0.2   #q值cutoff
)

Asparagales_co_ko_result=enrichKEGG(Asparagales_co_ko_list,
                                organism = "ko", #物种选择ko
                                pvalueCutoff = 0.05,  #p值cutoff
                                pAdjustMethod = "fdr", #fdr矫正p值
                                qvalueCutoff = 0.2   #q值cutoff
)

barplot(Asparagales_ko_result)
barplot(Asparagales_co_ko_result)

#保存富集分析结果
write.csv(file="./stats/KO_enrichment_Asparagales_pro_coevolved_ko.csv", data.frame(Asparagales_co_ko_result), row.names=F)
write.csv(file="./stats/KO_enrichment_Asparagales_pro_rf_ko.csv", data.frame(Asparagales_ko_result), row.names=F)

###循环top10 所有植物目
# 创建目录保存结果
dir.create("enrichment_results", showWarnings = FALSE)
dir.create("enrichment_plots", showWarnings = FALSE)

# 所有植物目列表
orders <- c("Arecales", "Asparagales", "Fabales", "Gentianales", 
            "Lamiales", "Malpighiales", "Malvales", "Myrtales", 
            "Others", "Rosales", "Sapindales")

# 创建空数据框收集所有结果
order_results <- data.frame()

# 循环处理每个植物目
for (order in orders) {
  cat("Processing order:", order, "\n")
  
  # 提取当前目的基因列表
  order_ko_list <- subset(pro_lamda_rf_enrich, Tree_Order == order)$KOID
  order_co_ko_list <- subset(pro_lamda_rf_enrich, Tree_Order == order & padj < 0.05)$KOID
  
  # 准备存储结果的列表
  results_list <- list()
  
  # 1. 分析指示基因 (indicators)
  if (length(order_ko_list) > 0) {
    ko_result <- tryCatch({
      enrichKEGG(
        gene = order_ko_list,
        organism = "ko",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("Error in enrichKEGG for ", order, " indicators: ", e$message)
      NULL
    })
    
    if (!is.null(ko_result) && !is.null(ko_result@result)) {
      ko_df <- as.data.frame(ko_result)
      
      # 仅在有结果时添加标识列
      if (nrow(ko_df) > 0) {
        ko_df$Category <- paste0(order, "_indicators")
        results_list[["indicators"]] <- ko_df
        
        # 保存图表
        png_file <- file.path("enrichment_plots", paste0(order, "_indicators.png"))
        png(png_file, width = 800, height = 600)
        print(barplot(ko_result))
        dev.off()
        cat("Saved indicator plot:", png_file, "\n")
      }
    }
  }
  
  # 2. 分析共进化基因 (coevolved)
  if (length(order_co_ko_list) > 0) {
    co_ko_result <- tryCatch({
      enrichKEGG(
        gene = order_co_ko_list,
        organism = "ko",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("Error in enrichKEGG for ", order, " coevolved: ", e$message)
      NULL
    })
    
    if (!is.null(co_ko_result) && !is.null(co_ko_result@result)) {
      co_ko_df <- as.data.frame(co_ko_result)
      
      # 仅在有结果时添加标识列
      if (nrow(co_ko_df) > 0) {
        co_ko_df$Category <- paste0(order, "_coevolved")
        results_list[["coevolved"]] <- co_ko_df
        
        # 保存图表
        png_file <- file.path("enrichment_plots", paste0(order, "_coevolved.png"))
        png(png_file, width = 800, height = 600)
        print(barplot(co_ko_result))
        dev.off()
        cat("Saved coevolved plot:", png_file, "\n")
      }
    }
  }
  
  # 3. 合并当前目的所有结果
  if (length(results_list) > 0) {
    order_df <- do.call(rbind, results_list)
    
    # 保存当前目的结果
    csv_file <- file.path("enrichment_results", paste0(order, "_all_results.csv"))
    write.csv(order_df, file = csv_file, row.names = FALSE)
    cat("Saved results file:", csv_file, "\n")
    
    # 添加到总结果集
    order_results <- rbind(order_results, order_df)
  }
}

# 4. 保存所有结果的合并表
if (nrow(order_results) > 0) {
  combined_file <- file.path("enrichment_results", "all_orders_combined_results.csv")
  write.csv(order_results, file = combined_file, row.names = FALSE)
  cat("Saved combined results:", combined_file, "\n")
}

########################
order_enriched <- read.csv(file = "./enrichment_results/all_orders_combined_results.csv")

order_enriched$Order <- gsub("_.*", "", order_enriched$Category)
order_enriched$Order <- factor(order_enriched$Order, levels = plant_colors$Order)
order_enriched$Group <- gsub(".*_", "", order_enriched$Category)
order_enriched$Group <- factor(order_enriched$Group, levels = c("indicators", "coevolved"))
###去除不属于微生物的富集pathway
table(order_enriched$Description)
order_enriched <- order_enriched[!order_enriched$Description %in% 
                                   c("Yersinia infection", "Th17 cell differentiation", "Th1 and Th2 cell differentiation",
                                     "T cell receptor signaling pathway","PD-L1 expression and PD-1 checkpoint pathway in cancer",
                                     "Osteoclast differentiation", "Notch signaling pathway", "NF-kappa B signaling pathway",
                                     "Adipocytokine signaling pathway"),]


order_enriched_g <- ggplot(order_enriched, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched$p.adjust)),
                       median(-log2(order_enriched$p.adjust)),
                       max(-log2(order_enriched$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order + Group, scales = "fixed", space = "fixed", switch = "x") +  #, drop = FALSE显示所有因子，包括没有富集pathway的分组
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

windows(12,6)
order_enriched_g

ggsave(order_enriched_g,
       filename = "./plots/KO_enrichment_each_plant_orders_pro_coevolved_ko.pdf",
       width = 24, height = 12, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_g,
       filename = "./plots/KO_enrichment_each_plant_orders_pro_coevolved_ko.png",
       width = 24, height = 12, units = "cm", dpi = 900, bg = "white")


####################################
familys <- family_used$Family

family_results <- data.frame()
# 循环处理每个植物目
for (family in familys) {
  cat("Processing family:", family, "\n")
  
  # 提取当前目的基因列表
  family_ko_list <- subset(pro_lamda_rf_enrich, Tree_Family == family)$KOID
  family_co_ko_list <- subset(pro_lamda_rf_enrich, Tree_Family == family & padj < 0.05)$KOID
  
  # 准备存储结果的列表
  results_list <- list()
  
  # 1. 分析指示基因 (indicators)
  if (length(family_ko_list) > 0) {
    ko_result <- tryCatch({
      enrichKEGG(
        gene = family_ko_list,
        organism = "ko",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("Error in enrichKEGG for ", family, " indicators: ", e$message)
      NULL
    })
    
    if (!is.null(ko_result) && !is.null(ko_result@result)) {
      ko_df <- as.data.frame(ko_result)
      
      # 仅在有结果时添加标识列
      if (nrow(ko_df) > 0) {
        ko_df$Category <- paste0(family, "_indicators")
        results_list[["indicators"]] <- ko_df
        
        # 保存图表
        png_file <- file.path("enrichment_plots", paste0(family, "_indicators.png"))
        png(png_file, width = 800, height = 600)
        print(barplot(ko_result))
        dev.off()
        cat("Saved indicator plot:", png_file, "\n")
      }
    }
  }
  
  # 2. 分析共进化基因 (coevolved)
  if (length(family_co_ko_list) > 0) {
    co_ko_result <- tryCatch({
      enrichKEGG(
        gene = family_co_ko_list,
        organism = "ko",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("Error in enrichKEGG for ", family, " coevolved: ", e$message)
      NULL
    })
    
    if (!is.null(co_ko_result) && !is.null(co_ko_result@result)) {
      co_ko_df <- as.data.frame(co_ko_result)
      
      # 仅在有结果时添加标识列
      if (nrow(co_ko_df) > 0) {
        co_ko_df$Category <- paste0(family, "_coevolved")
        results_list[["coevolved"]] <- co_ko_df
        
        # 保存图表
        png_file <- file.path("enrichment_plots", paste0(family, "_coevolved.png"))
        png(png_file, width = 800, height = 600)
        print(barplot(co_ko_result))
        dev.off()
        cat("Saved coevolved plot:", png_file, "\n")
      }
    }
  }
  
  # 3. 合并当前目的所有结果
  if (length(results_list) > 0) {
    family_df <- do.call(rbind, results_list)
    
    # 保存当前目的结果
    csv_file <- file.path("enrichment_results", paste0(family, "_all_results.csv"))
    write.csv(family_df, file = csv_file, row.names = FALSE)
    cat("Saved results file:", csv_file, "\n")
    
    # 添加到总结果集
    family_results <- rbind(family_results, family_df)
  }
}

# 4. 保存所有结果的合并表
if (nrow(family_results) > 0) {
  combined_file <- file.path("enrichment_results", "all_familys_combined_results.csv")
  write.csv(family_results, file = combined_file, row.names = FALSE)
  cat("Saved combined results:", combined_file, "\n")
}

########################
family_enriched <- read.csv(file = "./enrichment_results/all_familys_combined_results.csv")

family_enriched$Family <- gsub("_.*", "", family_enriched$Category)
family_enriched$Family <- factor(family_enriched$Family, levels = family_used$Family)
family_enriched$Group <- gsub(".*_", "", family_enriched$Category)
family_enriched$Group <- factor(family_enriched$Group, levels = c("indicators", "coevolved"))
###去除不属于微生物的富集pathway
table(family_enriched$Description)
family_enriched <- family_enriched[!family_enriched$Description %in% 
                                     c("Huntington disease","Yersinia infection", "Th17 cell differentiation", "Th1 and Th2 cell differentiation",
                                       "T cell receptor signaling pathway","PD-L1 expression and PD-1 checkpoint pathway in cancer",
                                       "Osteoclast differentiation", "Notch signaling pathway", "NF-kappa B signaling pathway",
                                       "Adipocytokine signaling pathway","Autophagy - animal","Coronavirus disease - COVID-19"),]


family_enriched_g <- ggplot(family_enriched, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(family_enriched$p.adjust)),
                       median(-log2(family_enriched$p.adjust)),
                       max(-log2(family_enriched$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Family + Group, scales = "fixed", space = "fixed", switch = "x") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

windows(12,6)
family_enriched_g

ggsave(family_enriched_g,
       filename = "./plots/KO_enrichment_each_plant_familys_pro_coevolved_ko.pdf",
       width = 36, height = 36, units = "cm", dpi = 900, bg = "white")

ggsave(family_enriched_g,
       filename = "./plots/KO_enrichment_each_plant_familys_pro_coevolved_ko.png",
       width = 36, height = 36, units = "cm", dpi = 900, bg = "white")


####################################################################################
classs <- class_used$Class

class_results <- data.frame()
# 循环处理每个植物目
for (class in classs) {
  cat("Processing class:", class, "\n")
  
  # 提取当前目的基因列表
  class_ko_list <- subset(pro_lamda_rf_enrich, Tree_Class == class)$KOID
  class_co_ko_list <- subset(pro_lamda_rf_enrich, Tree_Class == class & padj < 0.05)$KOID
  
  # 准备存储结果的列表
  results_list <- list()
  
  # 1. 分析指示基因 (indicators)
  if (length(class_ko_list) > 0) {
    ko_result <- tryCatch({
      enrichKEGG(
        gene = class_ko_list,
        organism = "ko",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("Error in enrichKEGG for ", class, " indicators: ", e$message)
      NULL
    })
    
    if (!is.null(ko_result) && !is.null(ko_result@result)) {
      ko_df <- as.data.frame(ko_result)
      
      # 仅在有结果时添加标识列
      if (nrow(ko_df) > 0) {
        ko_df$Category <- paste0(class, "_indicators")
        results_list[["indicators"]] <- ko_df
        
        # 保存图表
        png_file <- file.path("enrichment_plots", paste0(class, "_indicators.png"))
        png(png_file, width = 800, height = 600)
        print(barplot(ko_result))
        dev.off()
        cat("Saved indicator plot:", png_file, "\n")
      }
    }
  }
  
  # 2. 分析共进化基因 (coevolved)
  if (length(class_co_ko_list) > 0) {
    co_ko_result <- tryCatch({
      enrichKEGG(
        gene = class_co_ko_list,
        organism = "ko",
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr",
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message("Error in enrichKEGG for ", class, " coevolved: ", e$message)
      NULL
    })
    
    if (!is.null(co_ko_result) && !is.null(co_ko_result@result)) {
      co_ko_df <- as.data.frame(co_ko_result)
      
      # 仅在有结果时添加标识列
      if (nrow(co_ko_df) > 0) {
        co_ko_df$Category <- paste0(class, "_coevolved")
        results_list[["coevolved"]] <- co_ko_df
        
        # 保存图表
        png_file <- file.path("enrichment_plots", paste0(class, "_coevolved.png"))
        png(png_file, width = 800, height = 600)
        print(barplot(co_ko_result))
        dev.off()
        cat("Saved coevolved plot:", png_file, "\n")
      }
    }
  }
  
  # 3. 合并当前目的所有结果
  if (length(results_list) > 0) {
    class_df <- do.call(rbind, results_list)
    
    # 保存当前目的结果
    csv_file <- file.path("enrichment_results", paste0(class, "_all_results.csv"))
    write.csv(class_df, file = csv_file, row.names = FALSE)
    cat("Saved results file:", csv_file, "\n")
    
    # 添加到总结果集
    class_results <- rbind(class_results, class_df)
  }
}

# 4. 保存所有结果的合并表
if (nrow(class_results) > 0) {
  combined_file <- file.path("enrichment_results", "all_classs_combined_results.csv")
  write.csv(class_results, file = combined_file, row.names = FALSE)
  cat("Saved combined results:", combined_file, "\n")
}

########################
class_enriched <- read.csv(file = "./enrichment_results/all_classs_combined_results.csv")
class_enriched <- class_enriched[!is.na(class_enriched$Class), ]
class_enriched$Class <- gsub("_.*", "", class_enriched$Category)
class_enriched$Class <- factor(class_enriched$Class, levels = class_used$Class)
class_enriched$Group <- gsub(".*_", "", class_enriched$Category)
class_enriched$Group <- factor(class_enriched$Group, levels = c("indicators", "coevolved"))
###去除不属于微生物的富集pathway
table(class_enriched$Description)   

####没有富集到pathway

class_enriched <- class_enriched[!class_enriched$Description %in% 
                                   c("Yersinia infection", "Th17 cell differentiation", "Th1 and Th2 cell differentiation",
                                     "T cell receptor signaling pathway","PD-L1 expression and PD-1 checkpoint pathway in cancer",
                                     "Osteoclast differentiation", "Notch signaling pathway", "NF-kappa B signaling pathway",
                                     "Adipocytokine signaling pathway","Phototransduction - fly", "Synaptic vesicle cycle",
                                     "MAPK signaling pathway - plant","Endocrine and other factor-regulated calcium reabsorption",
                                     "Cutin, suberine and wax biosynthesis"),]


class_enriched_g <- ggplot(class_enriched, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(class_enriched$p.adjust)),
                       median(-log2(class_enriched$p.adjust)),
                       max(-log2(class_enriched$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Class + Group, scales = "fixed", space = "fixed", switch = "x") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

windows(8,6)
class_enriched_g

ggsave(class_enriched_g,
       filename = "KO_enrichment_each_plant_classs_pro_coevolved_ko.pdf",
       width = 18, height = 8, units = "cm", dpi = 900, bg = "white")

ggsave(class_enriched_g,
       filename = "KO_enrichment_each_plant_classs_pro_coevolved_ko.png",
       width = 18, height = 8, units = "cm", dpi = 900, bg = "white")
