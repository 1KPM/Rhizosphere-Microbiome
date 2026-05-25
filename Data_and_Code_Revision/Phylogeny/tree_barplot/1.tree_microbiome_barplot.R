### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()[["path"]])
setwd(pwd)

# Import package
library(ggtree)
library(ggplot2)
library(reshape2)
library(agricolae)
library(picante)
library(ape)
library(vegan)
library(treeio)
library(RColorBrewer)
library(cowplot)

#################################################################################################################
### 微生物分类颜色
color_manual <- c(colorRampPalette(brewer.pal(9, "Set1"))(11), "#666666")

### 数据导入
class_16S_ra <- read.csv(file = "../metadata/16S_top10Class_barplot_data.csv")
class_ITS_ra <- read.csv(file = "../metadata/ITS_top10Class_barplot_data.csv")
class_18S_ra <- read.csv(file = "../metadata/Protist_top10Class_barplot_data.csv")

# 将相对丰度转换为绝对丰度
fin_copies_all <- read.csv('../metadata/All_sample_copies.csv', row.names = 1)
tree_file <- read.tree("../metadata/tree_metadata_merge_info_final_align_tree.nwk")
metadata_rs <- read.csv("../metadata/rhizosphere_metadata_merge_info.csv")
metadata_tree <- read.csv("../metadata/tree_metadata_merge_info.csv")

# 读取top10植物目信息
top_tree <- read.csv("../metadata/Tree_top_order_color.csv", fileEncoding = "GBK")
top10orders <- top_tree$Order[1:10]

### 根据tree ID取均值
fin_copies_rs <- merge.data.frame(metadata_rs[,1:2], fin_copies_all, by.x = "FileID", by.y = "row.names", all.y = T)
fin_copies_tree <- aggregate(fin_copies_rs[,-c(1:2)], by = list(Tree = fin_copies_rs$TreeID), FUN = mean)

# 转换为绝对丰度
class_16S_qa <- merge.data.frame(class_16S_ra, fin_copies_tree[,c("Tree", "Copies_16S")], by = "Tree", all.x = T)
class_16S_qa$QA <- class_16S_qa$Value * class_16S_qa$Copies_16S
names(class_16S_qa)[3] <- "RA" 

class_ITS_qa <- merge.data.frame(class_ITS_ra, fin_copies_tree[,c("Tree", "Copies_ITS")], by = "Tree", all.x = T)
class_ITS_qa$QA <- class_ITS_qa$Value * class_ITS_qa$Copies_ITS
names(class_ITS_qa)[3] <- "RA" 

class_18S_qa <- merge.data.frame(class_18S_ra, fin_copies_tree[,c("Tree", "Copies_Protist")], by = "Tree", all.x = T)
class_18S_qa$QA <- class_18S_qa$Value * class_18S_qa$Copies_Protist
names(class_18S_qa)[3] <- "RA" 

# 过滤NA值
class_16S <- class_16S_qa[!is.na(class_16S_qa$QA), ]
class_ITS <- class_ITS_qa[!is.na(class_ITS_qa$QA), ]
class_18S <- class_18S_qa[!is.na(class_18S_qa$QA), ]

########################################################################################
### 循环可视化所有 top10orders 

# 确保输出目录存在
output_dir <- "./tree_barplots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 循环处理每个目
for (order_name in top10orders) {
  cat("正在处理:", order_name, "\n")
  
  ###########################################
  # 1. 筛选该目的植物ID
  ###########################################
  order_id <- subset(metadata_tree, Order == order_name)$TreeID
  
  # 取三个数据集都有的植物ID
  order_16S_id <- intersect(order_id, class_16S$Tree)
  order_ITS_id <- intersect(order_id, class_ITS$Tree)
  order_18S_id <- intersect(order_id, class_18S$Tree)
  order_id_used <- intersect(order_16S_id, intersect(order_18S_id, order_ITS_id))
  
  if (length(order_id_used) == 0) {
    cat("  ", order_name, "没有可用数据，跳过\n")
    next
  }
  
  ###########################################
  # 2. 过滤系统发育树
  ###########################################
  # 需要移除的tip
  order_to_drop <- tree_file$tip.label[-match(order_id_used, tree_file$tip.label)]
  order_tree <- treeio::drop.tip(tree_file, order_to_drop)
  
  ###########################################
  # 3. 准备元数据
  ###########################################
  metadata_order <- metadata_tree[metadata_tree$TreeID %in% order_tree$tip.label, ] 
  row.names(metadata_order) <- metadata_order$TreeID
  
  ###########################################
  # 4. 创建分组信息（按Family或Genus）
  ###########################################
  # 对于Fabales和Arecales，按Genus分组，其他目按Family分组
  if (order_name %in% c("Fabales", "Arecales")) {
    # 按Genus分组
    # 计算每个Genus的数量
    genus_counts <- table(metadata_order$Genus)
    
    # 如果Genus数量大于10，取前10，其他归类为Others
    if (length(genus_counts) > 10) {
      # 按数量降序排序
      sorted_genus <- sort(genus_counts, decreasing = TRUE)
      top_genus <- names(sorted_genus[1:10])
      
      # 创建分组变量
      metadata_order$Group <- ifelse(metadata_order$Genus %in% top_genus, 
                                     metadata_order$Genus, "Others")
    } else {
      metadata_order$Group <- metadata_order$Genus
    }
    
    # 创建分组信息
    groupInfo <- split(metadata_order$TreeID, metadata_order$Group)
    
  } else {
    # 按Family分组
    metadata_order$Group <- metadata_order$Family
    groupInfo <- split(metadata_order$TreeID, metadata_order$Family)
  }
  
  order_tree <- groupOTU(order_tree, groupInfo)
  
  # 创建标签
  metadata_order$Label <- ifelse(metadata_order$Species != "", metadata_order$Species, 
                                 ifelse(metadata_order$Genus != "", metadata_order$Genus,
                                        ifelse(metadata_order$Family != "", metadata_order$Family, 
                                               metadata_order$Order)))
  metadata_order$Label <- gsub("×", "", metadata_order$Label)
  metadata_order$Label <- paste(metadata_order$TreeID, metadata_order$Label, sep = " | ")
  
  ###########################################
  # 5. 准备颜色映射
  ###########################################
  # 获取唯一分组
  unique_groups <- unique(metadata_order$Group)
  n_colors <- length(unique_groups)
  
  # 设置颜色
  if (order_name %in% c("Fabales", "Arecales")) {
    # 对于Fabales和Arecales，使用Set3调色板
    if (n_colors <= 12) {
      group_colors <- brewer.pal(ifelse(n_colors >= 3, n_colors, 3), "Set3")[1:n_colors]
    } else {
      group_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
    }
    
    # 将Others设置为灰色
    if ("Others" %in% unique_groups) {
      # 找到Others的位置
      others_idx <- which(unique_groups == "Others")
      # 设置Others为灰色
      group_colors[others_idx] <- "gray70"
    }
    
    names(group_colors) <- unique_groups
  } else {
    # 其他目按原逻辑
    if (n_colors <= 8) {
      group_colors <- brewer.pal(n_colors, "Set2")
    } else if (n_colors <= 12) {
      group_colors <- brewer.pal(n_colors, "Set3")
    } else {
      group_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_colors)
    }
    names(group_colors) <- unique_groups
  }
  
  ###########################################
  # 6. 准备微生物数据
  order_16S_class <- class_16S[class_16S$Tree %in% order_id_used, ]
  order_16S_class$Taxonomy <- factor(order_16S_class$Taxonomy, levels = unique(class_16S_ra$Taxonomy))
  
  order_ITS_class <- class_ITS[class_ITS$Tree %in% order_id_used, ]
  order_ITS_class$Taxonomy <- factor(order_ITS_class$Taxonomy, levels = unique(class_ITS_ra$Taxonomy))
  
  order_18S_class <- class_18S[class_18S$Tree %in% order_id_used, ]
  order_18S_class$Taxonomy <- factor(order_18S_class$Taxonomy, levels = unique(class_18S_ra$Taxonomy))
  
  ###########################################
  # 7. 绘制系统发育树
  ###########################################
  # 设置图例名称
  legend_name <- ifelse(order_name %in% c("Fabales", "Arecales"), "Genus", "Family")
  
  order_tree_g <- ggtree(order_tree, aes(color=group), 
                         layout = "rectangular", branch.length = "none", 
                         size = 0.3, ladderize = F) + 
    theme(legend.position = c(-0.5, 0.5),
          legend.justification = c(-0.5, 0.5),
          legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(size = 6, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          legend.background = element_blank(),
          legend.margin = margin(2, 2, 2, 2)) +
    xlab("") +
    xlim(-10, 35) +
    scale_color_manual(values = group_colors, name = paste0(order_name, "\n(", legend_name, ")")) +
    theme(axis.title.x = element_text(size = 8, color = "black")) +
    theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
  
  # 添加标签
  order_tree_label_g <- order_tree_g %<+% metadata_order +
    geom_tiplab(aes(label = Label), size = 2, offset = 0, 
                hjust = 0, fontface = "italic") +
    theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
  
  ###########################################
  # 8. 绘制条形图
  ###########################################
  # 创建函数来绘制条形图
  create_bar_plot <- function(data, y_label, fill_manual) {
    ggplot(data, aes(x = Tree, y = QA, fill = Taxonomy)) +
      geom_bar(stat = "identity", position = position_stack(reverse = F), width = 1) +
      labs(x = "", y = y_label) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = fill_manual, name = y_label) +
      coord_flip() +
      theme_minimal() +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0.2, 0.1, 0, -0.3), "cm"), #上右下左
            panel.spacing = unit(0, "cm")) +
      theme(axis.title.x = element_text(size = 6, color = "black")) +  
      scale_x_discrete(expand = expansion(0))  +
      scale_y_discrete(expand = expansion(0))
  }
  
  # 16S条形图
  order_16S_class_g <- create_bar_plot(order_16S_class, "Bacteria", color_manual)
  
  # ITS条形图
  order_ITS_class_g <- create_bar_plot(order_ITS_class, "Fungi", color_manual)
  
  # 18S条形图
  order_18S_class_g <- create_bar_plot(order_18S_class, "Protists", color_manual)
  
  ###########################################
  # 9. 组合图形
  ###########################################
  order_cowplot <- cowplot::plot_grid(
    order_tree_label_g, 
    order_16S_class_g,
    order_ITS_class_g, 
    order_18S_class_g, 
    ncol = 4, 
    rel_widths = c(2, 0.8, 0.8, 0.8),
    align = "h",
    axis = "tb",
    greedy = FALSE
  )
  
  ###########################################
  # 10. 保存图形
  ###########################################
  # 动态计算高度（基于Sapindales的参考高度）
  base_height <- 9
  sapindales_count <- 33  # Sapindales的参考数量
  current_count <- length(order_id_used)
  height <- base_height * current_count / sapindales_count
  
  # 生成文件名
  file_base <- file.path(output_dir, paste0(order_name, "_phylogenetic_tree_and_bacteria_fungi_protist_bar"))
  
  # 保存PNG
  ggsave(paste0(file_base, ".png"), order_cowplot, 
         width = 9, height = height, dpi = 600, units = "cm", bg = "white")
  
  # 保存PDF
  ggsave(paste0(file_base, ".pdf"), order_cowplot, 
         width = 9, height = height, units = "cm", bg = "white")
  
  cat("  ", order_name, "处理完成，保存了", current_count, "个植物\n")
}

cat("\n所有目处理完成！\n")