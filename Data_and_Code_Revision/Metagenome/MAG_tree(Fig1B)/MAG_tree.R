# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "result"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(igraph)
library(ggraph)
library(ggtree)
library(ggtreeExtra)
library(grid)
library(tidytree)
library(treeio)
library(cowplot)

tree <- treeio::read.newick("./data/bacteria.nwk",node.label = "support")
group_file<-read.csv("./data/group.csv",row.names = 1)
label<-tree@phylo$tip.label
group_file <- group_file[rownames(group_file) %in% label,]
bin_tax <- read.csv("./data/MAG_taxonomy.csv",row.names = 1)
group_file$Class <- NULL
group_file$Class <- bin_tax[rownames(group_file), "Class"]
group_file <- group_file %>%
  rownames_to_column("bin") %>%
  add_count(Class, name = "freq") %>% 
  mutate(is_top10 = Class %in% head(levels(reorder(Class, -freq)), 10) & Class != "Unassign") %>%
  mutate(Class_new = case_when(
    is_top10 ~ as.character(Class),
    Class == "Unassign" ~ "Unassign",
    TRUE ~ "Others"
  )) %>%
  mutate(Class_new = factor(
    Class_new,
    levels = c(
      unique(Class[is_top10]),
      "Others",
      "Unassign"
    )
  )) %>%
  
  # 移除临时列
  select(-freq, -is_top10) %>%
  column_to_rownames("bin")
groupInfo <- split(row.names(group_file), group_file$Class_new)
tree <- groupOTU(tree, groupInfo)
group <- group_file

p1 <- ggtree(tree,layout = "fan",right = TRUE,open.angle = 15)+
  geom_aline(aes(color=group),linetype = "solid",size=0.2,key_glyph = 'rect',show.legend = TRUE)+
  #  scale_color_manual(values=c("#8DD3C7","#F8F8B6","#CAAEC5","#D68E8F","#D8C965","#DED3B3","#E3D5DC","#C191C2","#CAE0C4","#FFED6F","grey80"), breaks = levels(group_file$Class_new)) +
  scale_color_manual(values=c("#E41A1C","#596A98","#449B75","#6B886D","#AC5782","#FF7F00","#FFE528","#C9992C","#C66764","#E485B7","#666666","#999999"), breaks = levels(group_file$Class_new))+
  guides(
    colour = guide_legend(
      title = "Class",
      title.position = "top",
      title.theme = element_text(
        face = "bold",
        size = 12,
        margin = margin(b = 20),
        family = "sans",
        colour = "black"
      )
    )
  )+
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key.spacing.y = unit(1.5,"mm"),
    legend.key.size = unit(1,"line")
  )

SGB_plot <- group %>%
  rownames_to_column("label") %>%
  select(label,SGB) %>%
  mutate(count=1)
p2<-rotate_tree(p1 + guides(colour ='none'),angle = 170)+
  geom_fruit(
    data= SGB_plot,
    geom = geom_bar,
    stat ='identity',
    width =0.7,
    aes(y=label,x=count,fill =SGB),
    pwidth = 0.05
  ) +
  guides(fill=guide_legend('SGB'))+
  theme(
    legend.position ='left',
    legend.key.size =unit(1,'line'),
    legend.title =element_text(
      face ='bold',
      size = 12,
      margin =margin(b=20),
      family ='sans',
      colour ='black'
    )
  )+
  scale_fill_manual(values =c("#254689","#B5382C","#D18C39"))
size_plot <- group %>%
  rownames_to_column("label") %>%
  select(label,size)
p3 <- p2 + geom_fruit(
  data = size_plot,
  geom = geom_bar,
  stat = 'identity',
  width = 0.7,
  aes(y = label, x=size),
  pwidth = 0.25
)

tree.legend <- visbuilder::get_legend(p1)
SGB.legend <- visbuilder::get_legend(p2)


title.legend <- grobTree(
  textGrob(
    label = 'Tree scales 0.1:',
    x = 0.1, y = 0.5,
    hjust = 0, vjust = 0.5,
    gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  segmentsGrob(
    x0 = 0.7, x1 = 0.8,
    y0 = 0.5, y1 = 0.5,
    gp = gpar(lwd = 2)
  ),
  segmentsGrob(
    x0 = 0.7, x1 = 0.7,
    y0 = 0.48, y1 = 0.52,
    gp = gpar(lwd = 2)
  ),
  segmentsGrob(
    x0 = 0.8, x1 = 0.8,
    y0 = 0.48, y1 = 0.52,
    gp = gpar(lwd = 2)
  )
)


size.legend <- grobTree(
  textGrob(
    label ='Genome size(MB)',
    x = 0.1, y = 0.5,
    hjust = 0, vjust = 0.5,
    gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  rectGrob(
    x = 0.19, y = 0.3, width = 0.15, height = 0.075,
    gp = gpar(fill = '#000000', col = NA)
  ),
  textGrob(
    label = '8799 Bacteria MAGs',
    x = 0.31, y = 0.315,
    hjust = 0, vjust = 1,
    gp = gpar(fontsize = 9)
  )
)

p4 <- ggdraw() +
  draw_plot(p3 + guides(fill = 'none'), x = 0.15, y = 0, width = 0.90, height = 1.2) +
  draw_plot(tree.legend, x = 0.018, y = 0.85, width = 0.2, height = 0.2, vjust = 1) +
  draw_line(x = c(0.055, 0.160), y = c(0.87, 0.87), color = 'black', size = .6) +
  draw_plot(SGB.legend,x = 0.032, y = 0.615, width = 0.1, height = 0.2, vjust = 1) + 
  draw_line(x = c(0.055, 0.160), y = c(0.54, 0.54), color = 'black', size = .6) +
  draw_plot(size.legend,x = 0.045, y = 0.52, width = 0.1, height = 0.2, vjust = 1) +
  draw_line(x = c(0.055, 0.160), y = c(0.4, 0.4), color = 'black', size = .6) +
  draw_label("MAGs of Bacteria",x = 0.5, y = 0.98,hjust = 0.5, vjust = 1,size = 18,fontface = "bold")    

width = 12; height = 12
ggsave(filename = "./MAG_tree.pdf", p4, width = width, height = height, dpi = 600, device = 'pdf', bg = 'white')
ggsave(filename = "./MAG_tree.jpg", p4, width = width, height = height, units = "cm",dpi=600)