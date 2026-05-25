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
# BiocManager::install('clusterProfiler')
# BiocManager::install("DOSE")
library("clusterProfiler")

### 设置工作目录
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

### Figure S26
##################################################################################################
#####################富集分析
### 富集分析 ----------------------------------------------------------------
bac_core_ko <- read.csv(file="./data/KO_enrichment_all_bac_core_ko.csv")
bac_lamda_ko <- read.csv(file="./data/KO_enrichment_all_bac_lamda_ko.csv")

fun_core_ko <- read.csv(file="./data/KO_enrichment_all_fun_core_ko.csv")
fun_lamda_ko <- read.csv(file="./data/KO_enrichment_all_fun_lamda_ko.csv")

pro_core_ko <- read.csv(file="./data/KO_enrichment_all_pro_core_ko.csv")
pro_lamda_ko <- read.csv(file="./data/KO_enrichment_all_pro_lamda_ko.csv")

# 绘制气泡图

####整合图
inter_enriched <- rbind.data.frame(data.frame(bac_core_ko),
                                   data.frame(bac_lamda_ko),
                                   data.frame(fun_core_ko), 
                                   data.frame(fun_lamda_ko),
                                   data.frame(pro_core_ko),
                                   data.frame(pro_lamda_ko))

# Phylogenetically conserved genes = PCG  
inter_enriched$Group <- c(rep("Bacterial core", dim(data.frame(bac_core_ko))[1]), 
                          rep("Bacterial PCG", dim(data.frame(bac_lamda_ko))[1]),
                          rep("Fungal core", dim(data.frame(fun_core_ko))[1]),
                          rep("Fungal PCG", dim(data.frame(fun_lamda_ko))[1]),
                          rep("Protist core", dim(data.frame(pro_core_ko))[1]),
                          rep("Protist PCG", dim(data.frame(pro_lamda_ko))[1]))

inter_enriched$Group <- factor(inter_enriched$Group, levels = c("Bacterial core", "Bacterial PCG", 
                                                                "Fungal core", "Fungal PCG",
                                                                "Protist core", "Protist PCG"))
inter_enriched <- inter_enriched[order(inter_enriched$p.adjust),]
inter_enriched$Description <- factor(inter_enriched$Description, levels = unique(inter_enriched$Description))

###去除非微生物来源pathway
inter_enriched <- inter_enriched[!inter_enriched$Description %in% c("Huntington disease", "Diabetic cardiomyopathy","Amyotrophic lateral sclerosis (ALS)",
                                                                    "Synaptic vesicle cycle","Parkinson disease","Prion disease",
                                                                    "Alzheimer disease","Spinocerebellar ataxia","Ribosome biogenesis in eukaryotes",
                                                                    "Autophagy - animal", "Insulin signaling pathway", "Oocyte meiosis",
                                                                    "Vascular smooth muscle contraction", "Cardiac muscle contraction", "Thyroid hormone synthesis/signaling",
                                                                    "Adrenergic signaling in cardiomyocytes", "Proximal tubule bicarbonate reclamation", "Collecting duct acid secretion",
                                                                    "Alcoholic liver disease", "Lipid and atherosclerosis", "Phototransduction - fly",
                                                                    "Spliceosome", "Phagosome","Human infectious disease pathways"),]



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

###只展示50个pathway
inter_enriched <- subset(inter_enriched, Description %in% unique(inter_enriched$Description)[1:50])



p <- ggplot(inter_enriched, aes(x=-log2(qvalue), y=Description)) +
  geom_point(aes(size=log(Count,base=2), color=-log2(qvalue))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("white", "red"))(1000),
    values = rescale(c(min(-log2(inter_enriched$qvalue)),
                       median(-log2(inter_enriched$qvalue)),
                       max(-log2(inter_enriched$qvalue))))) +
  scale_size_continuous(name="log2(gene counts)") +
  labs(title="", x="-log2(qvalue)", y="") +
  facet_grid(~ Group) +
  theme_bw() + theme(
      text = element_text(color = "black", size = 6),
      plot.title = element_text(size = 7, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(color = "black", size = 7, margin = margin(0.1, 0.1, 0.1, 0.1, "cm")),
      panel.spacing = unit(0.1, "cm"),
      legend.box.spacing = unit(0.1,"cm"),
      legend.key.size = unit(0.25, "cm"),
      legend.position = "right"
  )


width <- 17.5
height <- 20

name <- paste0("FigureS26")
ggsave(paste0(name, ".png"), p, width = width, height = height, dpi = 600, units = "cm")
ggsave(paste0(name, ".pdf"), p, width = width, height = height, units = "cm")


