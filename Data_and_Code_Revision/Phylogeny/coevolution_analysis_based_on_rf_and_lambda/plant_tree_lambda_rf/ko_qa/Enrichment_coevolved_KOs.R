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
library("clusterProfiler")


##################################################################################################
### 富集分析 ----------------------------------------------------------------

bac_core_ko <- read.csv(file="./bacteria/stats/KO_enrichment_all_bac_core_ko.csv")
bac_lamda_ko <- read.csv(file="./bacteria/stats/KO_enrichment_all_bac_lamda_ko.csv")

fun_core_ko <- read.csv(file="./fungi/stats/KO_enrichment_all_fun_core_ko.csv")
fun_lamda_ko <- read.csv(file="./fungi/stats/KO_enrichment_all_fun_lamda_ko.csv")

pro_core_ko <- read.csv(file="./protist/stats/KO_enrichment_all_pro_core_ko.csv")
pro_lamda_ko <- read.csv(file="./protist/stats/KO_enrichment_all_pro_lamda_ko.csv")

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
# 将列表中明显不属于细菌、真菌、原生生物等微生物核心功能的通路筛选出来。这些通路主要涉及高等动植物特有的生理系统、器官功能、复杂疾病及宿主特异性过程。
# 不属于微生物功能的通路列表如下：
# [7] Pathways of neurodegeneration (人类神经退行性疾病大类)
# [8] Huntington disease (人类遗传性神经退行性疾病)
# [12] Amyotrophic lateral sclerosis (人类运动神经元病)
# [14] Alzheimer disease (人类神经退行性疾病)
# [16] Parkinson disease (人类神经退行性疾病)
# [18] Spinocerebellar ataxia (人类遗传性神经疾病)
# [34] Synaptic vesicle cycle (动物神经元突触功能)
# [44] Axon regeneration (动物神经再生)
# [45] Meiosis (高等真核生物有性生殖，与微生物的简单分裂不同)
# [56] Glucagon signaling pathway (哺乳动物激素调节)
# [57] Insulin signaling pathway (哺乳动物激素调节)
# [60] Diabetic cardiomyopathy (人类糖尿病并发症)
# [70] Thermogenesis (哺乳动物/高等动物产热)
# [82] Chemical carcinogenesis (人类化学致癌)
# [84] Coronavirus disease - COVID-19 (人类传染病)
# [90] Gap junction (动物细胞连接通讯)
# [94] Tight junction (动物上皮细胞连接)
# [99] Human papillomavirus infection (人类病毒感染)
# [102] Shigellosis (侧重宿主对细菌感染的病理反应)
# [108] Dopaminergic synapse (哺乳动物神经递质系统)
# [112] Thyroid hormone signaling pathway (哺乳动物内分泌调节)
# [113] GnRH signaling pathway (哺乳动物生殖内分泌)
# [115] Longevity regulating pathway (主要基于动物模型研究)
# [116] Oxytocin signaling pathway (哺乳动物神经激素)
# [122] Proteoglycans in cancer (人类癌症相关分子)
# [123] Oocyte meiosis (动物卵母细胞减数分裂)
# [126] HIF-1 signaling pathway (动物细胞缺氧应答核心通路)
# [127] Insulin resistance (人类代谢综合征)
# [128] Sphingolipid signaling pathway (动物细胞复杂信号转导)
# [131] Adrenergic signaling in cardiomyocytes (哺乳动物心脏细胞信号)
# [132] Focal adhesion (动物细胞与基质粘附和信号枢纽)
# [133] Glutamatergic synapse (哺乳动物兴奋性突触)
# [136] Central carbon metabolism in cancer (人类肿瘤细胞代谢重编程)
# [141] Estrogen signaling pathway (哺乳动物性激素调节)
# [142] Morphine addiction (人类药物成瘾与神经适应性)
# [146] Apelin signaling pathway (哺乳动物心血管和体液平衡调节)
# [148] Salmonella infection (侧重宿主-病原体互作中的宿主反应)
# [149] IgSF CAM signaling (动物免疫球蛋白超家族细胞粘附分子)
# [150] Hepatocellular carcinoma (人类肝癌)
# [152] Alcoholism (人类行为与疾病)
# [156] Endocrine and other factor-regulated calcium reabsorption (哺乳动物肾脏精密调节)
# [157] FoxO signaling pathway (动物关键转录因子调控网络)
# [158] Adherens junction (动物细胞间粘着连接)
# [159] Salivary secretion (哺乳动物消化腺分泌)
# [160] cGMP-PKG signaling pathway (动物重要第二信使系统)
# [163] Cholinergic synapse (哺乳动物胆碱能神经系统)
# [164] Choline metabolism in cancer (人类肿瘤代谢)
# [165] Axon guidance (动物神经系统发育导向)
# [166] Long-term potentiation (动物神经元突触可塑性)
# [168] Human immunodeficiency virus 1 infection (人类病毒感染)
# [170] Drug metabolism - cytochrome P450 (人类肝脏药物代谢)
# [174] GABAergic synapse (哺乳动物抑制性神经递质系统)
# [176] Olfactory transduction (哺乳动物嗅觉感知)
# [177] Neutrophil extracellular trap formation (哺乳动物免疫细胞防御机制)
# [178] Neurotrophin signaling pathway (动物神经营养因子)
# [179] Progesterone-mediated oocyte maturation (动物生殖细胞成熟)
# [180] Vascular smooth muscle contraction (哺乳动物血管生理)
# [181] Relaxin signaling pathway (哺乳动物肽类激素)
# [182] Retrograde endocannabinoid signaling (动物神经系统反馈调节)
# [184] Growth hormone synthesis, secretion and action (哺乳动物内分泌轴)
# [186] Melanogenesis (动物黑色素合成与沉积)
# [188] Viral carcinogenesis (病毒致人类癌症机制)
# [190] Circadian entrainment (动物生物钟与环境同步)
# [192] ErbB signaling pathway (动物生长因子受体信号)
# [193] Pancreatic secretion (哺乳动物外分泌腺功能)
# [194] Endometrial cancer (人类癌症)
# [195] Gastric acid secretion (哺乳动物胃生理)
# [196] Platelet activation (哺乳动物凝血与血栓)
# [197] Aldosterone synthesis and secretion (哺乳动物肾上腺内分泌)
# [198] Lipid and atherosclerosis (人类心血管疾病)
# [202] Yersinia infection (侧重宿主对细菌感染的病理反应)
# [203] Alcoholic liver disease (人类肝脏疾病)
# [204] Neuroactive ligand-receptor interaction (动物神经递质/激素受体)
# [205] Pathogenic Escherichia coli infection (侧重宿主-病原体互作中的宿主反应)
# [206] Phagosome (动物免疫细胞的吞噬过程)
# [207] Thyroid hormone synthesis (哺乳动物内分泌腺功能)
# [208] Phototransduction (动物视觉光信号转换)
# [210] VEGF signaling pathway (动物血管生成)
# [212] Hedgehog signaling pathway (动物胚胎发育模式形成)
# [213] Fc gamma R-mediated phagocytosis (哺乳动物抗体介导的免疫吞噬)
# [216] Long-term depression (动物神经元突触可塑性)
# [217] Renin secretion (哺乳动物肾素-血管紧张素系统)
# [218] Colorectal cancer (人类癌症)
# [220] Integrin signaling (动物细胞粘附与信号转导)
# [221] Epstein-Barr virus infection (人类病毒感染)
# [224] Fc epsilon RI signaling pathway (哺乳动物过敏反应)
# [225] Insulin secretion (哺乳动物胰岛β细胞功能)
# [226] Serotonergic synapse (哺乳动物5-羟色胺能神经系统)
# [228] Wnt signaling pathway (动物发育与干细胞维持)
# [229] T cell receptor signaling pathway (哺乳动物适应性免疫核心)
# [230] Glioma (人类脑肿瘤)
# [231] Collecting duct acid secretion (哺乳动物肾脏精细调节)
# [234] AGE-RAGE signaling pathway in diabetic complications (人类糖尿病并发症分子机制)
# [235] Bacterial invasion of epithelial cells (侧重宿主细胞对细菌入侵的反应)
# [239] Proximal tubule bicarbonate reclamation (哺乳动物肾小管生理)
# [240] Renal cell carcinoma (人类肾脏癌症)
# [241] EGFR tyrosine kinase inhibitor resistance (人类癌症靶向治疗耐药)
# [242] GnRH secretion (哺乳动物下丘脑神经内分泌)
# [247] Non-alcoholic fatty liver disease (人类代谢性疾病)
# [248] Pathways in cancer (人类癌症通路的汇总)
# [249] Pancreatic cancer (人类癌症)
# [250] Type II diabetes mellitus (人类代谢性疾病)
# [252] Hippo signaling pathway (动物器官大小与增殖调控)
# [253] Human cytomegalovirus infection (人类病毒感染)
# [259] Non-small cell lung cancer (人类癌症)
# [260] Parathyroid hormone synthesis, secretion and action (哺乳动物钙磷代谢调节)
# [261] Protein digestion and absorption (哺乳动物胃肠道生理)
# [267] Epithelial cell signaling in Helicobacter pylori infection (侧重宿主细胞对感染的信号反应)
# [268] Bile secretion (哺乳动物肝脏消化液分泌)
# [269] Retinol metabolism (动物维生素A代谢)
# [270] Inflammatory mediator regulation of TRP channels (动物疼痛与炎症感知)
# [271] Chemokine signaling pathway (哺乳动物免疫细胞趋化)
# [272] Amphetamine addiction (人类精神药物成瘾)
# [273] Chronic myeloid leukemia (人类血液肿瘤)
# [275] Leukocyte transendothelial migration (哺乳动物免疫细胞浸润)
# [276] Prostate cancer (人类癌症)
# [277] Regulation of lipolysis in adipocytes (哺乳动物脂肪细胞代谢调节)
# [278] Hepatitis B (人类病毒感染)
# [280] Endocrine resistance (人类癌症激素治疗耐药)
# [281] Gastric cancer (人类癌症)
# [285] Cushing syndrome (人类内分泌疾病)
# [286] PD-L1 expression and PD-1 checkpoint pathway in cancer (人类肿瘤免疫逃逸与治疗)
# [287] Cellular senescence (主要研究动物细胞衰老)
# [288] Acute myeloid leukemia (人类血液肿瘤)
# [292] Mineral absorption (哺乳动物肠道离子吸收)
# [294] Fanconi anemia pathway (人类遗传性DNA修复疾病)
# [295] Arachidonic acid metabolism (动物类二十烷酸炎症介质代谢)
# [296] Aldosterone-regulated sodium reabsorption (哺乳动物肾脏精密调节)
# [297] Human T-cell leukemia virus 1 infection (人类病毒感染)
# [299] Secondary bile acid biosynthesis (主要在哺乳动物肝脏和肠道微生物协作完成，但归类于消化系统)
# [300] Breast cancer (人类癌症)
# [301] Cardiac muscle contraction (哺乳动物心脏生理)
# [302] Viral life cycle - HIV-1 (人类病毒)
# [303] Nicotine addiction (人类物质成瘾)
# [304] Carbohydrate digestion and absorption (哺乳动物消化生理)
# [305] Hypertrophic cardiomyopathy (人类心脏疾病)
# [306] Vasopressin-regulated water reabsorption (哺乳动物肾脏水平衡调节)
# [307] Glycosaminoglycan biosynthesis (动物结缔组织基质)
# [308] Platinum drug resistance (人类癌症化疗耐药)
# [309] Dorso-ventral axis formation (动物胚胎发育)
# [310] PPAR signaling pathway (哺乳动物代谢与细胞分化调控)
# [311] Cytoskeleton in muscle cells (动物肌肉细胞特异性结构)
# [312] Chagas disease (人类寄生虫病)
# [313] Adipocytokine signaling pathway (哺乳动物脂肪组织内分泌)
# [318] Cornified envelope formation (哺乳动物皮肤角质化)
# [320] Dilated cardiomyopathy (人类心脏疾病)
# [321] C-type lectin receptor signaling pathway (哺乳动物固有免疫识别)
# 
delete_pathways <- c(
  "Pathways of neurodegeneration", "Huntington disease", "Amyotrophic lateral sclerosis", 
  "Alzheimer disease", "Parkinson disease", "Spinocerebellar ataxia", 
  "Ribosome biogenesis in eukaryotes", "mRNA surveillance pathway", "Nucleotide excision repair", 
  "Synaptic vesicle cycle", "Mitophagy", "Axon regeneration", "Meiosis", "Nucleocytoplasmic transport", 
  "Proteasome", "Glucagon signaling pathway", "Insulin signaling pathway", "Diabetic cardiomyopathy", 
  "Thermogenesis", "mTOR signaling pathway", "Chemical carcinogenesis", "Coronavirus disease - COVID-19", 
  "Gap junction", "Tight junction", "AMPK signaling pathway", "Endocytosis", "Human papillomavirus infection", 
  "Shigellosis", "Dopaminergic synapse", "Phosphatidylinositol signaling system", 
  "Protein processing in endoplasmic reticulum", "Thyroid hormone signaling pathway", 
  "GnRH signaling pathway", "One carbon pool by folate", "Longevity regulating pathway", 
  "Oxytocin signaling pathway", "Ubiquitin mediated proteolysis", "Proteoglycans in cancer", 
  "Oocyte meiosis", "HIF-1 signaling pathway", "Insulin resistance", "Sphingolipid signaling pathway", 
  "Adrenergic signaling in cardiomyocytes", "Focal adhesion", "Glutamatergic synapse", 
  "Central carbon metabolism in cancer", "Motor proteins", "Estrogen signaling pathway", 
  "Morphine addiction", "Apelin signaling pathway", "Salmonella infection", "IgSF CAM signaling", 
  "Hepatocellular carcinoma", "Alcoholism", "Endocrine and other factor-regulated calcium reabsorption", 
  "FoxO signaling pathway", "Adherens junction", "Salivary secretion", "cGMP-PKG signaling pathway", 
  "Homologous recombination", "Cholinergic synapse", "Choline metabolism in cancer", 
  "Axon guidance", "Long-term potentiation", "Human immunodeficiency virus 1 infection", 
  "Drug metabolism", "GABAergic synapse", "Olfactory transduction", "Neutrophil extracellular trap formation", 
  "Neurotrophin signaling pathway", "Progesterone-mediated oocyte maturation", 
  "Vascular smooth muscle contraction", "Relaxin signaling pathway", "Retrograde endocannabinoid signaling", 
  "Growth hormone synthesis, secretion and action", "Melanogenesis", "Viral carcinogenesis", 
  "Circadian entrainment", "ErbB signaling pathway", "Pancreatic secretion", 
  "Endometrial cancer", "Gastric acid secretion", "Platelet activation", "Aldosterone synthesis and secretion", 
  "Lipid and atherosclerosis", "Yersinia infection", "Alcoholic liver disease", 
  "Neuroactive ligand-receptor interaction", "Pathogenic Escherichia coli infection", 
  "Phagosome", "Thyroid hormone synthesis", "Phototransduction", "cAMP signaling pathway", 
  "VEGF signaling pathway", "Hedgehog signaling pathway", "Fc gamma R-mediated phagocytosis", 
  "Long-term depression", "Renin secretion", "Colorectal cancer", "Integrin signaling", 
  "Epstein-Barr virus infection", "Fc epsilon RI signaling pathway", "Insulin secretion", 
  "Serotonergic synapse", "Wnt signaling pathway", "T cell receptor signaling pathway", 
  "Glioma", "Collecting duct acid secretion", "AGE-RAGE signaling pathway in diabetic complications", 
  "Bacterial invasion of epithelial cells", "Proximal tubule bicarbonate reclamation", 
  "Renal cell carcinoma", "EGFR tyrosine kinase inhibitor resistance", "GnRH secretion", 
  "Phospholipase D signaling pathway", "Non-alcoholic fatty liver disease", "Pathways in cancer", 
  "Pancreatic cancer", "Type II diabetes mellitus", "Hippo signaling pathway", 
  "Human cytomegalovirus infection", "Folate transport and metabolism", "Non-small cell lung cancer", 
  "Parathyroid hormone synthesis, secretion and action", "Protein digestion and absorption", 
  "ATP-dependent chromatin remodeling", "Apoptosis", "Epithelial cell signaling in Helicobacter pylori infection", 
  "Bile secretion", "Retinol metabolism", "Inflammatory mediator regulation of TRP channels", 
  "Chemokine signaling pathway", "Amphetamine addiction", "Chronic myeloid leukemia", 
  "Leukocyte transendothelial migration", "Prostate cancer", "Regulation of lipolysis in adipocytes", 
  "Hepatitis B", "Endocrine resistance", "Gastric cancer", "Cushing syndrome", 
  "PD-L1 expression and PD-1 checkpoint pathway in cancer", "Cellular senescence", 
  "Acute myeloid leukemia", "Mineral absorption", "Fanconi anemia pathway", 
  "Arachidonic acid metabolism", "Aldosterone-regulated sodium reabsorption", 
  "Human T-cell leukemia virus 1 infection", "SNARE interactions in vesicular transport", 
  "Secondary bile acid biosynthesis", "Breast cancer", "Cardiac muscle contraction", 
  "Viral life cycle - HIV-1", "Nicotine addiction", "Carbohydrate digestion and absorption", 
  "Hypertrophic cardiomyopathy", "Vasopressin-regulated water reabsorption", 
  "Glycosaminoglycan biosynthesis", "Platinum drug resistance", "Dorso-ventral axis formation", 
  "PPAR signaling pathway", "Cytoskeleton in muscle cells", "Chagas disease", 
  "Adipocytokine signaling pathway", "Cornified envelope formation", "Non-homologous end-joining", 
  "Dilated cardiomyopathy", "C-type lectin receptor signaling pathway"
)
  
inter_enriched <- inter_enriched[!inter_enriched$Description %in% delete_pathways,]
unique(inter_enriched$Description) ###154 pathways

###只展示50个pathway
inter_enriched <- subset(inter_enriched, Description %in% unique(inter_enriched$Description)[1:50])



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
  theme(plot.title = element_text(size = 8, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 8, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 8, color = 'black'),
        axis.line = element_blank(),
        axis.text.x= element_text(size = 6, color = 'black'),
        axis.text.y = element_text(size = 8, color = "black"),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')

windows(8,12)
inter_enriched_g

ggsave(inter_enriched_g,
       filename = "KO_enrichment_Core_and_Phylogenetically_conserved_ko-50-pathway.pdf",
       width = 240, height = 300, units = "mm", dpi = 900, bg = "white")
ggsave(inter_enriched_g,
       filename = "KO_enrichment_Core_and_Phylogenetically_conserved_ko-50-pathway.png",
       width = 240, height = 300, units = "mm", dpi = 900, bg = "white")


#############################################################################################
#############################################################################################
# 1. 分别筛选每个组别中 qvalue 最小的前20个 pathway
# Bacterial PCG 组
bacterial_pathways <- inter_enriched %>%
  filter(Group == "Bacterial PCG") %>%
  arrange(qvalue) %>%
  slice_head(n = 20) %>%
  pull(Description)

# Fungal PCG 组  
fungal_pathways <- inter_enriched %>%
  filter(Group == "Fungal PCG") %>%
  arrange(qvalue) %>%
  slice_head(n = 20) %>%
  pull(Description)

# Protist PCG 组
protist_pathways <- inter_enriched %>%
  filter(Group == "Protist PCG") %>%
  arrange(qvalue) %>%
  slice_head(n = 20) %>%
  pull(Description)

# 2. 合并所有选中的 pathway（去重，因为不同组可能有相同的 pathway）
selected_pathways <- unique(c(bacterial_pathways, fungal_pathways, protist_pathways))

# 3. 筛选 inter_enriched 数据框，只保留这些 pathway
inter_enriched_selected <- inter_enriched %>%
  filter(Description %in% selected_pathways)

# cell cycle 有三个，选ko04112代表微生物的
inter_enriched_selected <- inter_enriched_selected %>%
  filter(!(ID %in% c("ko04110","ko04111")))

inter_enriched_g2 <- ggplot(inter_enriched_selected, aes(x=-log2(qvalue), y=Description)) +
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
  theme(plot.title = element_text(size = 8, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 8, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 8, color = 'black'),
        axis.line = element_blank(),
        axis.text.x= element_text(size = 6, color = 'black'),
        axis.text.y = element_text(size = 8, color = "black"),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')

windows(8,12)
inter_enriched_g2

ggsave(inter_enriched_g2,
       filename = "KO_enrichment_Core_and_Phylogenetically_conserved_ko-each-kingdom-top20.pdf",
       width = 240, height = 240, units = "mm", dpi = 900, bg = "white")
ggsave(inter_enriched_g2,
       filename = "KO_enrichment_Core_and_Phylogenetically_conserved_ko-each-kingdom-top20.png",
       width = 240, height = 240, units = "mm", dpi = 900, bg = "white")




#####################################################################################################################
#####################################################################################################################
############################三界分开，植物order水平
theme_used <-
  theme(plot.title = element_text(size = 8, color = 'black', hjust = 0.5),
        plot.subtitle = element_text(size = 8, color = 'black', hjust = 0.5),
        axis.title = element_text(size = 8, color = 'black'),
        axis.line = element_blank(),
        axis.text = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8,  color = 'black'),
        legend.key.size = unit(0.25, 'cm'),
        # legend.margin = margin(0, 0, 0, -0.2, 'cm'),
        legend.position = 'right')



########################bacteria
order_enriched_bacterial <- read.csv(file = "./bacteria/enrichment_results/all_orders_combined_results.csv")

order_enriched_bacterial$Order <- gsub("_.*", "", order_enriched_bacterial$Category)
order_enriched_bacterial$Order <- factor(order_enriched_bacterial$Order, levels = plant_colors$Order)
order_enriched_bacterial$Group <- gsub(".*_", "", order_enriched_bacterial$Category)
order_enriched_bacterial$Group <- factor(order_enriched_bacterial$Group, levels = c("indicators", "coevolved"))
###去除不属于微生物的富集pathway
order_enriched_bacterial <- order_enriched_bacterial[!order_enriched_bacterial$Description %in% delete_pathways,]


order_enriched_bacterial_g <- ggplot(order_enriched_bacterial, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched_bacterial$p.adjust)),
                       median(-log2(order_enriched_bacterial$p.adjust)),
                       max(-log2(order_enriched_bacterial$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order + Group, scales = "fixed", space = "fixed", switch = "x") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

windows(12,6)
order_enriched_bacterial_g

ggsave(order_enriched_bacterial_g,
       filename = "KO_enrichment_each_plant_orders_bacterial_indicators_and_coevolved_ko.pdf",
       width = 24, height = 16, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_bacterial_g,
       filename = "KO_enrichment_each_plant_orders_bacterial_indicators_and_coevolved_ko.png",
       width = 24, height = 16, units = "cm", dpi = 900, bg = "white")



#################################################
###只绘制co-evolved， 没有的也画一个框
order_enriched_bacterial2 <- subset(order_enriched_bacterial, Group == "coevolved")
order_enriched_bacterial2$Order <- factor(order_enriched_bacterial2$Order, levels = plant_colors$Order)
length <- length(unique(order_enriched_bacterial2$Description))
height <- length/46 * 16

order_enriched_bacterial_g2 <- ggplot(order_enriched_bacterial2, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched_bacterial$p.adjust)),
                       median(-log2(order_enriched_bacterial$p.adjust)),
                       max(-log2(order_enriched_bacterial$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order, scales = "fixed", space = "fixed", switch = "x", drop = FALSE) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

ggsave(order_enriched_bacterial_g2,
       filename = "KO_enrichment_each_plant_orders_bacterial_coevolved_ko.pdf",
       width = 24, height = height, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_bacterial_g2,
       filename = "KO_enrichment_each_plant_orders_bacterial_coevolved_ko.png",
       width = 24, height = height, units = "cm", dpi = 900, bg = "white")



########################
### fungi
order_enriched_fungal <- read.csv(file = "./fungi/enrichment_results/all_orders_combined_results.csv")
plant_colors <- read.csv("../../../metadata/Tree_top_order_color.csv")


order_enriched_fungal$Order <- gsub("_.*", "", order_enriched_fungal$Category)
order_enriched_fungal$Order <- factor(order_enriched_fungal$Order, levels = plant_colors$Order)
order_enriched_fungal$Group <- gsub(".*_", "", order_enriched_fungal$Category)
order_enriched_fungal$Group <- factor(order_enriched_fungal$Group, levels = c("indicators", "coevolved"))

###去除不属于微生物的富集pathway
table(order_enriched_fungal$Description)
order_enriched_fungal <- order_enriched_fungal[!order_enriched_fungal$Description %in% delete_pathways,]


order_enriched_fungal_g <- ggplot(order_enriched_fungal, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched_fungal$p.adjust)),
                       median(-log2(order_enriched_fungal$p.adjust)),
                       max(-log2(order_enriched_fungal$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order + Group, scales = "fixed", space = "fixed", switch = "x") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

windows(12,6)
order_enriched_fungal_g

ggsave(order_enriched_fungal_g,
       filename = "KO_enrichment_each_plant_orders_fungal_indicators_and_coevolved_ko.pdf",
       width = 24, height = 8, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_fungal_g,
       filename = "KO_enrichment_each_plant_orders_fungal_indicators_and_coevolved_ko.png",
       width = 24, height = 8, units = "cm", dpi = 900, bg = "white")

#################################################
###只绘制co-evolved， 没有的也画一个框
order_enriched_fungal2 <- subset(order_enriched_fungal, Group == "coevolved")
order_enriched_fungal2$Order <- factor(order_enriched_fungal2$Order, levels = plant_colors$Order)
length <- length(unique(order_enriched_fungal2$Description))
height <- length/46 * 16


order_enriched_fungal_g2 <- ggplot(order_enriched_fungal2, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched_fungal$p.adjust)),
                       median(-log2(order_enriched_fungal$p.adjust)),
                       max(-log2(order_enriched_fungal$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order, scales = "fixed", space = "fixed", switch = "x", drop = FALSE) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

ggsave(order_enriched_fungal_g2,
       filename = "KO_enrichment_each_plant_orders_fungal_coevolved_ko.pdf",
       width = 24, height = height + 2, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_fungal_g2,
       filename = "KO_enrichment_each_plant_orders_fungal_coevolved_ko.png",
       width = 24, height = height + 2, units = "cm", dpi = 900, bg = "white")



########################protists
order_enriched_protistan <- read.csv(file = "./protist/enrichment_results/all_orders_combined_results.csv")

order_enriched_protistan$Order <- gsub("_.*", "", order_enriched_protistan$Category)
order_enriched_protistan$Order <- factor(order_enriched_protistan$Order, levels = plant_colors$Order)
order_enriched_protistan$Group <- gsub(".*_", "", order_enriched_protistan$Category)
order_enriched_protistan$Group <- factor(order_enriched_protistan$Group, levels = c("indicators", "coevolved"))
###去除不属于微生物的富集pathway
table(order_enriched_protistan$Description)
order_enriched_protistan <- order_enriched_protistan[!order_enriched_protistan$Description %in% delete_pathways,]


order_enriched_protistan_g <- ggplot(order_enriched_protistan, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched_protistan$p.adjust)),
                       median(-log2(order_enriched_protistan$p.adjust)),
                       max(-log2(order_enriched_protistan$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order + Group, scales = "fixed", space = "fixed", switch = "x") +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

windows(12,6)
order_enriched_protistan_g

ggsave(order_enriched_protistan_g,
       filename = "KO_enrichment_each_plant_orders_protistan_indicators_and_coevolved_ko.pdf",
       width = 24, height = 12, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_protistan_g,
       filename = "KO_enrichment_each_plant_orders_protistan_indicators_and_coevolved_ko.png",
       width = 24, height = 12, units = "cm", dpi = 900, bg = "white")


#################################################
###只绘制co-evolved， 没有的也画一个框
order_enriched_protistan2 <- subset(order_enriched_protistan, Group == "coevolved")
order_enriched_protistan2$Order <- factor(order_enriched_protistan2$Order, levels = plant_colors$Order)
length <- length(unique(order_enriched_protistan2$Description))
height <- length/46 * 16


order_enriched_protistan_g2 <- ggplot(order_enriched_protistan2, aes(x=-log2(p.adjust), y=Description)) +
  geom_point(aes(size=Count, color=-log2(p.adjust))) +
  scale_colour_gradientn(
    colours = colorRampPalette(c("darkgreen", "red"))(1000),
    values = rescale(c(min(-log2(order_enriched_protistan$p.adjust)),
                       median(-log2(order_enriched_protistan$p.adjust)),
                       max(-log2(order_enriched_protistan$p.adjust))))) +
  scale_size_continuous(name="gene counts") +
  labs(title="", x="-log2(p.adjust)", y="") +
  facet_grid( ~ Order, scales = "fixed", space = "fixed", switch = "x", drop = FALSE) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 5, color = "black"))+
  theme_used 

ggsave(order_enriched_protistan_g2,
       filename = "KO_enrichment_each_plant_orders_protistan_coevolved_ko.pdf",
       width = 24, height = height+2, units = "cm", dpi = 900, bg = "white")

ggsave(order_enriched_protistan_g2,
       filename = "KO_enrichment_each_plant_orders_protistan_coevolved_ko.png",
       width = 24, height = height+2, units = "cm", dpi = 900, bg = "white")





