### Settings -------------------------------------------------------------------
# Set Work Path
pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

# Set seed
set.seed(1994)

# Create directory
dir_name <- "03-absolute_abundance"
if (!file.exists(dir_name)) {dir.create(dir_name, recursive = T)}

# Import package
library(tidyverse)
library(RColorBrewer)

# Define function

## SynSpike拷贝数计算原则
### 用于扩增的DNA为微生物组DNA加上内标DNA，混合比例为30μl: 10μl

### 10μl SynSpike DNA含有的内标数如下：
# 16S 细菌 799F-1193R 内标144.48pg 拷贝数为42332891
# ITS 真菌 fITS7-ITS4 内标1pg 拷贝数为274984
# rpoB 根瘤菌 1479F-1831R 内标100pg 拷贝数为28864595
# AMF 丛枝菌根真菌 18S NS31-AML2 内标0.1pg 拷贝数为22577
# protist 原生生物 18S fw-rv 内标0.1pg 拷贝数为22712

### ITS检测的是真菌的内源转录间隔区，位于真菌18S、5.8S和28S rRNA基因之间，基于此，真菌ITS与18S的拷贝数是一致的
### 原生生物添加的SynSpike DNA过少，大多检测不到，因此使用ITS的SynSpike量按真菌的定量数据进行校正，如图方法如下：
# ITS loads = ITS reads × (ITS spike copies / spike ITS reads)
# Protist loads = ITS reads × (ITS spike copies / spike ITS reads) × (18S Protist reads / 18S Fungi reads)

# ------------------------------------------------------------------------------


### Define variable -----------------------------------------------------------
copies_16s_number <- 42332891
copies_its_number <- 274984
# ------------------------------------------------------------------------------


### Import data ----------------------------------------------------------------
spike_16s <- read.csv('../00-rawdata/data_stats/16S_stats_spike.csv', row.names = 1)
spike_its <- read.csv('../00-rawdata/data_stats/ITS_stats_spike.csv', row.names = 1)
spike_pro <- read.csv('../00-rawdata/data_stats/Protist_stats_spike.csv', row.names = 1)

metadata_rs <- read.csv('../00-rawdata/metadata/rhizosphere_metadata_merge_info.csv')
# ------------------------------------------------------------------------------


### Sort data ------------------------------------------------------------------

# ------------------------------------------------------------------------------


### Get results ----------------------------------------------------------------
## 基于样本总reads数计算copies数
### 16S
copies_16s <- merge(metadata_rs[c('FileID', 'Weight')], spike_16s, by.x = 'FileID', by.y = 'row.names')
copies_16s$Copies <- round((copies_16s$Bacteria / copies_16s$SynSpike) * copies_16s_number * (1000 / copies_16s$Weight) 
                           * (100 / 30), digits = 0)

### ITS
copies_its <- merge(metadata_rs[c('FileID', 'Weight')], spike_its, by.x = 'FileID', by.y = 'row.names')
copies_its$Copies <- round((copies_its$Fungi / copies_its$SynSpike) * copies_its_number * (1000 / copies_its$Weight) 
                           * (100 / 30), digits = 0)

### Protist
sub_copies_its <- copies_its[c('FileID', 'Fungi', 'Copies')]
names(sub_copies_its)[2:3] <- c('ITSFungi', 'ITSCopies')
copies_pro <- merge(sub_copies_its, spike_pro, by.x = 'FileID', by.y = 'row.names')

copies_pro$Copies <- round(copies_pro$ITSCopies * (copies_pro$Protist / copies_pro$Fungi), digits = 0)



### 整理用于构建跨界网络的数据集
copies_16s_sub <- copies_16s[c('FileID', 'Weight', 'SynSpike', 'Bacteria', 'Copies')]
names(copies_16s_sub)[3:5] <- c('SynSpike_16S', 'Bacteria_16S', 'Copies_16S')

copies_its_sub <- copies_its[c('FileID', 'SynSpike', 'Fungi', 'Copies')]
names(copies_its_sub)[2:4] <- c('SynSpike_ITS', 'Fungi_ITS', 'Copies_ITS')

copies_pro_sub <- copies_pro[c('FileID', 'Fungi', 'Protist', 'Copies')]
names(copies_pro_sub)[2:4] <- c('Fungi_18S', 'Protist_18S', 'Copies_Protist')


# 合并所有数据框
copies_all <- merge(copies_16s_sub, copies_its_sub)
copies_all <- merge(copies_all, copies_pro_sub)

# 筛选数据（去除Spike含量过低、序列数量过低、样本质量过低的样本）
fin_copies_all <- subset(
    copies_all, 
    Weight >= 20 & Bacteria_16S >= 10000 & Fungi_ITS >= 10000 & Fungi_18S >= 10000 & Protist_18S >= 10000 &
        SynSpike_16S / Bacteria_16S >= 0.005 & SynSpike_ITS / Fungi_ITS >= 0.005
)

name <- paste0(dir_name, "/all_sample_copies")
write.csv(fin_copies_all, paste0(name, ".csv"), quote = F, row.names = F)

# ------------------------------------------------------------------------------