pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(vegan)

##Bacteria
#Alpha Diversity
abun_table <- read.csv(file="../kegg_abundance/absolute/Bacteria_RS_quantitative_abundance.csv",check.names = F,row.names = 1,header=T)
shannon<- diversity(abun_table,index="shannon", MARGIN = 2)
simpson<- diversity(abun_table,index="simpson", MARGIN = 2)
richness<-specnumber(abun_table, MARGIN = 2)
evenness <- shannon / log(richness)
index <- as.data.frame(cbind(shannon, simpson, richness, evenness))
write.csv(index, 'Bacteria.quantitative.diversity.index.csv', quote = F)

#Beta Diversity
df <- t(abun_table)
distance <- as.matrix(vegdist(df,method = "bray"))
write.csv(distance,file="Bacteria.quantitative.bray.csv")
distance <- as.matrix(vegdist(df,method = "jaccard",binary = T))
write.csv(distance,file="Bacteria.quantitative.jaccard.csv")

##Fungi
#Alpha Diversity
abun_table <- read.csv(file="../kegg_abundance/absolute/Fungi_RS_quantitative_abundance.csv",check.names = F,row.names = 1,header=T)
shannon<- diversity(abun_table,index="shannon", MARGIN = 2)
simpson<- diversity(abun_table,index="simpson", MARGIN = 2)
richness<-specnumber(abun_table, MARGIN = 2)
evenness <- shannon / log(richness)
index <- as.data.frame(cbind(shannon, simpson, richness, evenness))
write.csv(index,'Fungi.quantitative.diversity.index.csv', quote = F)

#Beta Diversity
df <- t(abun_table)
distance <- as.matrix(vegdist(df,method = "bray"))
write.csv(distance,file="Fungi.quantitative.bray.csv")
distance <- as.matrix(vegdist(df,method = "jaccard",binary = T))
write.csv(distance,file="Fungi.quantitative.jaccard.csv")

##Protist
#Alpha Diversity
abun_table <- read.csv(file="../kegg_abundance/absolute/Protist_RS_quantitative_abundance.csv",check.names = F,row.names = 1,header=T)
shannon<- diversity(abun_table,index="shannon", MARGIN = 2)
simpson<- diversity(abun_table,index="simpson", MARGIN = 2)
richness<-specnumber(abun_table, MARGIN = 2)
evenness <- shannon / log(richness)
index <- as.data.frame(cbind(shannon, simpson, richness, evenness))
write.csv(index,'Protist.quantitative.diversity.index.csv', quote = F)

#Beta Diversity
df <- t(abun_table)
distance <- as.matrix(vegdist(df,method = "bray"))
write.csv(distance,file="Protist.quantitative.bray.csv")
distance <- as.matrix(vegdist(df,method = "jaccard",binary = T))
write.csv(distance,file="Protist.quantitative.jaccard.csv")