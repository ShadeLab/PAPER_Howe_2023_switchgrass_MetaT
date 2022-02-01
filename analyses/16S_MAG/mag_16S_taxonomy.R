library(tidyverse)
library(RColorBrewer)

phyllo_16_core=unique(c(tempCore,tempCore_sw17,tempCore_m16))
saveRDS(phyllo_16_core, '16S.core.taxa.list.RDS')
phyllo_16_core=readRDS('16S.core.taxa.list.RDS')

phyllo_core_taxonomy=tax_filtered[rownames(tax_filtered) %in% phyllo_16_core,]
saveRDS(phyllo_core_taxonomy, '16S.core.taxa.taxonomy.RDS')
phyllo_core_taxonomy=readRDS('16S.core.taxa.taxonomy.RDS')

magTax=as.data.frame(str_split(glbrc_gtdbtk$classification, ";", simplify = T))
names(magTax)=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
rownames(magTax)=glbrc_gtdbtk$user_genome

magTax=magTax %>%
  mutate(Kingdom=str_remove(Kingdom, 'd__'),
         Phylum=str_remove(Phylum, 'p__'),
         Class=str_remove(Class, 'c__'),
         Order=str_remove(Order, 'o__'),
         Family=str_remove(Family, 'f__'),
         Genus=str_remove(Genus, 'g__'),
         Species=str_remove(Species, 's__'))

saveRDS(magTax, 'MAG.gtdbtk.taxonomy.RDS')
magTax=readRDS('MAG.gtdbtk.taxonomy.RDS')

mag16Score=left_join(magTax, phyllo_core_taxonomy[,c(1,7)])
length(unique(mag16Score$otu)) 
# 33 OTUs taxonomically (at Genus level) similar to 21 MAGs.

summary(magTax$Genus %in% phyllo_core_taxonomy$Genus)

saveRDS(mag16Score, 'mag_16S_core_taxonomy.RDS')
mag16Score=readRDS('mag_16S_core_taxonomy.RDS')

magTax$dataset='metagenomics'
phyllo_core_taxonomy$dataset='amplicon'

magTax$otu=rownames(magTax)
names(magTax)
names(phyllo_core_taxonomy)

magTax=magTax %>% mutate(Genus=str_remove(Genus, "_A|_E"))
magTax$membership=if_else(magTax$Genus %in% phyllo_core_taxonomy$Genus, 'shared', 'other')

phyllo_core_taxonomy$membership=if_else(phyllo_core_taxonomy$Genus %in% magTax$Genus, 'shared', 'other')

longDF=rbind(magTax[-10,], phyllo_core_taxonomy[,c(1:8,11,12)])

longDF=longDF %>% 
  mutate(Genus=if_else(is.na(Genus), "", Genus))

plotDF_tax=longDF %>%
  group_by(Genus, membership, dataset) %>%
  summarise(nFeatures=length(unique(otu))) %>%
  mutate(genusNew=if_else(membership == 'shared', Genus, 'other'))

plotDF_tax_small=plotDF_tax[,c(3:5)]
tidyr::spread(plotDF_tax, dataset, nFeatures)

unique(plotDF_tax$genusNew)
plotDF_tax=plotDF_tax%>%
  mutate(genusNew= forcats::fct_relevel(genusNew,
                               'other', 'Acidovorax', "Chryseobacterium",
                               "Hymenobacter", "Massilia","Methylobacterium",
                               "Novosphingobium","Pedobacter" ,"Pseudomonas", 
                               "Quadrisphaera","Sphingomonas","Spirosoma"))
  
library("ggsci")
#library(wesanderson)
library(RColorBrewer)
colorCode=c('grey70',brewer.pal(n = 12, name = "Paired")[-11])
colorCoreCode=colorCode
names(colorCoreCode)=c('other', 'Acidovorax', "Chryseobacterium",
"Hymenobacter", "Massilia","Methylobacterium",
"Novosphingobium","Pedobacter" ,"Pseudomonas", 
"Quadrisphaera","Sphingomonas","Spirosoma")

ggplot()+
  geom_bar(data = plotDF_tax,
           aes(x = dataset, y =nFeatures, fill = genusNew),
           width =0.3, stat="identity")+
  theme_bw()+
  #scale_fill_jco()+
  scale_fill_manual(values = colorCode) +
  labs(y='# features', x=NULL, fill=NULL)
#ggsave('OneDrive - Michigan State University/MSU/phyllosphere/recentMAGdata/newFig3B.eps', 
#       device = "eps", width = 95, height=100, units = 'mm')

#ggsave('fig3B_noMAG22.eps', 
#       device = "eps", width = 95, height=100, units = 'mm')




