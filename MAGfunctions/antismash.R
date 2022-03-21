library(tidyverse)

#Import datasets
bgc_pos=read_delim('inputData/bgcPosition.tsv')
bgc_pos_orfNo=read_delim('inputData/bgcPosition_orfNo.txt') %>% mutate(contig_orf=paste(contig,gene_position, sep='_'))

magAnnot=read.delim('inputData/fullMAG.annotation.tsv', header = T)

magAnnot=magAnnot%>%mutate(MAG=paste(str_split(magAnnot$fasta, '[.]', simplify = T)[,1],
                                     str_split(magAnnot$fasta, '[.]', simplify = T)[,3], sep='')) %>%
  mutate(MAG=str_replace(MAG, 'misc', 'M'), 
         MAG=str_replace(MAG, 'switch', 'S'),
         orf_id=paste(scaffold, gene_position, sep='_'))

magAnnot_short= magAnnot[,c('scaffold', 'orf_id','gene_position', 'start_position', 'end_position')]

metatCov=read.delim('inputData/orfs-normalized-median-metat-no-fungal-no-ann.txt')
bgc_class=read.delim('inputData/antiSMASH_Annotations_Full.tsv')

map.metag=readRDS('inputData/metaG.map.full.RDS')
HKGlist.mags=read.table('inputData/list-of-hkg2.txt',header=F)

magTax=readRDS('inputData/MAG.gtdbtk.taxonomy.RDS')

bgcCore=read.delim('inputData/bgcCorePosition.txt') %>% mutate(contig_end=paste(contig, endCore, sep='_'))


## Start of the analysis:
longDF.metat=data.frame(contigID=rownames(metatCov), metatCov) %>%
  gather(name, relCov,-contigID)

splitDF=as.data.frame(str_split(longDF.metat$contigID, '_', simplify = T)[,-5])
longDF.metat$contigs_only=paste(splitDF$V1,splitDF$V2,splitDF$V3,splitDF$V4, sep='_')
longDF.metat$bin=paste(splitDF$V3, splitDF$V4, sep='_')

bgc.metat= longDF.metat %>% filter(contigID %in% unique(bgc_pos_orfNo$contig_orf))

plotDF=bgc.metat %>%
  mutate(date=str_split(name, '_', simplify = T)[,3],
         year=if_else(str_detect(date, '2017'), 2017, 2016),
         month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
  mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
  mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
  mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
  mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
  left_join(bgc_class, by = c('contigs_only' = 'contig')) %>% 
  left_join(magAnnot_short, by = c('contigID' = 'orf_id'))
  

bgcActivity=plotDF %>% group_by(contigID, mag, Product_Prediction, BiG.SCAPE_class) %>%
  summarise(activity=if_else(sum(relCov)>0, 'active', 'inactive'))

# plotDF %>%
#   filter(BGC.x == 'bgc', 
#          relCov > 0) %>%
#   left_join(HKG_meanCov_time) %>%
#   group_by(bin.x, name, contigID, tag_full, BiG.SCAPE_class, month, year) %>%
#   dplyr::summarise(ratioB_H=log2(mean(relCov)/meanHKG_time)) %>%
#   filter(str_detect(bin.x,'S_')) %>%
#   ggplot(aes(x=as.factor(month), y=ratioB_H, col=BiG.SCAPE_class, 
#              group=BiG.SCAPE_class, alpha=ratioB_H)) +
#   theme_bw()+
#   geom_hline(yintercept=0, col='grey60') +
#   geom_hline(yintercept=2, linetype="dashed", col='grey60') +
#   geom_hline(yintercept=-2, linetype="dashed", col='grey60') +
#   scale_y_continuous(breaks = seq(-6, 12, by = 4)) +
#   geom_jitter() +
#   facet_grid(year~bin.x) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   labs(x='Month', y='log2 ratio between BGC and mean HKG', 
#        alpha=NULL,col='Big SCAPE class')

# Find only BGC's biosynthethic genes and use those for plotting

# bgc_dram_hits=c()
# 
# for(contig in unique(bgcCore$contig)){
#   #contig = 'k127_10410083_M_bin.8'
#   targetContig=bgcCore[bgcCore$contig == contig,]
#   
#   for(i in 1:length(rownames(targetContig))){
#     
#     bgc1 = with(targetContig[i,], IRanges(startCore, endCore))
#     dram1 = with(magAnnot[magAnnot$scaffold == contig,], IRanges(start_position , end_position))
#     overlap = countOverlaps(dram1, bgc1) #!= 0
#     dram_short=magAnnot[magAnnot$scaffold == contig,]
#     
#     if(sum(overlap)>0){
#       
#       dram.hits=dram_short[overlap>0,]
#       bothTrue=cbind(targetContig[i,],dram.hits)
#       bothTrue$lineNoGAPSEQ=i
#       bgc_dram_hits = rbind(bgc_dram_hits,bothTrue)
#       
#     }
#     
#     else{
#       
#       dram.hits=dram_short[overlap>0,]
#       dram.hits[1,] = NA
#       bothTrue=cbind(targetContig[i,],dram.hits)
#       bothTrue$lineNoGAPSEQ=i
#       bgc_dram_hits = rbind(bgc_dram_hits,bothTrue)
#       
#     }
#     
#   }
# }

bgc_dram_hits=readRDS('inputData/antismash_coreGene_dram_matchup.rds')
antiSMASH_dram=left_join(bgc_dram_hits, bgc_class) %>%
  mutate(activity=if_else(orf_id %in% rownames(metatCov), 'active', 'inactive'))


antiSMASH_dram_tax_activity=left_join(antiSMASH_dram, magTax, by=c('fasta' = 'otu'))

antiSMASH_dram_tax_activity %>%
  mutate(PA=if_else(activity == 'active', 1, 0)) %>%
  group_by(contig_end, mag, Product_Prediction, BiG.SCAPE_class) %>%
  summarise(active=if_else(sum(PA)>0, 'active', 'inactive')) %>%
  ggplot(aes(x=Product_Prediction, y=mag, 
             fill=BiG.SCAPE_class, 
             shape=active)) +
  geom_point(size=2) +
  theme_bw()+
  scale_color_simpsons() +
  scale_shape_manual(values=c(19,1))+
  facet_grid(~BiG.SCAPE_class, scales = 'free_x',space = 'free_x') +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        strip.background = element_blank())+
  labs(x=NULL, y=NULL)

antismashFinds=antiSMASH_dram_tax_activity %>%
  mutate(PA=if_else(activity == 'active', 1, 0),
         tool='antiSMASH', 
         MAG=mag, 
         groupings=NA,
         group=BiG.SCAPE_class,
         functions=Product_Prediction) %>%
  group_by(MAG,groupings,group,functions, tool) %>%
  summarise(activity=if_else(sum(PA)>0, 'active', 'inactive'))


