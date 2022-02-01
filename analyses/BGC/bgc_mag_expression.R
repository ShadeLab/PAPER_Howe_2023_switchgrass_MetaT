library(tidyverse)

bgc_pos=read_delim('bgcPosition.tsv')
bgc_pos_orfNo=read_delim('bgcPosition_orfNo.txt') %>% mutate(contig_orf=paste(contig,gene_position, sep='_'))

head(bgc_pos)
magAnnot=readRDS('full.mag.annotation.RDS')
magAnnot_short= magAnnot[,c('scaffold', 'gene_position', 'start_position', 'endqstat -u stopnise
                            _position')]

magAnnot_short_bgc=magAnnot_short[magAnnot_short$scaffold %in% unique(bgc_pos_orfNo$contig),]


metatCov=read.delim('orfs-normalized-median-metat-no-fungal-no-ann (1).txt')
bgc_class=read.delim('antiSMASH_Annotations_Full.tsv')

map.metag=readRDS('metaG.map.full.RDS')
HKGlist.mags=read.table('list-of-hkg2.txt',header=F)

head(HKGlist.mags)


longDF.metat=data.frame(contigID=rownames(metatCov), metatCov) %>%
  gather(name, relCov,-contigID)

splitDF=as.data.frame(str_split(longDF.metat$contigID, '_', simplify = T)[,-5])
longDF.metat$contigs_only=paste(splitDF$V1,splitDF$V2,splitDF$V3,splitDF$V4, sep='_')
longDF.metat$bin=paste(splitDF$V3, splitDF$V4, sep='_')

longDF.metat$HKG=if_else(longDF.metat$contigID %in% unique(HKGlist.mags$V1), 'hkg', '')
longDF.metat$BGC=if_else(longDF.metat$contigID %in% unique(bgc_pos_orfNo$contig_orf), 'bgc', '')

filt.metat=longDF.metat %>% filter(HKG == 'hkg' | BGC == 'bgc')

plotDF=filt.metat %>%
  mutate(date=str_split(name, '_', simplify = T)[,3],
         year=if_else(str_detect(date, '2017'), 2017, 2016),
         month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
  mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
  mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
  mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
  mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
  left_join(bgc_class, by = c('contigs_only' = 'contig'))

hkgSumm=plotDF %>%
  dplyr::group_by(bin.x) %>%
  filter(HKG == 'hkg') %>%
  dplyr::summarise(nHKGs= length(unique(contigID)))
#' All 41 MAGs with 22 and up to 69 HKG per MAG

bgcSumm=plotDF %>%
  dplyr::group_by(bin.x, ) %>%
  filter(BGC.x == 'bgc') %>%
  dplyr::summarise(nBGCs_orfs= length(unique(contigID)),
                   nBGC_classes=length(unique(BiG.SCAPE_class)))
#' 39 MAGs contain BGCs. With up to 5 different classes (M_bin.8 and S_bin.80).
#' S_bin.28 with most BGC orfs (n=158). Taking with a grain of salt since for 
#' some BGC ORFs we can not align with those predicted by prodigal.

plotDF= plotDF %>%
  mutate(tagH=if_else(HKG == 'hkg', 1, 0),
         tagB=if_else(BGC.x == 'bgc', 2, 0),
         tagH_B=tagH+tagB,
         tag=if_else(tagB == 2,  'bgc', 'hkg'),
         tag_full=if_else(tagH_B == 3, 'bgc_hkg',tag)) 

# Classified as HKG and BGC:
BGC_HKG_orfs=unique(plotDF$contigID[plotDF$tag_full == 'bgc_hkg'])

HKG_meanCov_time=plotDF %>%
  filter(HKG == 'hkg') %>%
  group_by(bin.x, name) %>%
  dplyr::summarise(meanHKG_time=mean(relCov))       
         
plotDF %>%
  filter(BGC.x == 'bgc', 
         relCov > 0) %>%
  left_join(HKG_meanCov_time) %>%
  group_by(bin.x, name, contigID, tag_full, BiG.SCAPE_class, month, year) %>%
  dplyr::summarise(ratioB_H=log2(mean(relCov)/meanHKG_time)) %>%
  filter(str_detect(bin.x,'S_')) %>%
  ggplot(aes(x=as.factor(month), y=ratioB_H, col=BiG.SCAPE_class, 
             group=BiG.SCAPE_class, alpha=ratioB_H)) +
  theme_bw()+
  geom_hline(yintercept=0, col='grey60') +
  geom_hline(yintercept=2, linetype="dashed", col='grey60') +
  geom_hline(yintercept=-2, linetype="dashed", col='grey60') +
  scale_y_continuous(breaks = seq(-6, 12, by = 4)) +
  geom_jitter() +
  facet_grid(year~bin.x) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x='Month', y='log2 ratio between BGC and mean HKG', 
       alpha=NULL,col='Big SCAPE class')


# Include information about the BGC core location
bgcCore=read.delim('bgcCorePosition.txt') %>% mutate(contig_end=paste(contig, endCore, sep='_'))
bgc_pos_orfNo$contig_end=paste(bgc_pos_orfNo$contig, bgc_pos_orfNo$end, sep='_')

bgc_pos_orfNo$Core= if_else(bgc_pos_orfNo$contig_end %in% bgcCore$contig_end, 'core', '')

coreBGClist=bgc_pos_orfNo$contig_orf[bgc_pos_orfNo$Core == 'core']

head(plotDF)

plotDF %>%
  filter(BGC.x == 'bgc', 
         relCov > 0) %>%
  left_join(HKG_meanCov_time) %>%
  mutate(coreBGC=if_else(contigID %in% coreBGClist, 'core', 'other')) %>%
  group_by(bin.x, contigID, tag_full, BiG.SCAPE_class, month, year,coreBGC) %>%
  dplyr::summarise(ratioB_H=log2(mean(relCov)/meanHKG_time),
                   p_value = t.test(unlist(relCov), unlist(meanHKG_time))$p.value) %>%
  filter(str_detect(bin.x,'S_')) %>%
  mutate(outline=if_else(coreBGC == 'core' & ratioB_H>2, 'yes', 'no')) %>%
  ggplot(aes(x=as.factor(month), y=ratioB_H, #fill=BiG.SCAPE_class, 
             group=BiG.SCAPE_class, alpha=ratioB_H,
             size=coreBGC)) +
  geom_jitter(aes(fill=BiG.SCAPE_class, color=outline), pch=21) +
  theme_bw()+
  geom_hline(yintercept=0, col='grey60') +
  geom_hline(yintercept=2, linetype="dashed", col='grey60') +
  geom_hline(yintercept=-2, linetype="dashed", col='grey60') +
  scale_y_continuous(breaks = seq(-6, 12, by = 4)) +
  scale_color_manual(values = c('transparent','black')) +
  scale_size_manual(values=c(3,.5))+
  facet_grid(year~bin.x) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x='Month', y='log2 ratio between BGC and mean HKG', 
       alpha=NULL,fill='Big SCAPE class', color='Upregulated BGC')

samples2keep=unique(plotDF$name[plotDF$BGC.x == 'bgc'])


