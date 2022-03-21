library(IRanges)
library(tidyverse)
library(ggsci)
# 
# GAPSEQfiles=list.files("gapseq_outdir_MetaCyc/", pattern = "all-Pathways.tbl") %>%
#   str_remove(., "-all-Pathways.tbl")
# 
# combinedGAPSEQ=c()
# for (mag in GAPSEQfiles){
#   
#   # import pathway output
#   headers_path = read.csv(paste0('gapseq_outdir_MetaCyc/',mag,'-all-Pathways.tbl'), 
#                      skip = 2, header = F, nrows = 1, as.is = T, sep = '\t')
#   df_path = read.csv(paste0('gapseq_outdir_MetaCyc/',mag,'-all-Pathways.tbl'), 
#                      skip = 3, header = F, sep='\t')
#   colnames(df_path)= headers_path
#   
#   #import reaction output
#   headers_react = read.csv(paste0('gapseq_outdir_MetaCyc/',mag,'-all-Reactions.tbl'), 
#                      skip = 2, header = F, nrows = 1, as.is = T, sep = '\t')
#   df_reac = read.csv(paste0('gapseq_outdir_MetaCyc/',mag,'-all-Reactions.tbl'), 
#                      skip = 3, header = F, sep='\t')
#   colnames(df_reac)= headers_react
#   
#   #merge pathways with reaction output
#   merged_df=left_join(df_path, df_reac, by = c("ID" = "pathway")) %>% 
#     #filter(Prediction == "true") %>%
#     mutate(MAG=paste(str_split(mag, '[.]', simplify = T)[,1],
#                      str_split(mag, '[.]', simplify = T)[,3], sep='') %>% 
#              str_replace(., 'misc', 'M') %>% str_replace(., 'switch', 'S')
#            )
#   combinedGAPSEQ=rbind(combinedGAPSEQ,merged_df)
# }

#saveRDS(combinedGAPSEQ,'inputData/gapseq.combined.rds')

# Import all datasets needed for the anlaysis
combinedGAPSEQ=readRDS('inputData/gapseq.combined.rds')
metatCov=read.delim('inputData/orfs-normalized-median-metat-no-fungal-no-ann.txt')
map.metag=readRDS('inputData/metaG.map.full.RDS')        
magAnnot=read.delim('inputData/fullMAG.annotation.tsv', header = T)
magTax=readRDS('inputData/MAG.gtdbtk.taxonomy.RDS')
gapseqKey=read.delim("inputData/ranked_gapseq_pathways_metacyc_UPDATED.txt")
gapseq_dram_hits=readRDS('inputData/gapseq_dram_hits.rds')

magAnnot$X = NULL

magAnnot=magAnnot%>%mutate(MAG=paste(str_split(magAnnot$fasta, '[.]', simplify = T)[,1],
                                     str_split(magAnnot$fasta, '[.]', simplify = T)[,3], sep='')) %>%
  mutate(MAG=str_replace(MAG, 'misc', 'M'), 
         MAG=str_replace(MAG, 'switch', 'S'))

# Filter out only best hits.
# Certain ORF have multiple calls and for linking this with transcript coverage
# we need to get only one hit. The hits will be called by the overlapped genomic location.
# For example if gapseq called a gene between locations 1 and 1200 on a contig A 
# and on the same contig DRAM annotated 2 genes, first from 1-800 and next one
# from 900 - 1400, gapseq will be assigned these two genes since transcript tables
# are based on prodigal/DRAM annotations.

# filtDF=combinedGAPSEQ %>%
#   dplyr::mutate(eval = as.numeric(evalue)) %>%
#   dplyr::filter(!is.na(eval)) %>%
#   dplyr::group_by(MAG,stitle, ID,Name,name,ec) %>%
#   dplyr::mutate(bestMin=min(evalue)) %>%
#   ungroup() %>%
#   dplyr::group_by(MAG,stitle, ID,Name,name,ec) %>%
#   filter(eval == bestMin) %>%
#   ungroup() %>%
#   mutate(newStart=if_else(sstart>send, send, sstart),
#          end_position=if_else(send<sstart, sstart, send))
# 
# testDF=filtDF[filtDF$Prediction == 'true',]
# names(testDF)[31]='newEnd'
# testDF$seqLeng=NULL

# gapseq_dram_hits=c()
# 
# for(contig in unique(testDF$stitle)){
#   #contig = 'k127_6459063_M_bin.1'
#   targetContig=testDF[testDF$stitle == contig,]
# 
#   for(i in 1:length(rownames(targetContig))
#       ){
#     
#     gs1 = with(targetContig[i,], IRanges(newStart, newEnd))
#     dram1 = with(magAnnot[magAnnot$scaffold == contig,], IRanges(start_position , end_position))
#     overlap = countOverlaps(dram1, gs1) #!= 0
#     dram_short=magAnnot[magAnnot$scaffold == contig,]
#     
#     if(sum(overlap)>0){
#       
#       dram.hits=dram_short[overlap>0,]
#       bothTrue=cbind(targetContig[i,],dram.hits)
#       bothTrue$lineNoGAPSEQ=i
#       gapseq_dram_hits = rbind(gapseq_dram_hits,bothTrue)
#       
#     }
#     
#     else{
#       
#       dram.hits=dram_short[overlap>0,]
#       dram.hits[1,] = NA
#       bothTrue=cbind(targetContig[i,],dram.hits)
#       bothTrue$lineNoGAPSEQ=i
#       gapseq_dram_hits = rbind(gapseq_dram_hits,bothTrue)
#       
#     }
# 
#   }
# }
# # Result of the loop are all gapseq hits linked to a ORF from DRAM annotations. 
# # If no ORF was detected then NAs are added.

#saveRDS(gapseq_dram_hits, 'inputData/gapseq_dram_hits.rds')

gapseq_dram_hits=gapseq_dram_hits[-27]
gapseq_dram_transcript_hits = 
gapseq_dram_hits %>% mutate(orf_column=paste(scaffold, gene_position, sep='_'),
                            activity=if_else(orf_column %in% rownames(metatCov), 
                                             'active', 'inactive'))
#gapseq_dram_transcript_hits %>% group_by(MAG,name, activity) %>% summarise(nCount=length(activity))

gapseqShort=gapseq_dram_transcript_hits[c(56,58,59, 2, 10)] %>%
  left_join(gapseqKey)

gapseq_activity_fig=gapseqShort %>%
  filter(!is.na(group),
         !(group %in% c('amino acid biosynthesis',
                      'amino acid degradation',
                      'central pathway', 'respiration'))) %>%
  mutate(groupings=if_else(fun %in% c('stress','phytohormone', 'plant'), fun, "general")) %>%
  ggplot(aes(pathway_summary, y=MAG, 
             col=factor(groupings, levels=c('stress','phytohormone', 'plant', 'general')), 
             shape=activity)) +
  geom_point(size=2) +
  theme_bw() +
  scale_color_lancet()+
  scale_shape_manual(values=c(19,1))+
  #scale_size_manual(values = c(-1,4), breaks = c(FALSE, TRUE))+ 
  facet_grid(~group, scales = 'free_x',space = 'free_x') +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        #strip.text.x = element_blank()
  ) +
  labs(x=NULL, y=NULL, col=NULL)


gapseqFinds=gapseqShort %>%
  filter(!is.na(group),
         !(group %in% c('amino acid biosynthesis',
                        'amino acid degradation',
                        'central pathway', 'respiration'))) %>%
  mutate(groupings=if_else(fun %in% c('stress','phytohormone', 'plant'), fun, "general"),
         tool="gapseq",
         PA=if_else(activity == 'active', 1,0)) %>%
  group_by(MAG, groupings, group, pathway_summary, tool) %>%
  summarise(activity=if_else(sum(PA)>0, 'active', 'inactive'))


