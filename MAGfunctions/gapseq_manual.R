metatCov=read.delim('inputData/orfs-normalized-median-metat-no-fungal-no-ann.txt')
map.metag=readRDS('inputData/metaG.map.full.RDS')
map.metag.short=map.metag[c(5,6,15,39)] 
magAnnot=read.delim('inputData/fullMAG.annotation.tsv', header = T)
magAnnot=magAnnot%>%mutate(MAG=paste(str_split(magAnnot$fasta, '[.]', simplify = T)[,1],
                                     str_split(magAnnot$fasta, '[.]', simplify = T)[,3], sep='')) %>%
  mutate(MAG=str_replace(MAG, 'misc', 'M'), 
         MAG=str_replace(MAG, 'switch', 'S'),
         orf_id=paste(scaffold, gene_position, sep='_'))

magAnnot_short= magAnnot[,c('scaffold', 'orf_id','gene_position', 'start_position', 'end_position')]


#' #' Color palette for all MAGs
#' magAnnot %>%
#'   filter(MAG %in% MAGS_gapseq) %>%
#'   group_by(order, genus) %>%
#'   summarise(nMags=length(unique(MAG))) %>%
#'   arrange(desc(order))
#' 
#' colActino <- colorRampPalette(c("#52006A", "white")) #purple
#' colActino=colActino(10)[1:6] 
#' names(colActino)=unique(magAnnot$MAG[magAnnot$order == 'Actinomycetales' & magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colBurk=colorRampPalette(c("#DD4A48", "white")) # red 
#' colBurk=colBurk(10)[1:6] 
#' names(colBurk)=unique(magAnnot$MAG[magAnnot$order == 'Burkholderiales'& magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colRhiz=colorRampPalette(c("#FF7600", "white")) # orange
#' colRhiz=colRhiz(10)[1:6]
#' names(colRhiz)=unique(magAnnot$MAG[magAnnot$order == 'Rhizobiales'& magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colSphingom=colorRampPalette(c("#0DA907", "white")) # green
#' colSphingom=colSphingom(10)[1:5]
#' names(colSphingom)=unique(magAnnot$MAG[magAnnot$order == 'Sphingomonadales'& magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colPseu=colorRampPalette(c("#005792", "white")) #blue
#' colPseu=colPseu(10)[1:3]
#' names(colPseu)=unique(magAnnot$MAG[magAnnot$order == 'Pseudomonadales' & magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colCyto=colorRampPalette(c("#FFE400", "white")) #yellow
#' colCyto=colCyto(10)[1:2]
#' names(colCyto)=unique(magAnnot$MAG[magAnnot$order == 'Cytophagales' & magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colFlavo=colorRampPalette(c("darkgray", "white")) #gray
#' colFlavo=colFlavo(10)[1:2]
#' names(colFlavo)=unique(magAnnot$MAG[magAnnot$order == 'Flavobacteriales' & magAnnot$MAG %in% MAGS_gapseq])
#' 
#' colRick=colorRampPalette(c("#632626", "white")) #brown
#' colRick=colRick(10)[1:2]
#' names(colRick)=unique(magAnnot$MAG[magAnnot$order == 'Rickettsiales' & magAnnot$MAG %in% MAGS_gapseq])
#' 
#' oneRepMAGs=c("Enterobacterales", "Mycobacteriales","Propionibacteriales","Sphingobacteriales","XYA12-FULL-58-9")
#' onecolor=c('#FF4646','#FFEEAD','black','#AACDBE','#064635')
#' names(onecolor)=unique(magAnnot$MAG[magAnnot$order %in% oneRepMAGs & magAnnot$MAG %in% MAGS_gapseq ])
#' 
#' MAGpalette=c(colRick, colFlavo, colCyto,colPseu, colSphingom, colRhiz, colBurk, colActino, onecolor)
#' library(scales)
#' show_col(MAGpalette)

#'Pathways/genes involved in specific functions

#' Phytohormone homeostasis

#' ACC deaminase (1-aminocyclopropane-1-carboxylate deaminase)
accDeam=magAnnot %>%
  filter(str_detect(kegg_hit, 'carboxylate deaminase'))
accDeam$groupings = 'phytohormone'
accDeam$group= "specialized functions"
accDeam$pathway_summary= "ACC deaminase"
accDeam$tool='manual'

# accDeamPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% accDeam$orf_id) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#   geom_jitter() +
#   scale_color_manual(values = MAGpalette) +
#   theme_bw() +
#   xlim(4.5,9.5) +
#   scale_shape_manual(values = c(16,1)) +
#   theme(legend.key.size = unit(.4, 'cm')) +
#   labs(x='Sampling month', y='HKG normalized median coverage', title = 'ACC deaminase') 

#'IAA degradation
iaa = magAnnot %>%
  filter(str_detect(kegg_hit, 'aldehyde dehydrogenase')) %>%
  mutate(
    groupings = 'phytohormone',
    group= "specialized functions",
    pathway_summary= "IAA degradation",
    tool='manual'
  )

# iaaPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% iaa$orf_id) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#   geom_jitter(size=2) +
#   theme_bw() +
#   xlim(4.5,9.5)+
#   scale_color_manual(values = MAGpalette) +
#   scale_shape_manual(values = c(16,1)) +
#   theme(legend.key.size = unit(.4, 'cm')) +
#   labs(x='Sampling month', y='HKG normalized median coverage', title = 'IAA degradation')

#'Salicylic acid production and degradation
#'
# sa_genes=gapseq_annot[gapseq_annot$Name == 'methylsalicylate degradation' & 
#                         gapseq_annot$Prediction == 'true',]

SA=magAnnot %>%
  filter(str_detect(kegg_hit, 'salicylate'))  %>%
  mutate(
    groupings = 'phytohormone',
    group= "specialized functions",
    pathway_summary= "SA degradation",
    tool='manual'
  )

# #no transcriptional activity
# saPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in%unique(SA$orf_id)) 
#   #%>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#   geom_jitter(size=2) +
#   theme_bw() +
#   scale_color_manual(values = MAGpalette) +
#   scale_shape_manual(values = c(16,1)) +
#   labs(x='Sampling month', y='HKG normalized median coverage', title = 'SA degradation')



#' Secretion systems important in plant-microbe
#' https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-9-S1-S2
typeI=magAnnot %>%
  filter(str_detect(kegg_hit, 'type I secretion system'))%>%
  mutate(type='type I')
typeII=magAnnot %>%
  filter(str_detect(kegg_hit, 'type II secretion system'))%>%
  mutate(type='type II')
typeIII=magAnnot %>%
  filter(str_detect(kegg_hit, 'type III secretion system')) %>%
  mutate(type='type III') #type III secretion system (T3SS), are employed for effector translocation into the host plants
typeIV=magAnnot %>%
  filter(str_detect(kegg_hit, 'type IV secretion system'))%>%
  mutate(type='type IV')
typeVI=magAnnot %>%
  filter(str_detect(kegg_hit, 'type VI secretion system')) %>%
  mutate(type='type VI') #Type VI Secretion System (T6SS) contractile nanoweapons allows bacteria to inject toxins

secretionSystems=rbind(typeI, typeII, typeIII, typeIV, typeVI)%>%
  mutate(
    groupings = 'stress',
    group= "Secretion systems",
    pathway_summary= type,
    tool='manual'
  )

# SSplot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% secretionSystems$orf_id) %>%
#   left_join(secretionSystems, by= c('contigID'='orf_id')) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#   geom_jitter(size=2, alpha=.6) +
#   scale_color_manual(values = MAGpalette) +
#   scale_shape_manual(values = c(16,1)) +
#   theme_bw() +
#   xlim(4.5,9.5)+
#   theme(legend.key.size = unit(.4, 'cm')) +
#   labs(x='Sampling month', y='HKG normalized median coverage', title = 'Bacterial secretion systems') +
#   facet_grid(~type)

#' Oxalate degradation
oxalate=magAnnot %>%
  filter(str_detect(kegg_hit, 'oxalate'))  %>%
  mutate(
    groupings = 'plant',
    group= "specialized functions",
    pathway_summary= 'oxalate degradation',
    tool='manual'
  )

# oxalPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% oxalate$orf_id) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(month, y=relCov, col=MAG, shape=as.factor(year))) +
#   geom_jitter(size=2, alpha=.6) +
#   theme_bw() +
#   theme(legend.key.size = unit(.4, 'cm')) +
#   xlim(4.5,9.5)+
#   scale_color_manual(values = MAGpalette) +
#   scale_shape_manual(values = c(16,1)) +
#   labs(x='Sampling month', y='HKG normalized median coverage', title = 'Oxalate metabolism')



#' #' Osmotic stress 
#' #Glycine betaine
#' glyBet=magAnnot %>%
#'   filter(str_detect(kegg_hit, 'betaine'))
#' glyBet=gapseq_annot[str_detect(gapseq_annot$Name, 'betaine') &
#'                        gapseq_annot$Prediction == 'true',] %>% filter(!is.na(orf_id))
#' 
#' glyPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#'   gather(name, relCov,-contigID) %>%
#'   filter(contigID %in% glyBet$orf_id) %>%
#'   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#'   mutate(date=str_split(name, '_', simplify = T)[,3],
#'          year=if_else(str_detect(date, '2017'), 2017, 2016),
#'          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#'   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#'   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#'                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#'   filter(relCov>0) %>%
#'   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#'   geom_jitter(size=2, alpha=.6) +
#'   theme_bw() +
#'   xlim(4.5, 9.5)+
#'   theme(legend.key.size = unit(.4, 'cm')) +
#'   scale_color_manual(values = MAGpalette) +
#'   scale_shape_manual(values = c(16,1)) +
#'   labs(x='Sampling month', y='HKG normalized median coverage', title = 'Glycine betaine') 
#' 
#' 
#' #Taurine
#' taur=magAnnot %>%
#'   filter(str_detect(kegg_hit, 'taurine'))
#' 
#' tauPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#'   gather(name, relCov,-contigID) %>%
#'   filter(contigID %in% taur$orf_id) %>%
#'   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#'   mutate(date=str_split(name, '_', simplify = T)[,3],
#'          year=if_else(str_detect(date, '2017'), 2017, 2016),
#'          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#'   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#'   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#'                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#'   filter(relCov>0) %>%
#'   ggplot(aes(factor(month), y=relCov, col=MAG, shape=factor(year))) +
#'   geom_jitter(size=2, alpha=.6) +
#'   theme_bw() +
#'   scale_color_manual(values = MAGpalette) +
#'   scale_shape_manual(values = c(16,1)) +
#'   labs(x='Sampling month', shape=NULL,
#'        y='HKG normalized median coverage', title = 'Taurine') 
#' 
#' 
#' #choline
#' chol=magAnnot %>%
#'   filter(str_detect(kegg_hit, 'choline'))
#' choline=gapseq_annot[str_detect(gapseq_annot$Name, 'choline degradation') &
#'                        gapseq_annot$Prediction == 'true',] %>% filter(!is.na(orf_id))
#' 
#' choPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#'   gather(name, relCov,-contigID) %>%
#'   filter(contigID %in% choline$orf_id) %>%
#'   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#'   mutate(date=str_split(name, '_', simplify = T)[,3],
#'          year=if_else(str_detect(date, '2017'), 2017, 2016),
#'          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#'   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#'   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#'                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#'   filter(relCov>0) %>%
#'   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#'   geom_jitter(size=2, alpha=.6) +
#'   theme_bw() +
#'   xlim(4.5,9.5)+
#'   theme(legend.key.size = unit(.4, 'cm')) +
#'   scale_color_manual(values = MAGpalette) +
#'   scale_shape_manual(values = c(16,1)) +
#'   labs(x='Sampling month', y='HKG normalized median coverage', title = 'Choline') 
#' 
#' library(ggpubr)
#' ggarrange(glyPlot, tauPlot, choPlot)
#' 
#' #' Xylose utilization
#' xylose=magAnnot %>%
#'   filter(str_detect(kegg_hit, 'xylose'))
#' xylose=gapseq_annot[str_detect(gapseq_annot$Name, 'xylitol') &
#'                        gapseq_annot$Prediction == 'true',] %>% filter(!is.na(orf_id))
#' xyloPlot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#'   gather(name, relCov,-contigID) %>%
#'   filter(contigID %in% xylose$orf_id) %>%
#'   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#'   mutate(date=str_split(name, '_', simplify = T)[,3],
#'          year=if_else(str_detect(date, '2017'), 2017, 2016),
#'          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#'   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#'   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#'   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#'                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#'   filter(relCov>0) %>%
#'   ggplot(aes(month, y=relCov, col=MAG, shape=factor(year))) +
#'   geom_jitter(size=2, alpha=.6) +
#'   theme_bw() +
#'   xlim(4.5, 9.5)+
#'   theme(legend.key.size = unit(.4, 'cm')) +
#'   scale_color_manual(values = MAGpalette) +
#'   scale_shape_manual(values = c(16,1)) +
#'   labs(x='Sampling month', y='HKG normalized median coverage', title = "Xylitole metabolism") 

#' CO oxidation!
#' https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.15770

CO=magAnnot %>%
  filter(str_detect(kegg_hit, 'monoxide dehydrogenase')) %>%
  mutate(
    groupings = 'general',
    group= "specialized functions",
    pathway_summary= 'CO oxydation',
    tool='manual'
  )

# COplot=data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% CO$orf_id) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(factor(month), y=relCov, col=MAG, shape=factor(year))) +
#   geom_jitter(size=2, alpha=.6) +
#   theme_bw() +
#   scale_color_manual(values = MAGpalette) +
#   scale_shape_manual(values = c(16,1)) +
#   labs(x='Sampling month', y='HKG normalized median coverage', title="CO oxidation") 

# # Nitrate reduction
# nitratRed1=magAnnot %>%
#   filter(str_detect(kegg_hit, 'nitrate reductase')) %>%
#   mutate(geneName='nitrate reductase') 
# nitratRed2=magAnnot %>%
#   filter(str_detect(kegg_hit,'formate dehydrogenase')) %>%
#   mutate(geneName='formate dehydrogenase')
# 
# nitRed=rbind(nitratRed1, nitratRed2)
# 
# data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% nitRed$orf_id) %>%
#   left_join(nitRed, by = c('contigID' = 'orf_id')) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0,
#          geneName=='nitrate reductase') %>%
#   ggplot(aes(factor(month), y=relCov, col=MAG)) +
#   geom_jitter(size=2, alpha=.6) +
#   theme_bw() +
#   labs(x='Sampling month', y='HKG normalized median coverage', title="CO oxidation") +
#   facet_grid(geneName~year)


# Hydrogenases (H2 oxydation)
hydr=magAnnot %>%
    filter(str_detect(kegg_hit, ' hydrogenase')) %>%
  mutate(
    groupings = 'general',
    group= "specialized functions",
    pathway_summary= 'H2 oxydation',
    tool='manual'
  )

# data.frame(contigID=rownames(metatCov), metatCov) %>%
#   gather(name, relCov,-contigID) %>%
#   filter(contigID %in% fol$orf_id) %>%
#   mutate(plant=if_else(str_detect(name, 'G5'), 'swtichgrass', 'miscanthus')) %>%
#   mutate(date=str_split(name, '_', simplify = T)[,3],
#          year=if_else(str_detect(date, '2017'), 2017, 2016),
#          month= if_else(str_detect(date, 'JUN'), 6, 0)) %>%
#   mutate(month= if_else(str_detect(date, 'JUL'), 7, month)) %>%
#   mutate(month= if_else(str_detect(date, 'SEP'), 9, month)) %>%
#   mutate(month= if_else(str_detect(date, 'AUG'), 8, month)) %>%
#   mutate(month= if_else(str_detect(date, 'MAY'), 5, month)) %>%
#   mutate(MAG=paste(str_split(contigID, '_', simplify = T)[,3],
#                    str_split(contigID, '_', simplify = T)[,4], sep='') %>% str_remove(., 'bin.')) %>%
#   filter(relCov>0) %>%
#   ggplot(aes(factor(month), y=relCov, col=MAG)) +
#   geom_jitter(size=2, alpha=.6) +
#   theme_bw() +
#   labs(x='Sampling month', y='HKG normalized median coverage', title="glutamate dehydrogenase") +
#   facet_grid(~year)

secretionSystems$type=NULL

manualFinds=rbind(accDeam,iaa,SA,oxalate, CO, hydr, secretionSystems) %>%
  mutate(PA=if_else(orf_id %in% rownames(metatCov), 1, 0)) %>% #presence of transcripts
  group_by(MAG, groupings, group, pathway_summary, tool) %>%
  summarise(activity=if_else(sum(PA)>0, 'active', 'inactive'))
