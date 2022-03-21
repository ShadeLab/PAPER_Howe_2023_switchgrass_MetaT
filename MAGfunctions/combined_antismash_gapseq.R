magTax=readRDS('MAG.gtdbtk.taxonomy.RDS')
magTax_mod=magTax %>% 
  mutate(mag=if_else(str_detect(otu, 'misc'), "M", "S"),
         mag=paste(mag, str_split(.$otu, "[.]", simplify = T)[,3], sep=''))
names(magTax_mod)[11]='MAG'

MAGfunctions=rbind(gapseqFinds, manualFinds)

names(MAGfunctions)[4]='functions'
MAGfunctions=rbind(MAGfunctions,antismashFinds) %>%
  left_join(magTax_mod[,c(1:7,10,11)]) %>%
  mutate(taxonomy=paste(Phylum, Class, Order, Family, Genus, sep='.')) %>%
  group_by(taxonomy)%>% 
  mutate(newMAG=factor(MAG, levels=unique(MAG)))

MAGfunctions$groupings[is.na(MAGfunctions$groupings)]='other'

MAGfunctions %>%
  filter(MAG != 'M22') %>%
  ggplot(aes(functions, y=newMAG, 
           col=factor(groupings, levels=c('stress','phytohormone', 'plant', 'general', 'other')), 
           shape=activity)) +
  geom_point(size=2) +
  theme_bw() +
  scale_color_npg()+
  scale_shape_manual(values=c(19,1))+
  facet_grid(~tool+group, scales = 'free_x',space = 'free_x') +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        #strip.text.x = element_blank()
  ) +
  labs(x=NULL, y=NULL, col=NULL)
  

