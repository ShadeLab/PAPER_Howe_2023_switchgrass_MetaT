# gapseq analysis

GSoutput=read_delim('../gapseqOutput.txt')
GSkey=read_delim("../gapseqKey.txt") 
names(GSkey)[3] = 'pathway'
magTax=readRDS('MAG.gtdbtk.taxonomy.RDS')

magTax_mod=magTax %>% 
  mutate(mag=if_else(str_detect(otu, 'misc'), "M", "S"),
         mag=paste(mag, str_split(.$otu, "[.]", simplify = T)[,3], sep=''))

gapseqSummary=GSoutput %>%
  filter(Name %in% GSkey$pathway) %>%
  left_join(., GSkey, by = c( 'Name' = 'pathway')) %>%
  pivot_longer(cols=c(2:42)) %>%
  left_join(., magTax_mod, by = c('name' = 'mag'))

gapseqSummary %>%
  group_by(Phylum, Family, Genus) %>%
  filter(name != 'M22') %>%
  mutate(mag = factor(name, levels=unique(name))) %>%
  ggplot(aes(pathway_summary, mag, col=group, size=as.factor(value))) +
  geom_point(pch=20) +
  theme_bw() +
  scale_color_npg()+
  scale_size_manual(values = c(-1,4), breaks = c(FALSE, TRUE))+ 
  facet_grid(~group, scales = 'free_x',space = 'free_x') +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank()
        ) +
  labs(x=NULL, y=NULL)
  

