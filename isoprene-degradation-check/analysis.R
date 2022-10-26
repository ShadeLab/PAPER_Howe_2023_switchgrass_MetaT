library(stringr)
library(tidyr)
library(dplyr)

# Load metaT ORF abundances
load(file="melted_phy3.RData") #loads object ps, filtered
head(ps)

# Blast results for the three groups of genes
setwd('~/Box Sync/Papers/Phyllosphere-Function/review-response-isoprene/')


# Isoprene Genes from other bacteria
# Isoprene protein blast from Varivox 
blast <- read.delim(file='orfs.x.varivorax.blastnout.best', header=FALSE, sep="\t")
sub <- subset(blast, V11 < 1e-5)
split1 <- str_split_fixed(sub$V1, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, 'bin.', 2)[,2]
bin = paste(split2a, split4, sep='')
sub$bin <- bin
test <- sub %>% group_by(bin, V2) %>% dplyr::summarise(count=n())
test2 <- test %>% spread(bin, count)
#write.table(test2, file='iso-rhodo.blast.txt', sep='\t', quote=FALSE)
#test3 <- data.frame(test2)
#rownames(test3) <- test3$V2
#test3$V2 <- NULL
#test3[test3 > 0] = 1
#test3[is.na(test3)] = 0
foo <- subset(ps, OTU %in% sub$V1)
foo2 <- merge(foo, sub, by.x = "OTU", by.y = "V1")
all <- foo2 %>% select(Sample, year, OTU, Abundance, bin, month, year, V2)
f <- ddply(all, .(Sample, bin, year, month, V2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(bin, year, month, V2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9")
f2$bin <- factor(f2$bin, levels=c("M1", "M8", "M9", "M21",  "M44",  "M52",   "M60",  "M66",  "M67","M77",   "M86",  "M87",  "M94", "M102",  "M109", "M111", "S27",  "S28",  "S29",   "S36", "S56",  "S61", "S74"))
levels(f2$bin) <- c("M1_Actino", "M8_Rhizob", "M9_Cytoph", "M21_Burk",  "M44_Rhizob",  "M52_Rhizob",   "M60_Actino",  "M66_Burk",  "M67_Actino","M77_Burk",   "M86_Actino",  "M87_Rizhob",  "M94_Burk", "M102_Rhizob",  "M109_Sphingo", "M111_Rhizob", "S27_Burk",  "S28_Pseudo",  "S29_Actino",   "S36_Pseudo", "S56_Sphingo",  "S61_Burk", "S74_Pseudo")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2$year = as.factor(f2$year)
f2$zero_flag = "0"
f2[f2$MEAN != 0,]$zero_flag = "1"
p = ggplot(f2, aes_string(x="month", y="MEAN"))
p+geom_point(stat="identity", size=3, shape = 21, aes(fill=zero_flag)) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+scale_fill_manual(values=c("white","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_grid(year~bin)+theme(strip.text.x = element_text(angle = 90))
#MAGS_associated_with_varivorax.png


# Isoprene Induced Genes
blast <- read.delim(file='all+orfs.x.isometab.blastxout.best', header=FALSE, sep="\t")
#sub <- subset(blast, V3 > 30 & V4 > 300)
sub <- subset(blast, V11 < 1e-5)
split1 <- str_split_fixed(sub$V1, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, 'bin.', 2)[,2]
bin = paste(split2a, split4, sep='')
sub$bin <- bin
test <- sub %>% group_by(bin, V2) %>% dplyr::summarise(count=n())
test2 <- test %>% spread(bin, count)
#write.table(test2, file='isometab.blast.txt', sep='\t', quote=FALSE)
#test3 <- data.frame(test2)
#rownames(test3) <- test3$V2
#test3$V2 <- NULL
#test3[test3 > 0] = 1
#test3[is.na(test3)] = 0
foo <- subset(ps, OTU %in% sub$V1)
foo2 <- merge(foo, sub, by.x = "OTU", by.y = "V1")
all <- foo2 %>% select(Sample, year, OTU, Abundance, bin, month, year, V2)
f <- ddply(all, .(Sample, bin, year, month, V2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(bin, year, month, V2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9")
f2$bin <- factor(f2$bin, levels=c("M1", "M8", "M9", "M12",  "M17",  "M21",  "M32",  "M35",  "M44",  "M47",  "M52",  "M55",  "M60",  "M66",  "M67","M77",   "M86",  "M87",  "M94",  "M99",   "M100", "M102", "M105", "M109", "M111", "S8", "S9", "S27",  "S28",  "S29",  "S30",  "S36",  "S50",  "S56",  "S61",  "S71", "S74",   "S80",   "S117", "S120"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2$year = as.factor(f2$year)
f2$month = as.factor(f2$month)
p = ggplot(f2, aes_string(x="V2", y="MEAN",  shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Genes Previously Upregulated in Response to Isoprene') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(bin~month)#+theme(
    strip.text.x = element_blank()
  )
#MAGs_upregulated_isoprene_genes
blast <- read.delim(file='all_orfs.x.isopreneinduced.blastxout.best', header=FALSE, sep="\t")
sub <- subset(blast, V11 < 1e-5)
split1 <- str_split_fixed(sub$V1, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, 'bin.', 2)[,2]
bin = paste(split2a, split4, sep='')
sub$bin <- bin
test <- sub %>% group_by(bin, V2) %>% dplyr::summarise(count=n())
test2 <- test %>% spread(bin, count)
#write.table(test2, file='isoinduced.blast.txt', sep='\t', quote=FALSE)
foo <- subset(ps, OTU %in% sub$V1)
foo2 <- merge(foo, sub, by.x = "OTU", by.y = "V1")
all <- foo2 %>% select(Sample, year, OTU, Abundance, bin, month, year, V2)
f <- ddply(all, .(Sample, bin, year, month, V2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(bin, year, month, V2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9")
f2$bin <- factor(f2$bin, levels=c("M1", "M8", "M9", "M12",  "M17",  "M21",  "M32",  "M35",  "M44",  "M47",  "M52",  "M55",  "M60",  "M66",  "M67","M77",   "M86",  "M87",  "M94",  "M99",   "M100", "M102", "M105", "M109", "M111", "S8", "S9", "S27",  "S28",  "S29",  "S30",  "S36",  "S50",  "S56",  "S61",  "S71", "S74",   "S80",   "S117", "S120"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2_2016 <- subset(f2, year == "2016")
p = ggplot(f2_2016, aes_string(x="month", y="MEAN", color="bin"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~V2)
p = ggplot(f2_2016, aes_string(x="V2", y="MEAN", color="bin"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(bin~month)
f2_2017 <- subset(f2, year == "2017")
p = ggplot(f2_2017, aes_string(x="month", y="MEAN", color="bin"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~V2)
p = ggplot(f2_2017, aes_string(x="V2", y="MEAN", color="bin"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(bin~month)
f2$year = as.factor(f2$year)
p = ggplot(f2, aes_string(x="V2", y="MEAN", color="bin", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(bin~month)
# 
# blast <- read.delim(file='all_orf.x.methylcitrate.blastxout.best', header=FALSE, sep="\t")
# sub <- subset(blast, V11 < 1e-5)
# split1 <- str_split_fixed(sub$V1, "_", 3)[,3]
# split2a <- str_split_fixed(split1, '_', 2)[,1]
# split2b <- str_split_fixed(split1, '_', 2)[,2]
# split3 <- str_split_fixed(split2b, '_', 2)[,1]
# split4 <- str_split_fixed(split3, 'bin.', 2)[,2]
# bin = paste(split2a, split4, sep='')
# sub$bin <- bin
# test <- sub %>% group_by(bin, V2) %>% dplyr::summarise(count=n())
# test2 <- test %>% spread(bin, count)
# write.table(test2, file='methylcitrate.blast.txt', sep='\t', quote=FALSE)
# foo <- subset(ps, OTU %in% sub$V1)
# foo2 <- merge(foo, sub, by.x = "OTU", by.y = "V1")
# all <- foo2 %>% select(Sample, year, OTU, Abundance, bin, month, year, V2)
# f <- ddply(all, .(Sample, bin, year, month, V2), summarise, SUM=sum(Abundance))
# f2 <- ddply(f, .(bin, year, month, V2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
# temp = subset(f2, month == "9")
# f2$bin <- factor(f2$bin, levels=c("M1", "M8", "M9", "M12",  "M17",  "M21",  "M32",  "M35",  "M44",  "M47",  "M52",  "M55",  "M60",  "M66",  "M67","M77",   "M86",  "M87",  "M94",  "M99",   "M100", "M102", "M105", "M109", "M111", "S8", "S9", "S27",  "S28",  "S29",  "S30",  "S36",  "S50",  "S56",  "S61",  "S71", "S74",   "S80",   "S117", "S120"))
# limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
# f2_2016 <- subset(f2, year == "2016")
# p = ggplot(f2_2016, aes_string(x="month", y="MEAN", color="bin"))
# p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
#   xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~V2)
# p = ggplot(f2_2016, aes_string(x="V2", y="MEAN", color="bin"))
# p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
#   xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(bin~month)
# f2_2017 <- subset(f2, year == "2017")
# p = ggplot(f2_2017, aes_string(x="month", y="MEAN", color="bin"))
# p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
#   xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~V2)
# p = ggplot(f2_2017, aes_string(x="V2", y="MEAN", color="bin"))
# p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
#   xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(bin~month)
# f2$year = as.factor(f2$year)
# p = ggplot(f2, aes_string(x="V2", y="MEAN", color="bin", shape="year"))
# p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
#   xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(bin~month)


# Statistical differences?
blast <- read.delim(file='all+orfs.x.isometab.blastxout.best', header=FALSE, sep="\t")
#sub <- subset(blast, V3 > 30 & V4 > 300)
sub <- subset(blast, V11 < 1e-5)
split1 <- str_split_fixed(sub$V1, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, 'bin.', 2)[,2]
bin = paste(split2a, split4, sep='')
sub$bin <- bin
test <- sub %>% group_by(bin, V2) %>% dplyr::summarise(count=n())
test2 <- test %>% spread(bin, count)
#write.table(test2, file='isometab.blast.txt', sep='\t', quote=FALSE)
#test3 <- data.frame(test2)
#rownames(test3) <- test3$V2
#test3$V2 <- NULL
#test3[test3 > 0] = 1
#test3[is.na(test3)] = 0
foo <- subset(ps, OTU %in% sub$V1)
foo2 <- merge(foo, sub, by.x = "OTU", by.y = "V1")
all <- foo2 %>% select(Sample, year, OTU, Abundance, bin, month, year, V2)
f <- ddply(all, .(Sample, bin, year, month, V2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(bin, year, month, V2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))

head(f)
f_summed = ddply(f, .(Sample, bin, year, month), summarise, SUM2=sum(SUM))
f2 <- ddply(f_summed, .(bin, year, month), summarise, MEAN=mean(SUM2), SE=sd(SUM2)/sqrt(length(SUM2)))
f2$year <- as.factor(f2$year)
f2$zero_flag = "0"
f2[f2$MEAN != 0,]$zero_flag = "1"
p = ggplot(f2, aes_string(x="month", y="MEAN", color="year"))
p+geom_point(stat="identity", size=3, shape=21, aes(fill=zero_flag)) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+scale_fill_manual(values=c("white","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~bin)

ggplot(f_summed, aes(x=SUM2)) + geom_density() #not normal
f2 <- f_summed %>% mutate(month2 = if_else(month <= 6, "early", "late"))
#kruskal-wallis test 
ss <- unique(f2$bin)
ss_stats <- data.frame(matrix(NA, nrow = length(ss), ncol = 2))
rownames(ss_stats) <- ss
for (i in 1:length(ss)){
  ss_data <- subset(f2, bin == ss[i])
  res <- kruskal.test(SUM2 ~ month2, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}

# This is the ISO genes blast table from the BlastX
blast <- read.delim(file='orfs.x.isoproteins.blastxout.best.descr.out', header=FALSE, sep="\t")
sub <- blast #already filtered for evalues
split1 <- str_split_fixed(sub$V1, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, 'bin.', 2)[,2]
bin = paste(split2a, split4, sep='')
sub$bin <- bin
foo <- subset(ps, OTU %in% sub$V1)
foo2 <- merge(foo, sub, by.x = "OTU", by.y = "V1")
all <- foo2 %>% select(Sample, year, OTU, Abundance, bin, month, year, V4)
f <- ddply(all, .(Sample, bin, year, month, V4), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(bin, year, month, V4), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9")
f2$bin <- factor(f2$bin, levels=c("M1", "M8", "M9", "M12",  "M17",  "M21",  "M32",  "M35",  "M44",  "M47",  "M52",  "M55",  "M60",  "M66",  "M67","M77",   "M86",  "M87",  "M94",  "M99",   "M100", "M102", "M105", "M109", "M111", "S8", "S9", "S27",  "S28",  "S29",  "S30",  "S36",  "S50",  "S56",  "S61",  "S71", "S74",   "S80",   "S117", "S120"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2$year = as.factor(f2$year)
f2$MEAN2 = log(f2$MEAN)
p = ggplot(f2, aes_string(x="month", y="MEAN", color="V4", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(bin~.)

