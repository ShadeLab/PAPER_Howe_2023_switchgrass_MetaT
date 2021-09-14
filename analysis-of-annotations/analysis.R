library(tidyr)
library(dplyr)
library(plyr)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(corrplot)
library(Hmisc)
library(stringr)

# Set local workign directory where files are located
setwd("~/Box Sync/Papers/Phyllosphere-Function/Code/")

# Input Files needed
# sample-meta.txt = metadata describing samples
# checkm/checkm_output.txt = output from checkM describing MAGs
# melted_phy_average_coverage_mags_per_metag.rds = R data object, summary of mapped reads to MAGs
# mags-in-sg-and-misc2.csv = detection of mags in sg and misc
# read-count.txt = total reads in metagenomes
# transcipts.txt = total reads in metats
# summary-med-read-cover-meta-orfs.cleaned.txt = transcripts abundances to mags
# Meta3.txt = metadata for metatranscriptomes
# mag.x.refsoil.blastpout.columns.l4.plus.finalann2 = functional annotations of transcripts
# parsed-bins-subsystems.mappings3.txt = annotations of functional roles 
# hits-against-iowa-misc2.txt= transcripts hits in metagenomes from other perennials
# others-metags2.tsv = transcripts hits in metagenomes from other perennials
# enriched.txt = enriched functions in MAGs

# Metadata - Sample and MAGs associated data read in
meta <- read.delim(sep="\t", file="./sample-meta.txt", header=FALSE, strip.white=TRUE, row.names = 1)
colnames(meta) <- c("plot", "fert", "date", "id", "plant", "date2", "month")
meta$year <- format(as.Date(meta$date, format="%d-%b-%Y"), "%Y")
# CheckM Annotations with Selected Columns and filter by completeness and contamination
checkm <- read.delim(sep='\t', file="./checkm (adina@iastate.edu)/checkm_output.txt", header = TRUE, strip.white=TRUE)
checkm_qc <- subset(checkm, completeness > 70 & contamination < 20)
list_of_interest <- checkm_qc$plantbin

# Read in Average Estimated Median Read Coverage per MAG (i.e., plantbin)
cov <- readRDS('melted_phy_average_coverage_mags_per_metag.rds')
# Proportion of Samples that are Associated with MAGs, selecting bins associated with > 10% in each plant type
association <- read.csv(file="mags-in-sg-and-misc2.csv", na.strings=c("","NA"))
association$Association2 <- ifelse(association$M1 >= 0.1 & association$S1 >= 0.1, "B", "NB")
list_of_bins <- (subset(association, Association2 == "B"))$plantbin
select_association <- cov[cov$plantbin %in% list_of_interest,]
select_association2 <- select_association[select_association$plantbin %in% list_of_bins,]
length(unique(select_association2$plantbin))

# Creating the phyloseq MAGs in metagenomes object, phy
f1 <- select(select_association2, plantbin, MEAN, Sample)
f2 <- f1 %>% spread(Sample,MEAN)
rownames(f2) <- f2$plantbin
f2$plantbin<-NULL
og_names <- rownames(f2)
# Renaming MAGs based on CheckM to Paper Final Names
new_names <- c("MAG11_Acti",	"MAG5_Rhiz",	"MAG10_Burk", 	"MAG1_Rick",	"MAG16_Cyto", 	"MAG13_Unkn", 	"MAG4_Rhiz", 	"MAG15_Acti", 	"MAG8_Pseudo",	"MAG9_Acti", 	"MAG3_Rhiz",	"MAG12_Sphi",	"MAG2_Burk",	"MAG7_Pseudo", 	"MAG14_Acti",	"MAG6_Burk")
abundance <- otu_table(as.matrix(f2), taxa_are_rows=TRUE)
metadata <- sample_data(meta)
rownames(f2) <- mapvalues(rownames(f2), from=rownames(f2), to=c("MAG11_Acti",	"MAG5_Rhiz",	"MAG10_Burk", 	"MAG1_Rick",	"MAG16_Cyto", 	"MAG13_Unkn", 	"MAG4_Rhiz", 	"MAG15_Acti", 	"MAG8_Pseudo",	"MAG9_Acti", 	"MAG3_Rhiz",	"MAG12_Sphi",	"MAG2_Burk",	"MAG7_Pseudo", 	"MAG14_Acti",	"MAG6_Burk"))
abundance <- otu_table(as.matrix(f2), taxa_are_rows=TRUE)
metadata <- sample_data(meta)
phy <- phyloseq(metadata, abundance)
melted_phy <- psmelt(phy)

# Number of Metagenomes each MAG is represented
f <- ddply(melted_phy, .(OTU, Sample, plant), summarise, MEAN = mean(Abundance), SUM = length(Sample))
f$presence <- ifelse(f$MEAN > 0, 1, 0)
f2 <- ddply(f, .(OTU, plant), summarise, sum_samples = sum(presence))

# Normalize by abundances by total number of library metag reads, phy_norm 
counts <- read.csv("read-count.txt", sep="\t", header=FALSE)
head(counts)
colnames(counts) <- c("fullname", "counts", "samplename2", "date", "Sample")
counts$fullname <- NULL
counts <- unique(counts)
abundance_melt <- melt(abundance)
colnames(abundance_melt) <- c("plantbin", "Sample", "Abundance")
abund_merge <- merge(abundance_melt, counts, by = "Sample")
abund_merge$norm <- abund_merge$Abundance/abund_merge$counts
f1 <- select(abund_merge, plantbin, Sample, norm)
f2 <- f1 %>% spread(Sample, norm)
rownames(f2) <- f2$plantbin
f2$plantbin<-NULL
abundance <- otu_table(as.matrix(f2), taxa_are_rows=TRUE)
phy_norm <- phyloseq(metadata, abundance)
melted_phy <- psmelt(phy_norm)

# Figure production of Heatmap of MAGs per plant type [Figure 2A]
phy_heatmap <- prune_samples(sample_sums(phy)>0, phy)
otu_table(phy_heatmap) <- round(otu_table(phy_heatmap), 3)
samples_unordered <- rownames(sample_data(phy_heatmap))
order_of_samps <- order(as.Date(sample_data(phy_heatmap)$date,  format="%d-%B-%Y"))
sample_data(phy_heatmap)$plant <- str_replace(sample_data(phy_heatmap)$plant, "S", "Switchgrass")
sample_data(phy_heatmap)$plant <- str_replace(sample_data(phy_heatmap)$plant, "M", "Miscanthus")
taxa_order = c('MAG1_Rick',	'MAG2_Burk',	'MAG3_Rhiz',	'MAG4_Rhiz', 'MAG5_Rhiz',	'MAG6_Burk',	'MAG7_Pseudo',	'MAG8_Pseudo',	'MAG10_Burk',	'MAG11_Acti',	'MAG9_Acti',	'MAG12_Sphi',	'MAG13_Unkn',	'MAG14_Acti',	'MAG15_Acti',	'MAG16_Cyto')
p=plot_heatmap(phy_heatmap,  taxa.order = taxa_order, sample.order=samples_unordered[order_of_samps], sample.label = "month", low="#66CCFF", high="#000033", na.value="white")
p$data<- p$data[order(p$data$date),]
p$scales$scales[[1]]$name <- "Date of Sampling"
p$scales$scales[[2]]$name <- "Abundance"
p + theme(text = element_text(size=12), axis.text.x = element_text(angle = 90, hjust=0, size=10)) +ylab('') + facet_grid(~plant, scales="free")+scale_color_gradient(limits=c(0,4.0),expression(Abundance))
ggsave(file="Figure_2A.png", plot = last_plot(), width = 20, height = 10)


# Figure production of plot of MAGs over time in metagenomes [Figure S2]
colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
f <- ddply(melted_phy, .(OTU, plant, month), summarise, MEAN = mean(Abundance), SE=sd(Abundance)/sqrt(length(Abundance)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f$month <- factor(f$month, levels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
f_temp <- subset(f, plant == "S" & month == "Oct")
#f$OTU <- factor(f$OTU, levels=f_temp$OTU[order(-f_temp$MEAN)])
f$OTU <- factor(f$OTU, levels=taxa_order)
f$zero_flag = 0
f[f$MEAN != 0,]$zero_flag= 1
f$zero_flag <- as.character(f$zero_flag)
f$plant <- str_replace(f$plant, "S", "Switchgrass")
f$plant <-  str_replace(f$plant, "M", "Miscanthus")
p = ggplot(f, aes_string(x="OTU", y="MEAN", color="month"))
p+geom_point(stat="identity", size=2, shape=21, aes(fill=zero_flag))+theme_bw()+geom_errorbar(limits, width=0)+
  facet_grid(~plant)+theme(text=element_text(size=16))+scale_shape(solid=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust=0, size=10))+
  xlab("MAG ID") + ylab("Relative Abundance") + scale_fill_manual(values=c("white","black"))+
  scale_color_manual(values=colorBlindBlack8[1:7]) 
ggsave(file="Figure_S1.png", plot = last_plot(), width = 10, height = 10)

# Figure production of plot of MAGs over time in metagenomes [Figure S1A]
colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
f <- ddply(melted_phy, .(OTU, plant, month), summarise, MEAN = mean(Abundance), SE=sd(Abundance)/sqrt(length(Abundance)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f$month <- factor(f$month, levels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
f_temp <- subset(f, plant == "S" & month == "Oct")
#f$OTU <- factor(f$OTU, levels=f_temp$OTU[order(-f_temp$MEAN)])
f$OTU <- factor(f$OTU, levels=taxa_order)
f$zero_flag = 0
f[f$MEAN != 0,]$zero_flag= 1
f$zero_flag <- as.character(f$zero_flag)
p = ggplot(f, aes_string(x="OTU", y="MEAN", color="month"))
p+geom_point(stat="identity", size=2, shape=21, aes(fill=zero_flag))+theme_bw()+geom_errorbar(limits, width=0)+
  facet_grid(~plant)+theme(text=element_text(size=16))+scale_shape(solid=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust=0, size=10))+
  xlab("MAG ID") + ylab("Relative Abundance") + scale_fill_manual(values=c("white","black"))+
  scale_color_manual(values=colorBlindBlack8[1:7]) 
ggsave(file="Figure_S1A.png", plot = last_plot(), width = 10, height = 10)

# Table production of functional roles in MAGs vs RefSoil
mag_ss <- read.csv("parsed-bins-subsystems.mappings3.txt", sep="\t")
rownames(mag_ss )<-mag_ss$X
mag_ss$X <- NULL
mag_ss_norm<-data.frame(t(apply(mag_ss, 1, function(x) x / sum(x) )))
rownames(mag_ss_norm) <- gsub("_", "", rownames(mag_ss_norm))
#mag_ss_norm <- mag_ss_norm[rownames(mag_ss_norm) %in% core,]
refsoil_ss <- read.csv("all-refsoil-ss.txt", sep="\t")
rownames(refsoil_ss)<-refsoil_ss$X
refsoil_ss$X <- NULL
refsoil_ss2<-data.frame(t(apply(refsoil_ss, 1, function(x) x / sum(x) )))
merged_ss <- rbind(mag_ss_norm, refsoil_ss2)
merged_ss$X <- rownames(merged_ss)
merged_ss_melt <- melt(merged_ss)
meta1 <- data.frame(X=rownames(mag_ss_norm))
meta1$meta <- meta1$X
meta2 <- data.frame(X=rownames(refsoil_ss2), meta = "refsoil", meta2="refsoil")
meta1$meta2 = "phyllo"
meta_all <- rbind(meta1, meta2)
merged_gen_ann <- merge(merged_ss_melt, meta_all, by="X")
f_ann <-ddply(merged_gen_ann, .(variable, meta2), summarise, MEAN = mean(value), SE=sd(value)/sqrt(length(value)))
f_ann <- subset(f_ann, variable != "NA.")
f_ann <- subset(f_ann, variable != "Clustering.based.subsystems")
test <- f_ann
test$variable2 <- gsub("[.]", " ", test$variable)
f_ann$variable <- test$variable2
f_ann2 <- subset(f_ann, meta2 == "phyllo")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
temp <- f_ann2[order(-f_ann2$MEAN),]
f_ann$variable <- factor(f_ann$variable, levels=temp$variable)

# Statistical comparison of RefSoil vs MAG functions
library("ggpubr")
ggdensity(f_ann$MEAN) #not normal
head(merged_gen_ann)
#kruskal-wallis test 
ss <- unique(merged_gen_ann$variable)
ss_stats <- data.frame(matrix(NA, nrow = length(ss), ncol = 2))
rownames(ss_stats) <- ss
for (i in 1:length(ss)){
  ss_data <- subset(merged_gen_ann, variable == ss[i])
  res <- kruskal.test(value ~ meta2, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}
ss_stats$variable <- rownames(ss_stats)
colnames(ss_stats) <- c("H-stat", "p-value", "variable")
f_ann <-ddply(merged_gen_ann, .(variable, meta2), summarise, MEAN = mean(value), SE=sd(value)/sqrt(length(value)))
maxes <- ddply(f_ann, .(variable), summarise, Value = max(MEAN))
f_ann_max <- merge(f_ann, maxes, by = "variable")
f_ann_max$max_ss <- ifelse(f_ann_max$MEAN == f_ann_max$Value, 1, 0)
f_ann_max2 <- merge(f_ann_max, ss_stats, by = "variable")
stats_summary <- f_ann_max2 %>% select(variable, meta2, MEAN, "p-value") %>% spread(meta2, MEAN)
write.table(stats_summary, file="Table_S3.txt", sep="\t", quote=FALSE, row.names = FALSE)

#MetaTranscriptomes Analysis
x <- read.csv(file="summary-med-read-cover-meta-orfs.cleaned.txt", sep='\t')
dim(x)
x2 <- x[1:60641,1:74]
colnames(x2) = colnames(x)[2:length(colnames(x))]
abundance <- otu_table(as.matrix(x2), taxa_are_rows=TRUE)
meta <- read.delim(sep="\t", file="./meta3.txt", header=FALSE, strip.white=TRUE, row.names = 1)
colnames(meta) <- c( "plot", "rep", "fert", "date", "id", "date2", "month", "year")
metadata <- sample_data(meta)
transcripts <- phyloseq(metadata, abundance)
transcripts_ss <- prune_samples(sample_sums(transcripts)>0, transcripts)
transcripts_ss2 <- filter_taxa(transcripts_ss, function(x) sum(x) > 0, TRUE)
samples_unordered <- rownames(sample_data(transcripts_ss2))
order_of_samps <- order(as.Date(sample_data(transcripts_ss2)$date,  format="%d-%B-%Y"))

# Data wrangling of merging transcript abundances and standardizing by total transcripts
annotations_mags = read.csv("mag.x.refsoil.blastpout.columns.l4.plus.finalann2", sep="\t", header=FALSE, row.names=1)
head(annotations_mags)
head(og_names)
head(new_names)
head(annotations_mags$new_V6)
annotation_mags2 <- subset(annotations_mags, annotations_mags$V5 != "NA")
foo2<-as.matrix(annotation_mags2)
annotation <- tax_table(foo2)
trans_ann <- phyloseq(metadata, abundance, annotation)
trans_ann2 <- prune_samples(sample_sums(trans_ann)>0, trans_ann)
trans_ann2 <- filter_taxa(trans_ann2, function(x) sum(x) > 0, TRUE)
trans_ann2
trans_melt <- psmelt(trans_ann2)
head(trans_melt)
counts <- read.csv("transcripts.txt", sep=" ", header=FALSE)
colnames(counts) <- c("Sample", "counts")
head(counts)
ps_merged <- merge(trans_melt, counts, by = "Sample")
ps_merged$norm <- ps_merged$Abundance/ps_merged$counts #use this table
ps_merged$new_V6 <- mapvalues(ps_merged$V6, from=og_names, to=new_names)
ps_merged$V6 = ps_merged$new_V6


#Figure comparing functions
f=ddply(ps_merged, .(date, V5, Sample), summarise, SUM = sum(norm))
f2 = ddply(f, .(date, V5), summarise, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
f3 <- subset(f2, V5 %in% c("DNA Metabolism",  "RNA Metabolism"))
f3$date <- factor(f3$date, levels = c("31-May-16", "12-Jul-16", "12-Sep-16", "15-May-17","5-Jun-17", "26-Jun-17", "17-Jul-17", "7-Aug-17", "28-Aug-17", "18-Sep-17"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3, aes_string(x="date", y="MEAN", color="V5"))
p+geom_point(stat="identity", size=3)+scale_shape(solid=FALSE)+theme_bw()+geom_errorbar(limits, width=0)+
    xlab('Date of Sample') + ylab('Relateve Abundance of Transcripts')+theme_bw()+ labs(colour = "Functional Class")
ss_data <- subset(ps_merged, V5 %in%  c("DNA Metabolism",  "RNA Metabolism"))
f = ddply(ss_data, .(date, V5, Sample), summarise, SUM = sum(norm))
date_list <- unique(f$date)
ss_stats <- data.frame(matrix(NA, nrow = length(date_list), ncol = 2))
rownames(ss_stats) <- date_list
for (i in 1:length(date_list)){
  ss_data <- subset(f, date == date_list[i])
  res <- kruskal.test(SUM ~ V5, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}
write.table(ss_stats, file="dna_vs_rna.txt", sep="\t", quote=FALSE)

#Figure production of heatmap of switchgrass metaT transcripts mapped to MAGs over the two years [Figure 2B]
f=ddply(ps_merged, .(date, V6, Sample), summarise, SUM = sum(norm))
head(metadata)
abund_heatmap <- f %>% select(Sample, V6, SUM) %>% spread(Sample, SUM)
rownames(abund_heatmap) <- abund_heatmap$V6
abund_heatmap$V6 <- NULL
abund_heatmap2 <- otu_table(as.matrix(abund_heatmap), taxa_are_rows=TRUE)
heatmap_phy <- phyloseq(metadata, abund_heatmap2)
samples_unordered <- rownames(sample_data(heatmap_phy))
order_of_samps <- order(as.Date(sample_data(heatmap_phy)$date,  format="%d-%B-%Y"))
taxa_order = c('MAG1_Rick',	'MAG2_Burk',	'MAG3_Rhiz',	'MAG4_Rhiz', 'MAG5_Rhiz',	'MAG6_Burk',	'MAG7_Pseudo',	'MAG8_Pseudo',	'MAG10_Burk',	'MAG11_Acti',	'MAG9_Acti',	'MAG12_Sphi',	'MAG13_Unkn',	'MAG14_Acti',	'MAG15_Acti',	'MAG16_Cyto')
sample_data(heatmap_phy)$year <- str_replace(sample_data(heatmap_phy)$year, "16", "Switchgrass 2016")
sample_data(heatmap_phy)$year <- str_replace(sample_data(heatmap_phy)$year, "17", "Switchgrass 2017")
p=plot_heatmap(heatmap_phy,  taxa.order=taxa_order, sample.order=samples_unordered[order_of_samps], sample.label = "month", low="#66CCFF", high="#000033", na.value="white")
p$data$date <- as.Date(p$data$date, format="%d%B%Y")
p$data<- p$data[order(p$data$date),]
p = p + theme(axis.text.x = element_text(angle = 90, hjust=0, size=6)) + xlab("Sampling Month") + ylab("ORFs in MAG")+facet_grid(~year, scales="free_x")
p$scales$scales[[1]]$name <- "Date of Sample"
p + ylab('MAG ID')
ggsave(file="Figure_2B.png", plot = last_plot(), width = 20, height = 10)

# Figure production of plot of MAGs over time in metagenomes [Figure S1B]
colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ps_merged_subset <- ps_merged
f=ddply(ps_merged_subset, .(year, month, V5, V3, V4, V6, Sample), summarise, SUM = sum(norm))
f2=ddply(f, .(month, V6, Sample, year), summarise, SUM2 = sum(SUM))
f3=ddply(f2, .(month, V6, year), summarise, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))
f3$V6 <- factor(f3$V6, levels=taxa_order)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3$zero_flag = 0
f3[f3$MEAN != 0,]$zero_flag= 1
f3$zero_flag <- as.character(f3$zero_flag)
f3$year <- str_replace(f3$year, "16", "Switchgrass 2016")
f3$year <- str_replace(f3$year, "17", "Switchgrass 2017")
f3$month <- factor(f3$month, levels = c("May", "Jun", "Jul", "Aug", "Sep"))
p = ggplot(f3, aes_string(x="V6", y="MEAN", color="month"))
p+geom_point(stat="identity", size=3, shape=21, aes(fill=zero_flag))+scale_shape(solid=FALSE)+theme_bw()+geom_errorbar(limits, width=0)+
  theme(text=element_text(size=12))+theme(axis.text.x = element_text(angle = 90, hjust=1, size=12))+theme(strip.text.x = element_text(size = 10))+
  xlab("MAG ID") + ylab("Transcript Abundances") + scale_color_manual(values=colorBlindBlack8[1:5]) + scale_fill_manual(values=c("white","black"))+facet_grid(~year)
ggsave(file="Supp_Figure_1B.png", plot = last_plot(), width = 10, height = 10)


# Production of figure transcript functions by season [Figure S3, Figure 3]
ps_merged_subset <- ps_merged
f=ddply(ps_merged_subset, .(month, year, V5, V3, V4, new_V6, Sample), summarise, SUM = sum(norm))
f2=ddply(f, .(year, month, V5, Sample), summarise, SUM2 = sum(SUM))
f3=ddply(f2, .(year, month, V5), summarise, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))
temp = subset(f3, month == 'Sep')
f3$V5 <- factor(f3$V5, levels = temp$V5[order(-temp$MEAN)])
f3$month <- factor(f3$month, levels = c("May", "Jun", "Jul", "Aug", "Sep"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3$year <- str_replace(f3$year, "16", "Switchgrass 2016")
f3$year <- str_replace(f3$year, "17", "Switchgrass 2017")
p = ggplot(f3, aes_string(x="V5", y="MEAN", color="month"))
p+geom_point(stat="identity")+theme_bw()+geom_errorbar(limits, width=0)+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=12))+theme(strip.text.x = element_text(size = 10))+
  xlab("Subsystem Function") + ylab("Transcript Abundances")+ scale_color_manual(values=colorBlindBlack8[1:5]) + facet_grid(~year)
ggsave(file="Supp_Fig3.png", plot = last_plot(), width = 10, height = 10)

ps_merged_subset <- ps_merged
f=ddply(ps_merged_subset, .(month, V5, V3, V4, new_V6, Sample), summarise, SUM = sum(norm))
f2=ddply(f, .(month, V5, Sample), summarise, SUM2 = sum(SUM))
f3=ddply(f2, .(month, V5), summarise, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))
temp = subset(f3, month == 'Sep')
f3$V5 <- factor(f3$V5, levels = temp$V5[order(-temp$MEAN)])
f3$month <- factor(f3$month, levels = c("May", "Jun", "Jul", "Aug", "Sep"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3, aes_string(x="V5", y="MEAN", color="month"))
p+geom_point(stat="identity", size=3)+theme_bw()+geom_errorbar(limits, width=0)+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, size=12))+theme(strip.text.x = element_text(size = 10))+
  xlab("Subsystem Function") + ylab("Transcript Abundances")+ scale_color_manual(values=colorBlindBlack8[1:5]) 
ggsave(file="Fig_3.png", plot = last_plot(), width = 10, height = 10)


#stats
ps_merged_subset <- ps_merged
f=ddply(ps_merged_subset, .(date2, month, V5, V3, V4, new_V6, Sample), summarise, SUM = sum(norm))
f2=ddply(f, .(date2, V4, V5, V3), summarise, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
t1 <- f2 %>% select(date2,  V4, MEAN) %>% spread(date2, MEAN)
t1$ratio <- t1$`Aug - Sep`/t1$`May - July`
write.table(t1, file="ratio_of_subsystems.txt", sep="\t", quote=FALSE)

#Enriched and Core
ps_merged_subset <- ps_merged
f=ddply(ps_merged_subset, .(date2, V5, V3, V4, V6, Sample), summarise, SUM = sum(norm))
f2=ddply(f, .(date2, V4, Sample), summarise, SUM2 = sum(SUM))
f3=ddply(f2, .(date2, V4), summarise, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))
v1 <- f3 %>% select(date2,  V4, MEAN) %>% spread(date2, MEAN)
meta <- f %>% select(V4, V3, V5)
meta2 <- unique(meta)
t4 <- merge(v1, meta2, by = "V4")
t4$ratio <- t4$`Aug - Sep`/t4$`May - July`
t5 <- t4[t4$ratio != Inf,]
summary(t5$ratio) #MEAN is 35
t6 <- t4[t4$ratio == Inf | t4$ratio > 35,]
t6 <- t4[t4$ratio == Inf | t4$ratio > 1,]
enriched <- unique(t6$V4) #these are the V4 functions which are above average enriched between late and early

# Figure production Supp Figure 4 - Core MAGs
ps_merged_ubiq <- subset(ps_merged, Abundance > 0)
ps_merged_ubiq <- subset(ps_merged_ubiq, date2 == "Aug - Sep")
f=ddply(ps_merged_ubiq , .(V5, V4, V3, V6), summarise, MEAN = mean(norm))
f2 = ddply(f, .(V5, V4, V3), summarise, counts = length(V6))
f2b = ddply(f, .(V4), summarise, counts = length(V6))
summary(f2b$counts)
f3 <- subset(f2, counts >= 13) #present in more than 8 MAGS
#write.table(f3, file = "shared-transcripts.txt", quote=FALSE, sep="\t", row.names = FALSE)
f2$counts <- as.numeric(f2$counts)
#temp = subset(f3, date2 == "Aug - Sep")
f2$V5 <- factor(f2$V5, levels = temp$V5[order(-temp$MEAN)])
head(f2)
f2$flag <- 0
f2[f2$count >= 13,]$flag <- 1
p = ggplot(f2, aes_string(x="V5", y="counts", color="flag"))
p+theme_bw()+geom_jitter()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = "none")+
  xlab('Functional Subsystem') + ylab('Number of MAGs')
ggsave(file="Supp_Fig_4.png", plot = last_plot(), width = 10, height = 10)


flagged <- subset(f2, flag == 1)
v4_list <- unique(f3$V4) #these are the functions that are in > 8 MAGs
both <- intersect(v4_list, enriched)
dim(both)
t7 <- t6[t6$V4 %in% both,]
head(t7)
t8 <- merge(t7, f3, by = "V4")
t8$counts <- as.numeric(t8$counts)
t8$V3.y <- NULL
t8$V5.y <- NULL
colnames(t8) <- c("Function Role", "Aug - Sep Abundance", "May - Jul Abundance", "Functional Class", "Functional Group", "Ratio Late to Early Abundances", "Number of MAGs")
write.table(t8, file="Supp_Table_4-2.tsv", sep="\t", quote=FALSE, row.names = FALSE)
enriched_ps <- subset(ps_merged, V4 %in% t8$V4)
our_ann <- read.csv(file="enriched.txt", header=TRUE,sep="\t")
enriched_ps2 <- merge(enriched_ps, our_ann, by = "V4")

# Figure production for Supp Fig 5  - MAG clusters
library(vegan)
enrich_ps <- subset(ps_merged, date2 == "Aug - Sep")
enrich_ps <- ps_merged
f=ddply(enrich_ps, .(Sample, V6), summarise, SUM = sum(norm))
y <- f %>% spread(Sample, SUM)
rownames(y) <- y$V6
y$V6 <- NULL
dist.mat = vegdist(y, method="bray")
clusters <- hclust(dist.mat)
png('Supp_Figure_5.png')
plot(clusters)
dev.off()

# Figure production for Figure 4 - Abundance of transcripts on MAG clusters
grouping <- read.csv(file="grouping_2.txt", sep="\t", header=FALSE)
name_map <- read.csv(file="new-names.txt", sep="\t", header=FALSE)
colnames(name_map) <- c("new_V6", "old_V6")
colnames(grouping) <- c("old_V6", "Group")
grouping2 <- merge(grouping, name_map, by = "old_V6")
gpm <- merge(ps_merged, grouping2, by="new_V6")
f=ddply(gpm, .(month, Group, Sample, V5), summarise, SUM = sum(norm))
f2=ddply(f, .(month, Group, V5), summarise, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
temp = subset(f2, month == "Sep" & Group == "2" )
f2$V5 <- factor(f2$V5, levels = temp$V5[order(-temp$MEAN)])
f2$month <- factor(f2$month, levels = c("May", "Jun", "Jul", "Aug", "Sep"))
f2$Group <- as.character(f2$Group)
f2$zero_flag = 0
f2[f2$MEAN != 0,]$zero_flag= 1
f2$zero_flag <- as.character(f2$zero_flag)
p = ggplot(f2, aes_string(x="V5", y="MEAN"))
p+geom_point(stat="identity", size=2, shape=21, aes(fill=zero_flag))+theme_bw()+geom_errorbar(limits, width=0)+
  theme(text=element_text(size=12))+facet_grid(Group~month)+scale_shape(solid=FALSE)+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5, size=12))+scale_fill_manual(values=c("white","black"))+
  xlab("Subsystem Function") + ylab("Transcript Abundances")
ggsave(file="Fig_4.png", plot = last_plot(), width = 20, height = 10)

subset_aug <- subset(f, month == "Aug")
ss <- unique(subset_aug$V5)
ss_stats <- data.frame(matrix(NA, nrow = length(ss), ncol = 2))
rownames(ss_stats) <- ss
for (i in 1:length(ss)){
  ss_data <- subset(subset_aug, V5 == ss[i])
  res <- kruskal.test(SUM ~ Group, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}
subset_aug <- subset(f, month == "Sep")
ss <- unique(subset_aug$V5)
ss_stats <- data.frame(matrix(NA, nrow = length(ss), ncol = 2))
rownames(ss_stats) <- ss
for (i in 1:length(ss)){
  ss_data <- subset(subset_aug, V5 == ss[i])
  res <- kruskal.test(SUM ~ Group, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}


grouping <- read.csv(file="grouping_2.txt", sep="\t", header=FALSE)
name_map <- read.csv(file="new-names.txt", sep="\t", header=FALSE)
colnames(name_map) <- c("new_V6", "old_V6")
colnames(grouping) <- c("old_V6", "Group")
grouping2 <- merge(grouping, name_map, by = "old_V6")
gpm <- merge(ps_merged, grouping2, by="new_V6")
sulf <- subset(gpm,V5 == "Sulfur Metabolism")
f=ddply(sulf, .(month, Group, Sample, V3), summarise, SUM = sum(norm))
f2=ddply(f, .(month, Group, V3), summarise, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
temp = subset(f2, month == "Sep" & Group == "2" )
f2$V3 <- factor(f2$V3, levels = temp$V3[order(-temp$MEAN)])
f2$month <- factor(f2$month, levels = c("May", "Jun", "Jul", "Aug", "Sep"))
f2$Group <- as.character(f2$Group)
f2$zero_flag = 0
f2[f2$MEAN != 0,]$zero_flag= 1
f2$zero_flag <- as.character(f2$zero_flag)
f2 <- subset(f2, V3 != "NA")
p = ggplot(f2, aes_string(x="V3", y="MEAN"))
p+geom_point(stat="identity", size=2, shape=21, aes(fill=zero_flag))+theme_bw()+geom_errorbar(limits, width=0)+
  theme(text=element_text(size=12))+facet_grid(Group~month)+scale_shape(solid=FALSE)+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5, size=12))+scale_fill_manual(values=c("white","black"))+
  xlab("Subsystem Function") + ylab("Transcript Abundances")
ggsave(file="Supp_Fig_6.png", plot = last_plot(), width = 10, height = 10)


#stats on what is different beteween grouped MAGs
f=ddply(gpm, .(date2, Group, Sample, V5), summarise, SUM = sum(norm))
f2<-subset(f, date2 == "Aug - Sep")
all_V5 = unique(f2$V5)
for (i in 1:length(all_V5))
  {
  all_V5[i]
  subset_dat = subset(f2, V5 == all_V5[i])
  aov.result <- aov(SUM~as.factor(Group), data=subset_dat)
 
  #print(model.tables(aov.result,"means")) 
  #boxplot(SUM~group2, data=subset_dat) 
  if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05)
    {
    #print(all_V5[i])
    print(i)
    print(summary(aov.result))
    #print(TukeyHSD(aov.result))
    boxplot(SUM~Group, data=subset_dat) 
  }
}

# Figure production Supp Figure 7 Other Perennial Grasses - Other MetaGs/MetaTs

y <- read.csv(file="hits-against-iowa-misc2.txt", sep='\t')
dim(y)
y2 <- y[1:60641,1:7]
colnames(y2) = colnames(y)[2:length(colnames(y))]
y3 <- read.csv(file="others-metags2.tsv", sep='\t')
dim(y3)
y4 <- y3[1:60641,1:276]
colnames(y4) = colnames(y3)[2:length(colnames(y3))]
abundance <- otu_table(as.matrix(y2), taxa_are_rows=TRUE)
new_y <- merge(y2, y4, by="row.names")
rownames(new_y) <- new_y$Row.names
new_y$Row.names <- NULL
abundance <- otu_table(as.matrix(new_y), taxa_are_rows=TRUE)

meta <- read.csv("sample-list.txt", sep="\t", header=TRUE, row.names = 1)
phy_meta <- sample_data(meta)
annotations_mags = read.csv("mag.x.refsoil.blastpout.columns.l4.plus.finalann2", sep="\t", header=FALSE, row.names=1)
head(annotations_mags)
annotation_mags2 <- subset(annotations_mags, annotations_mags$V5 != "NA")
foo2<-as.matrix(annotation_mags2)
annotation <- tax_table(foo2)
phy_iowa <- phyloseq(abundance, annotation, phy_meta)
phy_iowa  <- prune_samples(sample_sums(phy_iowa)>0, phy_iowa)
phy_iowa  <- filter_taxa(phy_iowa , function(x) sum(x) > 0, TRUE)
phy_iowa_melt <- psmelt(phy_iowa)
phy_iowa_melt <- subset(phy_iowa_melt, Plant != "Corn")
phy_iowa_melt$Year <- as.character(phy_iowa_melt$Year)
phy_iowa_melt$new_V6 <- mapvalues(phy_iowa_melt $V6, from=og_names, to=new_names)
f=ddply(phy_iowa_melt, .(new_V6, Sample, Year.State, Plant), summarise, MEAN=mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f <- f[f$MEAN > 1,]
order_mags <- unique(f$new_V6)[order(nchar(unique(f$new_V6)), unique(f$new_V6))]
f$new_V6 <- factor(f$new_V6, levels = order_mags)
p = ggplot(f, aes(new_V6, y=MEAN, color=Year.State))+geom_point(stat="identity", size=2)+geom_errorbar(limits, width=0)+geom_jitter()
p = p +theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=0.5, size=10))+facet_grid(~Plant)
p + xlab("Phyllosphere MAG") + ylab("Average Abundance (Reads Mapped)")
ggsave(file="Supp_Fig_7.png", plot = last_plot(), width = 10, height = 10)


