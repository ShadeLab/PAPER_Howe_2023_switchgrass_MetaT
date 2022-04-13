# Packages needed for analysis
library(tidyr)
library(dplyr)
library(plyr)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(corrplot)
library(Hmisc)
library(stringr)
library(DESeq2)

# Set local workign directory where files are located
setwd("~/Downloads/phyllosphere-analysis/")
#Sys.setenv('R_MAX_VSIZE'=6000000000)

# Input Files needed
# sample-meta.https://www.facebook.com/txt = metadata describing samples
# checkm/checkm_output.txt = output from checkM describing MAGs
# melted_phy_average_coverage_mags_per_metag.rds = R data object, summary of mapped reads to all potential MAGs

# MetaG Files
# contigs-unnormalized-median-metag-no-fungal.txt = quality filtered median bp coverage of metag contigs to MAGs
# list-of-hkg.txt = list of contigs associated with housekeeping genes
# 2016-17_MetaG_map.txt = MetaG metadata of samples

# MetaT Files
# orfs-unnormalized-median-metat-no-fungal.txt = quality filtered median bp coverage of metat ORFs to MAGs
# draft_clustering.txt = assignment of MAGs to clusters based on pvclust
# list-of-hkg2.txt = list of ORFs associated with housekeeping genes
# orfs-normalized-median-metat-no-fungal-no-ann-uniquified.txt = quality filetered median bp coverage of metaT ORFs to MAGs (but with unique ID to each unique KEGG annotation)
# 2016-17_MetaT_map.txt = MetaT metadata of samples
# DRAM2.txt = annotation of ORFs using DRAM and brings in KEGG 
# mags-in-sg-and-misc2.csv = detection of mags in sg and misc
# orf-map.txt = mapping the ORFs to their originating MAG

# Other MetaGs for Comparisons
# other-metags-summary.txt = median bp coverage of MAGs in publicly availble metaGs
# sample-list.txt = metadata for these samples 


###################################################################################################
# Contamination of plant host and fungal genomes - Figure S1
##############################################################################################


meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
contam_metag = read.delim(sep="\t", file="./MetaG_HostReport2.tsv", header=TRUE, strip.white=TRUE, row.names = 1)
cm_merged = merge(contam_metag, meta, by = "row.names")
rownames(cm_merged) <- cm_merged$Row.names
cm_merged$pairs_aligned <- (cm_merged[,5]+cm_merged[,6])/(cm_merged[,3])
cm_merged$year <- as.factor(cm_merged$year)
f <- ddply(cm_merged, .(plant, year, month), summarise, MEAN=mean(pairs_aligned), SD=sd(pairs_aligned))
limits<-aes(ymin=MEAN-SD, ymax=MEAN+SD)
f$month <- factor(month.name[f$month], levels=month.name[1:12])
p = ggplot(f, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Reads Aligned to Host Plants (%)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(~plant)

meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
contam_metag = read.delim(sep="\t", file="./MetaG_FungalAln2.tsv", header=TRUE, strip.white=TRUE, row.names = 1)
cm_merged = merge(contam_metag, meta, by = "row.names")
rownames(cm_merged) <- cm_merged$Row.names
cm_merged$pairs_aligned <- (cm_merged[,5]+cm_merged[,6])/(cm_merged[,3])
cm_merged$year <- as.factor(cm_merged$year)
f <- ddply(cm_merged, .(plant, year, month), summarise, MEAN=mean(pairs_aligned), SD=sd(pairs_aligned))
f$month <- factor(month.name[f$month], levels=month.name[1:12])
limits<-aes(ymin=MEAN-SD, ymax=MEAN+SD)
p = ggplot(f, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Reads Aligned to Fungal Genomes (%)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(~plant)

meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
contam_metag = read.delim(sep="\t", file="./MetaT_HostAln2.tsv", header=TRUE, strip.white=TRUE, row.names = 1)
cm_merged = merge(contam_metag, meta, by = "row.names")
rownames(cm_merged) <- cm_merged$Row.names
cm_merged$pairs_aligned <- (cm_merged[,5]+cm_merged[,6])/(cm_merged[,3])
cm_merged$year <- as.factor(cm_merged$year)
f <- ddply(cm_merged, .(plant, year, month), summarise, MEAN=mean(pairs_aligned), SD=sd(pairs_aligned))
limits<-aes(ymin=MEAN-SD, ymax=MEAN+SD)
f$month <- factor(month.name[f$month], levels=month.name[1:12])
p = ggplot(f, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Reads Aligned to Host Plants (%)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
contam_metag = read.delim(sep="\t", file="./MetaT_FungalAln2.tsv", header=TRUE, strip.white=TRUE, row.names = 1)
cm_merged = merge(contam_metag, meta, by = "row.names")
rownames(cm_merged) <- cm_merged$Row.names
cm_merged$pairs_aligned <- (cm_merged[,5]+cm_merged[,6])/(cm_merged[,3])
cm_merged$year <- as.factor(cm_merged$year)
f <- ddply(cm_merged, .(plant, year, month), summarise, MEAN=mean(pairs_aligned), SD=sd(pairs_aligned))
f$month <- factor(month.name[f$month], levels=month.name[1:12])
limits<-aes(ymin=MEAN-SD, ymax=MEAN+SD)
p = ggplot(f, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Reads Aligned to Fungal Genomes (%)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###################################################################################################
# Characteristics of MAGs and selection of FOCAL MAGs
# Figure 3
###################################################################################################

meta <- read.delim(sep="\t", file="./sample-meta.txt", header=FALSE, strip.white=TRUE, row.names = 1)
colnames(meta) <- c("plot", "fert", "date", "id", "plant", "date2", "month")
meta$year <- format(as.Date(meta$date, format="%d-%b-%Y"), "%Y")
# CheckM Annotations with Selected Columns and filter by completeness and contamination - see https://github.com/ShadeLab/PAPER_Howe_2021_switchgrass_MetaT/tree/main/mag-evaluation
checkm <- read.delim(sep='\t', file="./checkm (adina@iastate.edu)/checkm_output.txt", header = TRUE, strip.white=TRUE)
checkm_qc <- subset(checkm, completeness >= 50 & contamination < 10)
checkm_qc <- checkm_qc[!checkm_qc$plantbin %in% c("M64"),]
list_of_interest <- checkm_qc$plantbin
list_of_interest <- list_of_interest[-25]
rownames(checkm) <- checkm$plantbin

# Read in Average Estimated Median Read Coverage per MAG (i.e., plantbin)
#sg <- rownames(subset(meta, plant == "S"))
cov <- readRDS('melted_phy_average_coverage_mags_per_metag.rds')
#cov2 <- cov[cov$Sample %in% sg,]
cov2 <- cov %>% select(Sample, plantbin, MEAN) %>% spread(Sample, MEAN)
#cov2 <- cov2 %>% select(Sample, plantbin, MEAN) %>% spread(Sample, MEAN)
rownames(cov2) <- cov2$plantbin
cov2$plantbin <- NULL
cov2.abund <- apply(cov2, 1, mean)
cov2.freq <- rowSums(cov2 != 0)
occ <- merge(cov2.abund, cov2.freq, by = 0)
rownames(occ) <- occ$Row.names
occ$Row.names <- NULL
colnames(occ) <- c("abundance", "prevalence")
occ_sorted <- occ[order(occ$abundance),]
ggplot(data = occ_sorted, aes(x=log(abundance), y = prevalence)) + geom_point() 
occ_sorted$inclusion = "not focal"
occ_sorted[rownames(occ_sorted) %in% list_of_interest,]$inclusion <- "focal"
#occ_sorted$misc <- 0
#occ_sorted$sg <- 0
#occ_sorted$miscandsg <- 0
#occ_sorted[rownames(sg_list),]$sg <- 1
#occ_sorted[rownames(misc_list),]$misc <- 1
#occ_sorted$miscandsg <- occ_sorted$sg + occ_sorted$misc
occ_sorted2 <- merge(occ_sorted, checkm, by="row.names")
# Some checkm contamination > 100, adjusted to 100.
occ_sorted2[occ_sorted2$contamination > 100,]$contamination = 100
ggplot(data = occ_sorted2, aes(x=log(abundance), y = prevalence, shape=inclusion, color = contamination, size = completeness)) + geom_point() +
  scale_shape_manual(values = c(19, 43))+scale_color_gradientn(colours = rev(rainbow(5)))+ xlab('Log Abundance (average median bp MAG coverage)') +
  ylab('Occupancy (n=136 metagenomes)')
ggsave('focalmag.png')
ggsave('focalmag.eps', dev="eps")
#136 metagenomes - 2016 switchgrass and miscanthus

###################################################################################################
# Characteristics of MAGs in Metagenomes - Clusters in Metagenomes by Month
# Figure 4
###################################################################################################

cov <- read.delim(file="contigs-unnormalized-median-metag-no-fungal.txt", sep="\t", row.names = 1)
dim(cov)
new_col_names = colnames(cov)[2:193]
cov[,193] <- NULL
dim(cov) #192 metags
new_col_names = str_split_fixed(new_col_names, "_L", 2)[,1]
new_col_names = str_split_fixed(new_col_names, '.sam', 2)[,1]
colnames(cov) <- new_col_names
meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
abundance <- otu_table(as.matrix(cov), taxa_are_rows=TRUE)
meta2 <- meta %>% mutate(month_range = ifelse(month < 7, as.character(1), as.character(2)))
metadata <- sample_data(meta2)
phy <- phyloseq(metadata, abundance)
phy2 = prune_taxa(taxa_sums(phy) > 0, phy)
phy3 = prune_samples(sample_sums(phy2) > 0, phy2) #139 samples
# Housekeeping genes trasnlated ORFs
x <- scan("./list-of-hkg-contigs.txt", what="", sep="\n")
hkg_phy <- prune_taxa(x, phy3)
phy_melt_hkg = psmelt(hkg_phy)
f <- ddply(phy_melt_hkg, .(Sample, plant, year, month), summarise, SUM=sum(Abundance))
f_hkg_metag <- f
f_hkg <- f %>% select(Sample, SUM)
zero_samples = subset(f, SUM==0)$Sample #25 samples
#f2 <- ddply(f, .(plant, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
#limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
#p = ggplot(f2, aes_string(x="month", y="MEAN", color="month"))
#p+geom_point(stat="identity") +theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=16), axis.text.x = element_text(angle = 90, hjust=0, size=10))
abundance <- as.data.frame(t(otu_table(phy3)))
abundance$Sample <- rownames(abundance)
hkg_foo <- merge(abundance, f_hkg, by="Sample")
rownames(hkg_foo) <- hkg_foo$Sample
hkg_foo$Sample <- NULL
hkg_norm <- hkg_foo/hkg_foo$SUM
hkg_norm$SUM <- NULL
abundance <- otu_table(as.matrix(t(hkg_norm)), taxa_are_rows=TRUE)
metadata2 = subset(metadata, !(rownames(metadata) %in% zero_samples))
metadata <- sample_data(metadata2)
phy <- phyloseq(metadata, abundance)  #this is the hkg normalized phyloseq object
phy_metag <- phy
#save(phy_metag, file="metag_physeq.RData")
#foo <- psmelt(phy_metag)
#save(foo, file="metag_physeq_melted.RData")

load('metag_physeq_melted.RData')
all <-  foo %>% select(Sample, OTU, Abundance, plant, month, year)
all$MAG_bin <- str_split_fixed(all$OTU, "_", 3)[,3]
f <- ddply(all, .(Sample, MAG_bin), summarise, MEAN=mean(Abundance))
f2 <- f %>% spread(MAG_bin, MEAN)
meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta$Sample <- rownames(meta)
f3 <- merge(f2, meta, by = "Sample")
#write.table(f3, file="metag_MAG_mean_median_bp_coverage.tsv", sep="\t", quote=FALSE)

f <- ddply(all, .(Sample, MAG_bin, month, plant, year), summarise, MEAN=mean(Abundance))
t<-str_split_fixed(f$MAG_bin, '_bin.',2)
new_bins <- paste(t[,1],t[,2], sep="")
f$bin <- new_bins
map_cluster = read.delim(file="draft_clustering.txt", sep="\t")
colnames(map_cluster) = c("bin", "phylum", "cluster")
f2 <- merge(f,map_cluster, by="bin" )
f2 <- ddply(f2, .(month, cluster, plant, year), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN)))
temp = subset(f2, month == "9")
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
f2$month <- factor(f2$month, levels=c("5", "6", "7", "8", "9", "10", "11"))
limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
f2 <- subset(f2, cluster != "-")
p = ggplot(f2, aes_string(x="month", y="MEAN2", color="cluster", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(~plant)
# Figure 3

###################################################################################################
# Characteristics of MAGs in Metatranscriptomes 
# Figure 4
###################################################################################################

cov <- read.delim(file="orfs-unnormalized-median-metat-no-fungal.txt", sep="\t", row.names = 1)
dim(cov)
new_col_names = colnames(cov)[2:79]
cov[,79] <- NULL
dim(cov)
new_col_names = str_split_fixed(new_col_names, "_L", 2)[,1]
colnames(cov) <- new_col_names
meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
abundance <- otu_table(as.matrix(cov), taxa_are_rows=TRUE)
meta2 <- meta %>% mutate(month_range = ifelse(month < 7, as.character(1), as.character(2)))
metadata <- sample_data(meta2)
#kegg_ann <- read.delim(sep="\t", file="./orf-kegg-annotations.txt", row.names = 1, header=FALSE)
#ann <- tax_table(as.matrix(kegg_ann))
phy <- phyloseq(metadata, abundance)
phy2 = prune_taxa(taxa_sums(phy) > 0, phy)
phy3 = prune_samples(sample_sums(phy2) > 0, phy2)
# Housekeeping genes trasnlated ORFs
x <- scan("./list-of-hkg2.txt", what="", sep="\n")
# Separate elements by one or more whitepace
y <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
hkg_orfs <- sapply(y, `[[`, 1)
hkg_phy <- prune_taxa(hkg_orfs, phy3)
phy_melt_hkg = psmelt(hkg_phy)
f <- ddply(phy_melt_hkg, .(Sample, plant,year, month), summarise, SUM=sum(Abundance))
f_hkg_metat <- f
f_hkg <- f %>% select(Sample, SUM)
#f2 <- ddply(f, .(plant, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
#limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
#p = ggplot(f2, aes_string(x="month", y="MEAN", color="month"))
#p+geom_point(stat="identity") +theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=16), axis.text.x = element_text(angle = 90, hjust=0, size=10))
abundance <- as.data.frame(t(otu_table(phy3)))
abundance$Sample <- rownames(abundance)
hkg_foo <- merge(abundance, f_hkg, by="Sample")
rownames(hkg_foo) <- hkg_foo$Sample
hkg_foo$Sample <- NULL
hkg_norm <- hkg_foo/hkg_foo$SUM
hkg_norm$SUM <- NULL
abundance <- otu_table(as.matrix(t(hkg_norm)), taxa_are_rows=TRUE)
metadata2 = subset(metadata, rownames(metadata) != "G5R1_NF_31MAY2016")
metadata <- sample_data(metadata2)
phy <- phyloseq(metadata, abundance)  #this is the hkg normalized phyloseq object
phy2 = prune_taxa(taxa_sums(phy) > 0, phy)
phy3 = prune_samples(sample_sums(phy2) > 0, phy2)
#ps <- psmelt(phy3)
summary(taxa_sums(phy3))
#first quartile of taxa sums and 10% of samples
phy4 = filter_taxa(phy3, function(x) sum(x > 0.000135) > (0.1*length(x)), TRUE)
ps = psmelt(phy4) #This is the melted filtered phy object used to do analyses
#save(ps, file="melted_phy3.RData") #no annotations and normalized
#save(ps, file="phy3_line192.RData") #no annotations and normalized
annotation = read.delim(file="../MAG_annotation/DRAM2.txt", header=FALSE)
rownames(annotation) = annotation$V1
annotation$V1 <- NULL
annotation_phy = tax_table(as.matrix(annotation))
phy_ann = phyloseq(metadata, abundance, annotation_phy)
#write.table(otu_table(phy3), file="orfs-normalized-median-metat-no-fungal-no-ann.txt", sep="\t", quote=FALSE)
load(file="melted_phy3.RData") #loads object ps
phy3


###################################################################################################
#MAG Seasonal Patterns in MetatT
###################################################################################################
load(file="melted_phy3.RData") #loads object ps, filtered
head(ps)
split1 <- str_split_fixed(ps$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
bin = paste(split2a, split3, sep='_')
ps$bin <- bin
all <- ps %>% select(Sample, year, OTU, Abundance, bin, month, year)
f <- ddply(all, .(Sample, year, bin, month), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(bin,  month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9")
f2$month <- as.character(f2$month)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Sampling Month') + ylab('Cumulative Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(.~bin, ncol = 6)

f <- ddply(all, .(Sample, year, bin, month), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(month, year), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9" & year == "2017")
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('Month') + ylab('Summed Coverage of 41 MAGs Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###################################################################################################
#Functional Analysis of MetaT ORFs against KEGG Subsystems
# Figure 5
###################################################################################################

cov <- read.delim(file="orfs-normalized-median-metat-no-fungal-no-ann-uniquified.txt", sep="\t", row.names = 1)
abundance <- otu_table(as.matrix(cov), taxa_are_rows=TRUE)
meta2 <- meta %>% mutate(month_range = ifelse(month < 7, as.character(1), as.character(2)))
metadata <- sample_data(meta2)
ann = read.delim(file="../MAG_annotation/DRAM2.txt", header=FALSE)
rownames(ann) <- ann$V1
ann$V1 <- NULL
ann2 = tax_table(as.matrix(ann))
phy = phyloseq(abundance, metadata, ann2)
#this phy object is ready for functional parsing above - !!!
phy2 = prune_taxa(taxa_sums(phy) > 0, phy)
phy3 = prune_samples(sample_sums(phy2) > 0, phy2)
summary(taxa_sums(phy3))
#first quartile of taxa sums and 10% of samples
phy4 = filter_taxa(phy3, function(x) sum(x > 0.000145) > (0.1*length(x)), TRUE)
ps <- psmelt(phy4)
#save(ps, file="melted_phy4-line224.RData") #no annotations and normalized
load("melted_phy4-line224.RData") #as ps

# # confirming bin M22 is absent
# split1 <- str_split_fixed(ps$OTU, "_", 3)[,3]
# split2a <- str_split_fixed(split1, '_', 2)[,1]
# split2b <- str_split_fixed(split1, '_', 2)[,2]
# split3 <- str_split_fixed(split2b, '_', 2)[,1]
# bin = paste(split2a, split3, sep='_')
# ps$bin <- bin
# bin = paste(split2a, split3, sep='_')
# subset(ps, bin == "M_bin.22")

all <- ps %>% select(Sample, year, OTU, Abundance, day, month, month_range, year, V7)
f <- ddply(all, .(Sample, year, month, month_range, V7), summarise, SUM=sum(Abundance))
fb <- subset(f, V7 != 'NA')
fb <- subset(f, V7 != 'Global and overview maps')
f2 <- ddply(fb, .(V7,  year, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9" & year == "2017")
f2$V7 <- factor(f2$V7, levels=temp$V7[order(-temp$MEAN)])
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="V7", y="MEAN", color="month", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=12), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('KEGG Metabolism Classification') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(~month)

p = ggplot(f2, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=12), axis.text.x = element_text(angle = 0, hjust=1, size=15))+
  xlab('KEGG Metabolism Classification') + ylab('Coverage Normalized by HKGs')+scale_shape_manual(values=c(16,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~V7, ncol = 7, labeller = labeller(V7 = label_wrap_gen(10)), scales="free_y")
#Figure 5                                                                              


library("ggpubr")
ggplot(f, aes(x=SUM)) + geom_density() #not normal
#kruskal-wallis test 
ss <- unique(f$V7[1:22])
ss_stats <- data.frame(matrix(NA, nrow = length(ss), ncol = 2))
rownames(ss_stats) <- ss
for (i in 1:length(ss)){
  ss_data <- subset(f, V7 == ss[i])
  res <- kruskal.test(SUM ~ month_range, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}

#> ss_stats[ss_stats$X2 > 0.01,]
#X1         X2
#Cell motility                    0.002730784 0.95832399
#Cellular community - prokaryotes 2.828408398 0.09261010
#Environmental adaptation         2.083931480 0.14885635
#Membrane transport               4.343205928 0.03715690
#Signal transduction              4.213565319 0.04010197
#Transcription                    1.069450666 0.30106956
#Translation                      4.630223677 0.03141347

###########

ps2 <- ps %>% select(Sample, OTU, Abundance)
split1 <- str_split_fixed(ps2$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
bin = paste(split2a, split3, sep='_')
ps2$bin <- bin
f <- ddply(ps2, .(Sample, bin), summarise, SUM=sum(Abundance))
f2 <- f %>% spread(bin, SUM)
meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta$Sample <- rownames(meta)
f3 <- merge(f2, meta, by = "Sample")
#write.table(f3, file="metat_MAG_sum_median_bp_coverage.tsv", sep="\t", quote=FALSE)


foo <- ps2 %>% group_by(OTU, month_range) %>% summarise_at(vars(Abundance), list(MEAN=mean))
foo2 <- foo %>% spread(month_range, MEAN)
kegg_ann$OTU <- rownames(kegg_ann)
foo3 <- merge(foo2, kegg_ann, by="OTU")
head(foo3)
colnames(foo3) = c("OTU", "early", "late", "evalue", "Kegg1", "Kegg2", "Kegg3", "Kegg4")
foo3$ratio <- foo3$late / foo3$early
#write.table(foo3, file = 'persistent_functions.txt', sep='\t', quote=FALSE)
genes_present_both = subset(foo3, ratio != "Inf")
genes_only_late = subset(foo3, ratio == "Inf")
genes_present_both$MEAN = (genes_present_both$early + genes_present_both$late)/2
genes_present_both$MEAN_log = log(genes_present_both$MEAN)
genes_present_both$ratio_log = log(genes_present_both$ratio)
genes_present_both2 = subset(genes_present_both, Kegg1 != "NA")
p = ggplot(genes_present_both2, aes_string(x="Kegg1", y="ratio", color="Kegg1"))
p+geom_point(stat="identity") +theme_bw()+theme(text=element_text(size=16), axis.text.x = element_text(angle = 90, hjust=0, size=10))

all <- ps %>% select(Sample, OTU, Abundance, day, month, year)
all2 <- merge(all, kegg_ann, by = "OTU")
f <- ddply(all, .(Sample, OTU,  month), summarise, SUM=sum(Abundance))
fb <- subset(f, V4 != 'NA')
fb <- subset(fb, V4 != "Alpha amylase, catalytic domain [PF00128.25]; Domain of unknown function (DUF3459) [PF11941.9]; Maltogenic Amylase, C-terminal domain [PF16657.6]")
f2 <- ddply(fb, .(V4,  month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9")
f2$V4 <- factor(f2$V4, levels=temp$V4[order(-temp$MEAN)])
f2$month <- as.character(f2$month)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="V4", y="MEAN", color="month"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
    xlab('KEGG Metabolism Classification') + ylab('Coverage Normalized by HKGs')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

orf_map <- read.delim(file = 'orf-map.txt', sep=' ', header=FALSE)
colnames(orf_map) = c("OTU", "bin")
foo <- merge(all, orf_map, by = "OTU")
f <- ddply(foo, .(Sample, bin, OTU), summarise, SUM=sum(Abundance))
fb <- subset(f, V4 != 'NA')
fb <- subset(fb, V4 != "Alpha amylase, catalytic domain [PF00128.25]; Domain of unknown function (DUF3459) [PF11941.9]; Maltogenic Amylase, C-terminal domain [PF16657.6]")
f2 <- ddply(fb, .(V4,  bin), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="V4", y="MEAN", color="bin"))
p+geom_point(stat="identity", size=3) + facet_grid(~V4)


# Characteristics of MAGs in Metatranscriptomes - Functions
cov =  read.delim(file="orfs-normalized-median-metat-no-fungal.txt", sep="\t", row.names = 1)
orf_list = read.delim(file="../MAG_annotation/unique-orfs-for-annotating.txt", header=FALSE)
orf_cov = merge(orf_list, cov, by.x = V1, by.y = "row.names")                      

###################################################################################################
# Clustering of Focal MAGs (filtered by Abundance of ORFs) - determination of clusters
# Figure S2
###################################################################################################

library(gplots)  
library(Heatplus)
library(pvclust)
library(RColorBrewer)
#phy #from line 57
#phy2 = prune_taxa(taxa_sums(phy) > 0, phy)
#phy3 = prune_samples(sample_sums(phy2) > 0, phy2)
#phy4 = filter_taxa(phy3, function(x) sum(x > 0.000142) > (0.1*length(x)), TRUE)
#ps = psmelt(phy4)
#save(ps, file="melted_phy4.RData")
load(file="melted_phy4.RData") #loads object ps
orf_map <- read.delim(file = 'orf-map.txt', sep=' ', header=FALSE)
colnames(orf_map) = c("OTU", "bin")
ps2 <- merge(ps, orf_map, by = "OTU")
#save(ps2, file="melted_phy4b.RData")
load(file="melted_phy4b.RData") #loads object ps2
f <- ddply(ps2, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
dim(f2)
f3 <- f2[,2:41] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
label_row <- meta[rownames(f3_log),]$month
result <- pvclust(f3, method.dist="euclidian", method.hclust="ward.D", nboot=10000, parallel=TRUE)
plot(result) #Clustering Dendrogram
pvrect(result, alpha=0.3)
clsig <- unlist(pvpick(result, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)

###################################################################################################
# Enrichment of MAGs by Season - Additional Heatmap figures (incuding clustering)
# 
###################################################################################################

#MetaT
# Figure 4
load(file="melted_phy4b.RData") #loads object ps2

f <- ddply(ps2, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
dim(f2)
f3 <- f2[,2:41] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
label_row <- meta[rownames(f3_log),]$month
result <- pvclust(f3, method.dist="euclidian", method.hclust="ward.D", nboot=10000, parallel=TRUE)
plot(result) #Clustering Dendrogram
pvrect(result, alpha=0.3)
clsig <- unlist(pvpick(result, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)
ordered_by_month <- rownames(meta[order(meta$month, meta$year),])
ordered_f3_log <- f3_log[ordered_by_month,]
ordered_f3_log2 <- ordered_f3_log[-5,]
#heatmap.2(as.matrix(ordered_f3_log2), col=brewer.pal(9,"Blues"), trace='none', Rowv=FALSE, labRow = meta[order(meta$month, meta$year),]$month, Colv=as.dendrogram(result$hclust), cexRow =1, distfun = function(x) dist(x,method = 'euclidean'))
heatmap.2(as.matrix(t(ordered_f3_log2)), col=brewer.pal(9,"Blues"), trace='none', Rowv=as.dendrogram(result$hclust), labCol = meta[order(meta$month, meta$year),]$month, Colv=FALSE, cexRow =1, distfun = function(x) dist(x,method = 'euclidean'))
order_of_dendrogram_labels = result$hclust$labels[order.dendrogram(as.dendrogram(result$hclust))]

ps2_subset = subset(ps2, year == "2016")
f <- ddply(ps2_subset, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
dim(f2)
f3 <- f2[,2:41] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta_subset <- subset(meta, year == "2016")
label_row <- meta_subset[rownames(f3_log),]$month
ordered_by_month <- rownames(meta_subset[order(meta_subset$month, meta_subset$year),])
ordered_f3_log <- f3_log[ordered_by_month,]
ordered_f3_log2 <- ordered_f3_log[-5,]
summary(melt(f3_log)$value)
colors <- c(0, .0027, 0.04588, 0.264, 4)
my_palette <- brewer.pal(4,"Blues")
heatmap.2(as.matrix(t(ordered_f3_log2)), col=my_palette, breaks=colors, trace='none', Colv = FALSE, Rowv=as.dendrogram(result$hclust), labCol = meta_subset[order(meta_subset$month, meta_subset$year),]$month,  cexRow =1, distfun = function(x) dist(x,method = 'euclidean'))

ps2_subset = subset(ps2, year == "2017")
f <- ddply(ps2_subset, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
dim(f2)
f3 <- f2[,2:41] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaT_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta_subset <- subset(meta, year == "2017")
label_row <- meta_subset[rownames(f3_log),]$month
ordered_by_month <- rownames(meta_subset[order(meta_subset$month, meta_subset$year),])
ordered_f3_log <- f3_log[ordered_by_month,]
ordered_f3_log2 <- ordered_f3_log[-5,]
summary(melt(f3_log)$value)
colors <- c(0, .0027, 0.04588, 0.264, 4)
my_palette <- brewer.pal(4,"Blues")
heatmap.2(as.matrix(t(ordered_f3_log2)), col=my_palette, breaks=colors, trace='none', Colv = FALSE, Rowv=as.dendrogram(result$hclust), labCol = meta_subset[order(meta_subset$month, meta_subset$year),]$month,  cexRow =1, distfun = function(x) dist(x,method = 'euclidean'))

#MetaG
load('metag_physeq_melted.RData') #loads as foo
all <-  foo %>% select(Sample, OTU, Abundance, plant, month, year)
split1 <- str_split_fixed(all$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, 'bin.', 2)[,2]
bin = paste(split2a, split3, sep='')
all$bin <- bin
sg <- subset(all, plant == "switchgrass" & year == "2016")
f <- ddply(sg, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
f3 <- f2[,2:42] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta_sg <- subset(meta, plant == "switchgrass" & year == "2016")
meta_sg2 <- meta_sg[rownames(meta_sg) %in% unique(sg$Sample),]
meta_sg <- meta_sg2
label_row <- meta_sg[rownames(f3_log),]$month
ordered_by_month <- rownames(meta_sg[order(meta_sg$month, meta_sg$year),])
ordered_f3_log <- f3_log[ordered_by_month,]
ordered_f3_log2  <- ordered_f3_log[,order_of_dendrogram_labels]
summary(melt(f3_log)$value)
colors <- c(0, .0027, 0.04588, 0.264, 4)
my_palette <- brewer.pal(4,"Blues")
matrix_foo = t(ordered_f3_log2)[nrow(t(ordered_f3_log2)):1,]
heatmap.2(as.matrix(matrix_foo), col=brewer.pal(9,"Blues"), trace='none', Colv=FALSE, Rowv= as.dendrogram(result$hclust), cexRow =1, distfun = function(x) dist(x,method = 'euclidean'), labCol = meta_sg[order(meta_sg$month, meta_sg$year),]$month)


sg <- subset(all, plant == "switchgrass" & year == "2017")
f <- ddply(sg, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
f3 <- f2[,2:42] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta_sg <- subset(meta, plant == "switchgrass" & year == "2017")
meta_sg2 <- meta_sg[rownames(meta_sg) %in% unique(sg$Sample),]
meta_sg <- meta_sg2
label_row <- meta_sg[rownames(f3_log),]$month
ordered_by_month <- rownames(meta_sg[order(meta_sg$month, meta_sg$year),])
ordered_f3_log <- f3_log[ordered_by_month,]
ordered_f3_log2  <- ordered_f3_log[,order_of_dendrogram_labels]
summary(melt(f3_log)$value)
colors <- c(0, .0027, 0.04588, 0.264, 4)
my_palette <- brewer.pal(4,"Blues")
matrix_foo = t(ordered_f3_log2)[nrow(t(ordered_f3_log2)):1,]
heatmap.2(as.matrix(matrix_foo), col=brewer.pal(9,"Blues"), trace='none', Colv=FALSE, Rowv= as.dendrogram(result$hclust), cexRow =1, distfun = function(x) dist(x,method = 'euclidean'), labCol = meta_sg[order(meta_sg$month, meta_sg$year),]$month)



sg <- subset(all, plant == "miscanthus" & year == "2016")
f <- ddply(sg, .(bin, Sample), summarise, SUM=sum(Abundance))
f2 = f %>% spread(bin, SUM)
f3 <- f2[,2:42] #M22 missing too low abundant
f3_log = log(f3+1)
rownames(f3_log) <- f2[,1]
meta <- read.delim(sep="\t", file="./2016-17_MetaG_map.txt", header=TRUE, strip.white=TRUE, row.names = 1)
meta_sg <- subset(meta, plant == "miscanthus" & year == "2016")
meta_sg2 <- meta_sg[rownames(meta_sg) %in% unique(sg$Sample),]
meta_sg <- meta_sg2
label_row <- meta_sg[rownames(f3_log),]$month
ordered_by_month <- rownames(meta_sg[order(meta_sg$month, meta_sg$year),])
ordered_f3_log <- f3_log[ordered_by_month,]
ordered_f3_log2  <- ordered_f3_log[,order_of_dendrogram_labels]
summary(melt(f3_log)$value)
colors <- c(0, .0027, 0.04588, 0.264, 4)
my_palette <- brewer.pal(4,"Blues")
matrix_foo = t(ordered_f3_log2)[nrow(t(ordered_f3_log2)):1,]
heatmap.2(as.matrix(matrix_foo), col=brewer.pal(9,"Blues"), trace='none', Colv=FALSE, Rowv= as.dendrogram(result$hclust), cexRow =1, distfun = function(x) dist(x,method = 'euclidean'), labCol = meta_sg[order(meta_sg$month, meta_sg$year),]$month)


###################################################################################################
#Cluster based analysis of MetaT and MetaG
###################################################################################################

head(ps2)
map_cluster = read.delim(file="draft_clustering.txt", sep="\t")
colnames(map_cluster) = c("bin", "phylum", "cluster")
map_cluster
ps3 <- merge(ps2, map_cluster, by = "bin")


f <- ddply(ps3, .(Sample, OTU,  month, year, V5), summarise, SUM=sum(Abundance))
fb <- subset(f, V5 != 'NA')
fb <- subset(fb, V4 != "Alpha amylase, catalytic domain [PF00128.25]; Domain of unknown function (DUF3459) [PF11941.9]; Maltogenic Amylase, C-terminal domain [PF16657.6]")
f2 <- ddply(fb, .(V5, year, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
temp = subset(f2, month == "9" & year == "2016")
f2$V5 <- factor(f2$V5, levels=temp$V5[order(-temp$MEAN)][1:20])
f2 <- subset(f2, V5 != "NA")
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="V5", y="MEAN", color="month", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme(text=element_text(size=15), axis.text.x = element_text(angle = 90, hjust=1, size=15))+
  xlab('KEGG Metabolism Classification') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(~year)

f <- ddply(ps3, .(Sample, cluster, month), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(cluster, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$month <- as.character(f2$month)
temp = subset(f2, month == "9")
f2$month <- as.character(f2$month)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", color="cluster"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

f <- ddply(ps3, .(Sample, cluster, year, month), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(cluster, year, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$month <- as.character(f2$month)
temp = subset(f2, month == "9")
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", color="cluster", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


f <- ddply(ps3, .(Sample, cluster, V5, year, month), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(cluster, V5, year, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$V5 <- factor(f2$V5, levels=temp$V5[order(-temp$MEAN)][1:20])
f2 <- subset(f2, V5 != "NA")
f2$month <- as.character(f2$month)
f2 <- subset(f2, V5 != "NA")
temp = subset(f2, month == "9" & year == "2016")
f2$zero_flag = "0"
f2[f2$MEAN != 0,]$zero_flag = "1"
f2$month <- as.character(f2$month)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", color="cluster"))
p+geom_point(stat="identity", size=3,  shape=21, aes(fill=zero_flag)) + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(cluster~V4)+ scale_fill_manual(values=c("white","black"))+
  theme(strip.text.x = element_text(angle = 90))+facet_grid(year~V5)

f <- ddply(ps3, .(Sample, cluster, V6, year, month), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(cluster, V6, year, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- subset(f2, V6 != "NA")
f2$month <- as.character(f2$month)
f2 <- subset(f2, V6 != "NA")
temp = subset(f2, month == "9" & year == "2016")
f2$V6 <- factor(f2$V6, levels=temp$V6[order(-temp$MEAN)][1:20])
f2$zero_flag = "0"
f2[f2$MEAN != 0,]$zero_flag = "1"
f2$month <- as.character(f2$month)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", color="cluster"))
p+geom_point(stat="identity", size=3,  shape=21, aes(fill=zero_flag)) + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(cluster~V4)+ scale_fill_manual(values=c("white","black"))+
  theme(strip.text.x = element_text(angle = 90))+facet_grid(year~V6)

###################################################################################################
# Core analysis - core PFAM annotations
###################################################################################################

kegg_map <- read.delim(file = '~/Box Sync/Papers/Phyllosphere-Function/MAG_annotation/orf-to-pfam-map.txt', sep='\t', header=FALSE)
colnames(kegg_map) <- c("OTU", "gene")
load(file="melted_phy3.RData")
head(ps) #line 213 loaded object - the melted phyloseq of metaT
ps2 <- merge(ps, kegg_map)
split1 <- str_split_fixed(ps2$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
bin = paste(split2a, split3, sep='_')
ps2$bin <- bin
#Checking to exclude M22
ps2 <- subset(ps2, bin != "M_bin.22")
ps3 <- ps2 %>% select(Sample, OTU, bin, gene, Abundance)
f <- ddply(ps3, .(Sample, bin, gene), summarise, MEAN=mean(Abundance))
f2 <- ddply(f, .(bin, gene), summarise, SUM_MEAN=sum(MEAN))
f3 <- f2 %>% spread(bin, SUM_MEAN)
rownames(f3) <- f3$gene
f3$gene <- NULL
f3[is.na(f3)] <- 0
f3[f3 > 0] <- 1
f3$SUM = rowSums(f3)
f4 <- f3[f3$SUM > 10,] #25 out of 40 bins
#> dim(f4)
#[1] 186  41
#186 out of 4070 PFAMs are in more than 10 bins, broad functions
hist(rowSums(f3))
ggplot(f4, aes(x=SUM)) + geom_histogram()
targets <- rownames(f4) # List of PFAMs that are considered core > 10

pfam_to_kegg = read.delim(file="~/Box Sync/Papers/Phyllosphere-Function/MAG_annotation/pfam-to-kegg.txt", sep='\t', , header=FALSE)
pfam_to_kegg2 <- pfam_to_kegg[pfam_to_kegg$V1 %in% targets,]
pfam_merge <- merge(pfam_to_kegg2, f4) #ORF annotations that match the targets
pfam_merge2 <- subset(pfam_merge, pfam_merge$V3 != "NA")
colnames(pfam_merge2)
MAG_sum = rowSums(pfam_merge2[4:43])
pfam_merge2$MAG_sum <- MAG_sum
ggplot(NULL, aes(x=pfam_merge2$MAG_sum, fill=pfam_merge2$V3)) + geom_histogram()+xlab("Number of MAGs out of 41 MAGs")+ ylab("Total PFams")
#summary(f3$SUM)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   2.762   3.000  39.000 

pfam_merge3 = pfam_merge2[pfam_merge2$SUM > 38,]
write.table(pfam_merge3, file="core_ORFS.txt", quote=FALSE, sep="\t")
sort(summary(as.factor(pfam_merge3$V3)))
head(pfam_merge3)
head(ps2)
f <- ddply(ps2, .(Sample, month_range, gene), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(month_range, gene), summarise, SUM_MEAN=mean(SUM))
f2b <- f2 %>% spread(month_range, SUM_MEAN)
foo <- merge(pfam_merge3, f2b, by.x = "V1", by.y = "gene")
write.table(foo,  file="core_ORFS2.txt", quote=FALSE, sep="\t")

#####################################################
#Supp Figure 4 - Core ORF functions 
#####################################################

mags=read.csv("fmags.csv", header=TRUE)
library(tidyverse)
mags=as_tibble(mags)


orfs=read.csv("core_ORFS2.csv", header=TRUE)
orfs=as_tibble(orfs)

#summarize by V3 classification (N=20)
head(orfs)

CCP<-filter(orfs, V3 == "Cellular community - prokaryotes")
terp<-filter(orfs, V3 == "Metabolism of terpenoids and polyketides")
hist(orfs$ratio)

highratios<-filter(orfs, ratio >20)

library(ggplot2)

ggplot(data = terp,aes(x=V1, y=ratio))+
  geom_bar(stat="identity")+
  ggtitle("Metabolism of terpenoids and polyketides")+
  coord_flip()

ggplot(data = orfs)+
  geom_bar(mapping=aes(x=V1, y=ratio), stat="identity")+
  #not working:  different color for ratios >20)
  #geom_segment(aes(x=highratios[,V1],y=ratio))+
  facet_wrap(~V3, nrow=4)+
  labs(title="ORF transcripts detected in 39/40 MAGs", y="Late:early normalized transcript ratio", x="ORF")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ggsave("FigureS2_ORFs_TranscriptRatios_v2.eps",
       plot=last_plot(),
       device="eps",
       dpi="print"
)


###################################################################################################
# Specific Functions of interest - pyruvate, terpenes, and isoprenes
# Figure S5
# Figure 7
###################################################################################################
load("melted_phy4-line224.RData") #as ps
pyruvate <- scan("../MAG_annotation/refinement/CoA-list.txt", what="", sep="\n") # ORFs with CoA associated in annotations
terpen_orfs = scan("terpen-orfs.txt", what ="", sep="\n") # ORFs associated with 'terpen' in annotations

foo_threonine <- subset(ps, V2 %in% pyruvate)
foo_threonine <- subset(ps, V2 %in% terpen_orfs)

ps3b <- foo_threonine
ps3b_gt0 <- subset(ps3b, ps3b$Abundance > 0)
unique(ps3b_gt0$OTU)
f <- ddply(ps3b, .(Sample, V4, V8, month, year), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(month, year, V8), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
#f2 <- subset(f2, V7 != "NA")
#f2 <- subset(f2, V7 != "Global and overview maps")
temp = subset(f2, month == "9")
#f2$zero_flag = "0"
#f2[f2$MEAN != 0,]$zero_flag = "1"
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", shape="year"))
p2 = p+geom_point(stat="identity") + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(text = element_text(angle = 0, size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_grid(year~V8, labeller = label_wrap_gen()) + theme(strip.text.x= element_text(angle = 90, size=12))+
  scale_shape_manual(values=c(16,2))+ylab('Coverage Normalized by HKGs')
ggsave('terpenes.png') #Figure 7A 
ggsave(p2, file="7A.eps", device="eps")

#checking into terpenoid backbone biosyntehsis within terpenes
ps3c <- subset(ps3b, V8 == "Terpenoid backbone biosynthesis")
ps3c_gt0 <- subset(ps3c, ps3c$Abundance > 0)
f <- ddply(ps3c, .(Sample, V4, month, year), summarise, SUM=sum(Abundance))
f <- subset(f, f$V4 != "NA")
f2 <- ddply(f, .(month, year, V4), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
#f2 <- subset(f2, V7 != "NA")
#f2 <- subset(f2, V7 != "Global and overview maps")
temp = subset(f2, month == "9")
#f2$zero_flag = "0"
#f2[f2$MEAN != 0,]$zero_flag = "1"
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity") + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(text = element_text(angle = 0, size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(angle = 90))+facet_grid(year~V4)+
  scale_shape_manual(values=c(16,2))+ylab('Coverage Normalized by HKGs')
ggsave('terpenoid_backbone_biosynthesis.png')

# Pulling out the iosprene synthesis genes by ORFs
f <- ddply(ps3c, .(Sample, OTU, V4, month, year), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(month, year, OTU, V4), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$month <- as.character(f2$month)
f2$year <- as.character(f2$year)
f2b <- subset(f2, V4 %in% c("LytB protein [PF02401.19]", "GcpE protein [PF04551.15]"))
#f2 <- subset(f2, V7 != "NA")
#f2 <- subset(f2, V7 != "Global and overview maps")
temp = subset(f2b, month == "9")
#f2$zero_flag = "0"
#f2[f2$MEAN != 0,]$zero_flag = "1"
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2b, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3) + theme_bw()+geom_errorbar(limits, width=0)+theme_bw()+xlab('Sampling Month') + ylab('Coverage Normalized by HKGs')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(angle = 90))+facet_grid(V4~OTU)+ylab('Coverage Normalized by HKGs')

# Pulling out the iosprene synthesis genes by MAGs
mag_phy = read.csv("../Revision1_29June2021/Datasets/mag_annotations.csv", row.names = 1)
ps3d <- subset(ps3c, ps3c$V4 %in% c("LytB protein [PF02401.19]", "GcpE protein [PF04551.15]"))
split1 <- str_split_fixed(ps3d$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, ".", 5)[,5]
bin = paste(split2a, split4, sep='')
ps3d$bin <- bin
ps3e <- merge(ps3d, mag_phy, by.x = "bin", by.y = "row.names")
f <- ps3e %>% select(Sample, bin, Genus, OTU, Abundance, year, V4,  month)
f2 <- ddply(f, .(bin, year, Genus, month, V4, Sample), summarise, SUM=sum(Abundance))
f3 <- ddply(f2, .(bin, month, Genus, V4, year), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f3$bin = factor(f3$bin, levels=c("M102", "M111", "M8", "M12", "M32", "M52", "M1", "M86", "M105", "M67", "M9", "M109", "S28"))
levels(f3$bin)=c("M102_Methyl", "M111_Methyl", "M8_Methyl", "M12_Methyl", "M32_Methyl", "M52_Methyl", "M1_Frigori", "M86_Pseudokineo", "M105_Microbac", "M67_Amnibac", "M9_Hymeno", "M109_Sphingo", "S28_Pseudo")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3$month <- as.character(f3$month)
f3$year <- as.character(f3$year)
p = ggplot(f3, aes_string(x="month", y="MEAN", shape="year"))
p2 = p+geom_point(stat="identity", size=3)+theme_bw()+geom_errorbar(limits, width=0)+theme(text = element_text(size=12), strip.text.x = element_text(angle = 90))+facet_grid(V4~bin)+
  scale_shape_manual(values=c(16,2))+ylab('Coverage Normalized by HKGs')
ggsave('isoprene_synthesis_by_mag.png')
ggsave(p2, file='7b.eps', device="eps")
mag_phy = read.csv("../Revision1_29June2021/Datasets/mag_annotations.csv", row.names = 1)
ps3d <- subset(ps3c, ps3c$V4 %in% c("LytB protein [PF02401.19]", "GcpE protein [PF04551.15]"))
split1 <- str_split_fixed(ps3d$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, ".", 5)[,5]
bin = paste(split2a, split4, sep='')
ps3d$bin <- bin
ps3e <- merge(ps3d, mag_phy, by.x = "bin", by.y = "row.names")
f <- ps3e %>% select(Sample, bin, Order, OTU, Abundance, year, V4,  month)
f2 <- ddply(f, .(bin, year, Order, month, V4, Sample), summarise, SUM=sum(Abundance))
f3 <- ddply(f2, .(bin, month, Order, V4, year), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f3$bin = factor(f3$bin, levels=c("M102", "M111", "M8", "M12", "M32", "M52", "M1", "M86", "M105", "M67", "M9", "M109", "S28"))
levels(f3$bin)=c("M102_Methyl", "M111_Methyl", "M8_Methyl", "M12_Methyl", "M32_Methyl", "M52_Methyl", "M1_Frigori", "M86_Pseudokineo", "M105_Microbac", "M67_Amnibac", "M9_Hymeno", "M109_Sphingo", "S28_Pseudo")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3$month <- as.character(f3$month)
f3$year <- as.character(f3$year)
p = ggplot(f3, aes_string(x="month", y="MEAN", shape="year"))
p+geom_point(stat="identity", size=3)+theme_bw()+geom_errorbar(limits, width=0)+facet_grid(V4~bin)

###################################################################################################
#Comparison of HKG MetaG and MetaT
###################################################################################################

f_hkg_metag
f_hkg_metat
f <- ddply(f_hkg_metat, .(plant, year, month), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f$month <- as.character(f$month)
f$month <- factor(f$month, levels=c("5", "6", "7", "8", "9", "10", "11"))
f$year <- as.character(f$year)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f, aes_string(x="month", y="MEAN", color="year"))
p+geom_point(stat="identity", size=3)+theme_bw()+geom_errorbar(limits, width=0)

###################################################################################################
#Functional Annotations of Clusters 
# Figure S3
###################################################################################################
load("melted_phy4-line224.RData") #as ps
split1 <- str_split_fixed(ps$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, ".", 5)[,5]

bin = paste(split2a, split4, sep='')
ps$bin <- bin
f <- ps %>% select(Sample, bin, OTU, Abundance, year, V7, month)
f2 <- ddply(f, .(V7, bin, year, month, Sample), summarise, SUM=sum(Abundance))
f3 <- ddply(f2, .(V7, bin, month, year), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3b <- subset(f3, V7 !=  "NA")
f3b$month <- as.character(f3b$month)
p = ggplot(f3b, aes_string(x="bin", y="MEAN", color="month"))
p+geom_point(stat="identity", size=3)+theme_bw()+geom_errorbar(limits, width=0)+facet_grid(~V7)

#bring in cluster
map_cluster = read.delim(file="draft_clustering.txt", sep="\t")
colnames(map_cluster) = c("bin", "phylum", "cluster")
map_cluster
ps3 <- merge(ps, map_cluster, by = "bin")
f <- ps3 %>% select(Sample, bin, OTU, Abundance, year, V7, month, cluster)
f2 <- ddply(f, .(V7, cluster, year, month, Sample), summarise, SUM=sum(Abundance))
f3 <- ddply(f2, .(V7, cluster, month, year), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f3b <- subset(f3, V7 !=  "NA")
f3b <- subset(f3b, V7 != "Global and overview maps")
f3b$month <- as.character(f3b$month)
f3b$year <- as.character(f3b$year)
f3b$zero_flag = "0"
f3b[f3b$MEAN > 0,]$zero_flag = "1"
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3b, aes_string(x="month", y="MEAN", color="cluster"))
p2= p+geom_point(stat="identity", size=3, shape=21, aes(fill=zero_flag))+theme_bw()+geom_errorbar(limits, width=0)+facet_grid(cluster~V7)+
    theme(strip.text.x = element_text(angle = 90)) + scale_fill_manual(values=c("white","black"))+ylab("Coverage Normalized by HKGs")
ggsave(p2, file="figure_s3.eps", dev = "eps", dpi = 600, width=12, height=6) #Supp Figure S3


###################################################################################################
# Comparison of late and early seasons (ratio analysis)
# Figure S4
###################################################################################################

# Stats on the KEGG Pathways of Difference between Months
#kruskal-wallis test 
#save(ps3, file="metat_melted_clusters_physeq.RData")
f <- ps3 %>% select(Sample, OTU, Abundance, year, V7, month, month_range, cluster)
save(ps3, file="metat_melted_clusters_physeq.RData")
f2 <- ddply(f, .(V7, cluster, year, month_range, Sample), summarise, SUM=sum(Abundance))

ss <- unique(f2$V7)
ss <- ss[1:22]
ss_stats <- data.frame(matrix(NA, nrow = length(ss), ncol = 2))
rownames(ss_stats) <- ss
for (i in 1:length(ss)){
  ss_data <- subset(f2, V7 == ss[i])
  res <- kruskal.test(SUM ~ month_range, data=ss_data)
  ss_stats[i,] <- c(res$statistic, res$p.value)
}

# 
# > ss_stats[ss_stats$X2 < 0.01,]
# X1           X2
# Amino acid metabolism                        81.06917 2.179542e-19
# Biosynthesis of other secondary metabolites  85.23758 2.645858e-20
# Carbohydrate metabolism                      64.27918 1.079810e-15
# Cell motility                                53.12555 3.128957e-13
# Cellular community - prokaryotes             20.78895 5.127795e-06
# Drug resistance: antimicrobial               73.36141 1.079564e-17
# Energy metabolism                            59.95990 9.680962e-15
# Folding, sorting and degradation             43.07078 5.279496e-11
# Global and overview maps                     57.49353 3.390937e-14
# Glycan biosynthesis and metabolism           46.91108 7.428152e-12
# Lipid metabolism                             94.58105 2.352594e-22
# Membrane transport                           43.38143 4.504435e-11
# Metabolism of cofactors and vitamins         60.74533 6.495701e-15
# Metabolism of other amino acids              93.88442 3.344971e-22
# Metabolism of terpenoids and polyketides    119.29176 9.040464e-28
# Nucleotide metabolism                       100.87253 9.809748e-24
# Replication and repair                       72.98320 1.307584e-17
# Signal transduction                          52.22339 4.953272e-13
# Transcription                                37.58999 8.729132e-10
# Translation                                  21.64182 3.286091e-06
# Xenobiotics biodegradation and metabolism    88.49130 5.105892e-21


f <- ps3 %>% select(Sample, OTU, Abundance, year, V7, month, month_range, cluster)
f2 <- ddply(f, .(V7, cluster, year, month_range, Sample), summarise, SUM=sum(Abundance))
f2b <- ddply(f2, .(V7, month_range), summarise, MEAN=mean(SUM))
f2b %>% spread(month_range, MEAN)

f <- ps3 %>% select(Sample, OTU, Abundance, year, V4, V7, month, month_range, cluster)
f$zero_flag = 0
f[f$Abundance > 0,]$zero_flag = 1
f2 <- ddply(f, .(OTU, month_range, V4, V7), summarise, num_samples_gtz = sum(zero_flag), num_samples_total = length(Abundance), MEAN = mean(Abundance), SE=sd(Abundance)/sqrt(length(Abundance)))
f2$percent <- f2$num_samples_gtz/f2$num_samples_total

f3 <- f2 %>% select(OTU, month_range, MEAN, percent, V4, V7)     
f4 <- f2 %>% select(OTU, month_range, MEAN)            
f5 <- f4 %>% spread(month_range, MEAN)
f5$ratio <- f5$`2`/f5$`1`
f6 <- subset(f5, ratio != "Inf")
f7 <- subset(f5, ratio == "Inf")
f6b <- merge(f6, f3, by = "OTU")

# > length(unique(ps$OTU))
#[1] 25244
# 25244 Unique OTUs in metaT (above abundance and sample threshold)
#> dim(f6)
#[1] 9823    4
#> dim(f7)
#[1] 15421     4
#>  count(f6$ratio > 1)
#[1] 6655
# 61% (15421) are only in late and not present early
# 26% (6655) are enriched in late more so than early

enriched = f6b[f6b$ratio > 1,]
enriched_uniq = f6[f6$ratio > 1,]
summary(factor(enriched$V7))[order(summary(factor(enriched$V7)))]
# Mainly CArb Tran aEnergy and amino 

only_late = merge(f7, f3, by = "OTU")
summary(factor(only_late$V7))[order(summary(factor(only_late$V7)))]
f1 = as.data.frame(summary(factor(only_late$V7)))
f2 = as.data.frame(summary(factor(enriched$V7)))
f3 <- merge(f1, f2, by = 'row.names')
colnames(f3) = c("row", "late", "enriched")
rownames(f3) = f3$row
f3$row <- NULL
f3b <- f3 %>% mutate(late_percentage=((late/sum(late))))
f3b <- f3b %>% mutate(enriched_percentage=((enriched/sum(enriched))))
f3b$ann <- rownames(f3b)
f3c <- f3b %>% select(late_percentage, enriched_percentage)
f3c$ann <- rownames(f3c)
f3c <- subset(f3c, ann != "NA's")
f3c <- subset(f3c, ann != "Global and overview maps")
p = ggplot(f3c)
p+geom_point( aes_string(x="ann", y="enriched_percentage"), colour="blue",stat="identity", size=3)+geom_point(aes_string(x="ann", y="late_percentage"), colour="red", stat="identity", size=3)+theme_bw()+
    xlab("KEGG Metabolic Pathway") + ylab("Proporation of ORFs") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


###################################################################################################
#Persistence and Enrichment - No Annotation based
###################################################################################################
load(file="phy3_line192.RData") #loads as ps
f <- ps  %>% select(Sample,  OTU, Abundance, year, month, month_range)
f$month <- as.character(f$month)
f$year <- as.character(f$year)
f$zero_flag = 0
f[f$Abundance > 0,]$zero_flag = 1

f2 <- ddply(f, .(OTU, month_range), summarise, num_samples_gtz = sum(zero_flag), num_samples_total = length(Abundance), MEAN = mean(Abundance), SE=sd(Abundance)/sqrt(length(Abundance)))
f2$percent <- f2$num_samples_gtz/f2$num_samples_total
f3 <- f2 %>% select(OTU, month_range, MEAN) 
f4 <- f2 %>% select(OTU, month_range, percent)
f5 <- f3 %>% spread(month_range, MEAN)
f5b <- f4 %>% spread(month_range, percent)
f5$ratio <- f5$`2`/f5$`1`
f6 <- subset(f5, ratio != "Inf")
f7 <- subset(f5, ratio == "Inf")
enriched_uniq = f6[f6$ratio > 1,]
head(enriched_uniq)
enriched2 <- merge(enriched_uniq, f5b, by = "OTU")
colnames(enriched2) = c("OTU", "early", "late", "ratio", "early_per", "late_per")
p = ggplot(enriched2, aes_string(x="ratio", y="late_per"))
p+geom_point()+xlab("Ratio Median Coverage of Late to Early")+ylab("Proportion of late season samples")+theme_bw()
selected <- enriched2[enriched2$ratio > 13.8 & enriched2$late_per > .564,]

enriched3 <- merge(f7, f5b, by = "OTU")
colnames(enriched3) = c("OTU", "early", "late", "ratio", "early_per", "late_per")
p = ggplot(enriched3, aes_string(x="late", y="late_per"))
p+geom_point()+xlab("Median bp Abundance in Late Season")+ylab("Proportion of late season samples")+theme_bw()

load("metat_melted_clusters_physeq.RData") #loads as ps3
ann <- ps3 %>% select(V2, OTU, V4, V5, V6, V7, V8, cluster, bin)
colnames(ann) = c("OTU", "OTU2", "V4", "V5", "V6", "V7", "V8", "cluster", "bin")
foo <- unique(merge(selected, ann, by ="OTU"))
p = ggplot(foo, aes_string(x="ratio", y="late_per", color="V7"))
p+geom_point()+xlab("Ratio Median Coverage of Late to Early")+ylab("Proportion of late season samples")+theme_bw()
write.table(foo, file="annotations_of_enriched_ORFS.txt", quote=FALSE, sep="\t")

summary(enriched3)
enriched3b <- subset(enriched3, late > 3.574e-4 & late_per > 0.3914)
foo2 <- unique(merge(enriched3b, ann, by = "OTU"))
p = ggplot(foo2, aes_string(x="late", y="late_per", color="V7"))
p+geom_point()+xlab("Ratio Median Coverage of Late to Early")+ylab("Proportion of late season samples")+theme_bw()
write.table(foo2, file="annotations_of_late_only_ORFS.txt", quote=FALSE, sep="\t")

###################################################################################################
#Other MetaGs
# Figure 8
###################################################################################################


y <- read.csv(file="other-metags-summary.txt", sep='\t', row.names=1)
dim(y)
abundance <- otu_table(as.matrix(y), taxa_are_rows=TRUE)
meta <- read.csv("sample-list.txt", sep="\t", header=TRUE, row.names = 1)
phy_meta <- sample_data(meta)
phy_iowa <- phyloseq(abundance,phy_meta)

x <- scan("./list-of-hkg-contigs.txt", what="", sep="\n")
hkg_phy <- prune_taxa(x, phy_iowa)
phy_melt_hkg = psmelt(hkg_phy)
f <- ddply(phy_melt_hkg, .(Sample), summarise, SUM=sum(Abundance))
f_hkg_metag <- f
f_hkg <- f %>% select(Sample, SUM)
abundance <- as.data.frame(t(otu_table(phy_iowa)))
abundance$Sample <- rownames(abundance)
hkg_foo <- merge(abundance, f_hkg, by="Sample")
rownames(hkg_foo) <- hkg_foo$Sample
hkg_foo$Sample <- NULL
hkg_norm <- hkg_foo/hkg_foo$SUM
hkg_norm$SUM <- NULL
abundance <- otu_table(as.matrix(t(hkg_norm)), taxa_are_rows=TRUE)
phy_iowa_norm <- phyloseq(phy_meta, abundance)  #this is the hkg normalized phyloseq object
phy_iowa <- phy_iowa_norm

#save(phy_iowa, file="other-metags-phy.RData")

phy_iowa  <- prune_samples(sample_sums(phy_iowa)>0, phy_iowa)
phy_iowa  <- filter_taxa(phy_iowa , function(x) sum(x) > 0, TRUE)
phy_iowa_melt <- psmelt(phy_iowa)
phy_iowa_melt <- subset(phy_iowa_melt, Plant != "Corn")
phy_iowa_melt$Year <- as.character(phy_iowa_melt$Year)
split1 <- str_split_fixed(phy_iowa_melt$OTU, "_", 3)[,3]
split2a <- str_split_fixed(split1, '_', 2)[,1]
split2b <- str_split_fixed(split1, '_', 2)[,2]
split3 <- str_split_fixed(split2b, '_', 2)[,1]
split4 <- str_split_fixed(split3, ".", 5)[,5]
bin = paste(split2a, split4, sep='')
phy_iowa_melt$bin <- bin
f=ddply(phy_iowa_melt, .(Sample, Year.State, Plant, bin), summarise, MEAN=mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f$zero_flag = "0"
f[f$MEAN > 0,]$zero_flag = "1"
f$bin <- factor(f$bin, levels=order_of_dendrogram_labels)
p = ggplot(f, aes(x=bin, y=MEAN,  color=Year.State))+geom_point(stat="identity", size=2)+geom_errorbar(limits, width=0)+geom_jitter()
p = p +theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0, vjust=0.5, size=10))+facet_grid(~Plant)
p2 = p + xlab("Phyllosphere MAG") + ylab("Average Abundance (Reads Mapped)")+scale_y_continuous(limits=c(0, 0.15))
ggsave(p2, file="figure_8.eps", dev = "eps", dpi = 600, width=12, height=6) #Supp Figure S3

# Figure 8
