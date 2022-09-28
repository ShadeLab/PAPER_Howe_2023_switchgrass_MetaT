
The scripts here include the quality control filtering.  In brief, JGI quality filtered pairs were filtered with Trimmomatic and plant host and fungal reads were subsequently removed.  

Methods in text:

We proceeded with bioinformatic analysis (Figure 2) of 192 metagenome (Dataset 1) and 78 metranscriptome (Dataset 2) observations that met JGI standards for raw data quality based on the Illumina proprietary software. We used Trimmomatic (v0.39) (Bolger et al., 2014) to remove adaptors and filter low-quality reads from fastq files using the following arguments: PE -phred33 ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36. After assembly, plant host reads were filtered out (for both metagenomes and metatranscriptomes) by removing all reads that mapped to the switchgrass genome (Panicum virgatum v1.0, DOE-JGI, http://phytozome.jgi.doe.gov/) or miscanthus genome (Miscanthus sinensis V7.1 http://phytozome.jgi.doe.gov/) using bowtie2 (v 2.4.1), samtools (v 1.13), and bedtools (v2.30.0) (Li et al., 2009; Quinlan and Hall, 2010; Langmead and Salzberg, 2012) (Figure S1). To remove the fungal reads and improve the prokaryotic signal in metagenome samples, we also filtered reads against 7 fungal genomes (Gostinčar et al., 2014; Druzhinina et al., 2018; Gill et al., 2019; Haridas et al., 2020) that represent close relatives of the most abundant fungal species in this system that we assessed and reported in our prior work (Bowsher et al., 2020) (Table S1). The genomes of these fungal species were retrieved from the JGI Genome Portal. 

----


Also here is an example of how the read counts that were mapped to each sample were parsed from Bowtie2 mapping log files.

Methods in text:

We also estimated the total reads associated with each MAG for metagenomes and metatranscriptomes (Dataset 4).
