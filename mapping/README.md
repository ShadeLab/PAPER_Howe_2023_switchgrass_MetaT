Methods in text:

Read recruitment was performed with MAGs that meet either medium- and/or high- quality standards (Swan et al., 2013; Bowers et al., 2017). The metagenome abundance of contigs in MAGs in each sample was estimated based on the median coverage of filtered metagenome reads associated to each MAG contig. Specifically, Bowtie2 (v2.2.2) was used to align reads to all focal MAG contigs (using default setting and allowing for a single read to map to only one reference).  Bedtools (v2.28) was used to estimate the coverage of each basepair within the contig.  The estimated abundance of contig was based on the median basepair coverage of all reads mapped to the contig, and the estimated abundance of a MAG was based on the average median coverage of all contig within its bins.  The metatranscriptome abundance was estimated based on protein-encoding genes identified in MAGs.  For each MAG contig, open reading frames (ORFs) and functional genes were identified using Prodigal (v2.6.3, default parameters).  Transcripts were mapped to ORFs associated with each MAG to estimate median base coverage of each ORF (Bowtie2, default parameters, no multiple mappings allowed).  

Scripts in this folder:

mapping.sh - This is the aggregated scripts to perform the analyses described in the method.  It references specific scripts contained in the folder in mapping-scripts.sh, which helps obtain the specific median bp coverage from the Bowtie2 output files.

This mapping strategy and coverage estimation are used to estimate also housekeeping gene abundances per sample.  Similarly, these are the same model scripts used to identify host and fungal reads that were removed to improve the prokaryotic signal in metagenomes.
