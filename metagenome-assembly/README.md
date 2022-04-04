Associated Methods in Paper:

Two metagenome assemblies, one for switchgrass and one for miscanthus, were created based on metagenomes collected in 2016.  These filtered metagenome reads were combined and used for co-assembly with MEGAHIT (v 1.2.9) using--kmin-1pass (low sequencing depth) and--presets meta-large (complex metagenome) (Li et al., 2015).   

metagenome-assembly.sh - Script to perform Megahit assembly with associated parameters from filtered reads.

Additionally, we curated metagenome assembled genomes (MAGs) from the 2016 switchgrass and miscanthus metagenome libraries (n = 136 metagenomes) using Metabat (v2.2.15) (Kang et al., 2019). MAG assemblies were performed using filtered reads from switchgrass and miscanthus sampled from 2016, separately, to maximize completeness and reliability, as also done in other studies (Nayfach et al., 2019).  To assess the MAGs quality and completeness, we used CheckM (v.1.13 with the lineage_wf option) estimates of quality and completeness (Parks et al., 2015).  

BinGenomes.sb - Script that shows how binning was performed on reads with Metabat and then annotated with CheckM.
