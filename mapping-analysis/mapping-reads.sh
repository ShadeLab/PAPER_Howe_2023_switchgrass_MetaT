#!/bin/bash -login
: '
to create an environment to run this pipeline using conda you can use the commands below
conda create -n microbe -c bioconda khmer,trimmomatic,samtools,bowtie2,microbecensus
conda activate microbe
This script also assumes that anaconda3 is installed at $USER/bin/anaconda3/bin. If it is not, 
change MappingHeader.sb to the location of your main conda install. If you already of some of these
tools installed, ensure all the packages above are installed and MappingHeader.sb will say "unable to activate microbe"
but everything will still run on your conda.
'

#################################### Paths for storage of files ####################################
READTYPE="metaG"
ASSEMBLY_TYPE='mags'
GLBRC=/home/GLBRCORG/dooley.shanek/GLBRC
META_DIR=$GLBRC/data/mapping/$READTYPE
SCRATCH_DIR=$META_DIR
DATA_DIR=$GLBRC/data

#bowtie2 output dirs
BAMS=$SCRATCH_DIR/$ASSEMBLY_TYPE/bams
SAMS=$SCRATCH_DIR/$ASSEMBLY_TYPE/sams

#fastqs dirs
CLEANED_FASTQS=$SCRATCH_DIR/cleaned_fastqs
PAIRED=$SCRATCH_DIR/paired
TRIMMED=$SCRATCH_DIR/trimmed
UNPAIRED=$META_DIR/unpaired
FUNGAL_REMOVED=$SCRATCH_DIR/fungalCleaned

#Log files dir
HCLEANING_STATS=$META_DIR/hostRemovalFlagstats
FCLEANING_STATS=$META_DIR/fungalRemovalFlagstats
MAGS_STATS=$META_DIR/magFlagstats
FLAGSTATS=$META_DIR/$ASSEMBLY_TYPE/flagstats
IDXSTATS=$META_DIR/$ASSEMBLY_TYPE/idxStats
COUNTS=$META_DIR/$ASSEMBLY_TYPE/counts
SINGLEGENES=$META_DIR/$ASSEMBLY_TYPE/singleCopyGeneCounts
TRIMSTATS=$META_DIR/$ASSEMBLY_TYPE/trimStats

#Executables dir
SAMPLE_SCRIPTS=$SCRATCH_DIR/scripts/$READTYPE

#################################### House keeping run Variables ####################################
HEADER=$GLBRC/scripts/hpc_scripts/MappingHeader.sb #Group Should be able to access this
cd $META_DIR/unpaired  #Starting Location

THREADS=20    ## of threads for trimming and alignment
# THREADS=1
MEM="100G"
# CONTIGS_FILE=$SCRATCH/bowtieDB/AnnotatedContigs.fa
# CONTIGS_FILE=$SCRATCH/bowtieDB/Final.contigs.fa

SWITCHGRASS_MetaG=$DATA_DIR/assemblies/SwitchgrassUnfiltered.fa
MISCANTHUS_MetaG=$DATA_DIR/assemblies/MiscanthusUnfiltered.fa

MISCANTS_MAG_ASSEMBLY=/home/GLBRCORG/dooley.shanek/GLBRC/data/assemblies/mags/MiscanthusMags_gt_50percComplete.fa
SWITCHG_MAG_ASSEMBLY=/home/GLBRCORG/dooley.shanek/GLBRC/data/assemblies/mags/SwitchgrassMags_gt_50percComplete.fa
MAG_ASSEMBLY=""
files=(*.gz)   #Get the sample files to process
nsamples=${#files[@]} #Number of samples to process
counter=0 #Counter to keep track of what sample we are processing

#If the dirs for output don't exist create them
mkdir -p $BAMS $SAMS $CLEANED_FASTQS $PAIRED $TRIMMED $FLAGSTATS $IDXSTATS $SAMPLE_SCRIPTS $TRIMSTATS $SINGLEGENES $COUNTS $MAGS_STATS $HCLEANING_STATS $FCLEANING_STATS 

cd $SCRATCH_DIR

HPC_SCRIPTS=/home/GLBRCORG/dooley.shanek/GLBRC/scripts/hpc_scripts
#################################### For Each Sample Build HPC Script and Launch ####################################
for fastq in "${files[@]}"; do 
	counter=$((counter + 1))
	sample=${fastq/\.fastq\.gz/}

	if [ -f $COUNTS/$sample.mags.krona.kegg.minpath.tab ]; then
		continue
	fi

	#################################### Step 1. Make a batch script file for each sample ####################################
	#Go to the scratch dir
	cat $HEADER >$SAMPLE_SCRIPTS/$sample.sh
	echo -en "cd $SCRATCH_DIR\n\n" >> $SAMPLE_SCRIPTS/$sample.sh
	echo "echo \"Processing sample: $sample\""  >> $SAMPLE_SCRIPTS/$sample.sh
	echo "echo \"Staring at: \$(date)\"" >> $SAMPLE_SCRIPTS/$sample.sh
	
	#################################### Step 2. Separate combined reads into separate files (PE1, PE2, SE) ####################################
	# echo -en "split-paired-reads.py --gzip -1 $PAIRED/$sample.fastq.pe1.gz -2 $PAIRED/$sample.fastq.pe2.gz $UNPAIRED/$fastq 2>/dev/null\n\n" >>$SAMPLE_SCRIPTS/$sample.sh
	# echo "echo \"\$(date) Completed Splitting Reads\"" >>$SAMPLE_SCRIPTS/$sample.sh

	# #################################### Step 3. Trim adapters and QC reads #################################### 
	# echo -en "trimmomatic PE -phred33 -threads $THREADS $PAIRED/$sample.fastq.pe1.gz $PAIRED/$sample.fastq.pe2.gz $TRIMMED/$sample.fastq.pe1.gz $TRIMMED/$sample.fastq.se1.gz $TRIMMED/$sample.fastq.pe2.gz $TRIMMED/$sample.fastq.se2.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>$TRIMSTATS/$sample.E.log \n\n" >>$SAMPLE_SCRIPTS/$sample.sh

	#################################### Step 4. Align reads to the respective metagenomic assembly #################################### 
	MAG_BED=""
	if [[ $sample == *"G5"* ]]; then
		PLANT="SWGRASS_MG"
 		SAMFILE=$SAMS/$sample.SWGRASS_MG.sam
 		BAMFILE=$BAMS/$sample.SWGRASS_MG.sorted.bam
 		CONTIGS_FILE=$SWITCHGRASS_MetaG
 		MAG_ASSEMBLY=$SWITCHG_MAG_ASSEMBLY
 		MAG_SOURCE="SwitchgrassMags"
 		MAG_BED=$GLBRC/data/annotations/prokkaAnnotation/SwitchgrassMags.AnnotationMap.bed
 		# cat $TRIMMED/switchgrass/$sample.fastq.se1.gz $TRIMMED/switchgrass/$sample.fastq.se2.gz >$TRIMMED/switchgrass/$sample.fastq.se12.gz
 		# echo -en "bowtie2 --threads $THREADS -x $SWITCHGRASS_MetaG -1 $TRIMMED/switchgrass/$sample.fastq.pe1.gz -2 $TRIMMED/switchgrass/$sample.fastq.pe2.gz -U $TRIMMED/switchgrass/$sample.fastq.se12.gz -S $SAMFILE >$FLAGSTATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sh
 	else
 		PLANT="MISCANTS_MG"
 		SAMFILE=$SAMS/$sample.MISCANTS_MG.sam
 		BAMFILE=$BAMS/$sample.MISCANTS_MG.sorted.bam
 		CONTIGS_FILE=$MISCANTHUS_MetaG
 		MAG_ASSEMBLY=$MISCANTS_MAG_ASSEMBLY
 		MAG_SOURCE="MiscanthusMags"
 		MAG_BED=$GLBRC/data/annotations/prokkaAnnotation/$MAG_SOURCE.AnnotationMap.bed
 		cat $TRIMMED/miscan/$sample.fastq.se1.gz $TRIMMED/miscan/$sample.fastq.se2.gz >$TRIMMED/miscan/$sample.fastq.se12.gz
 		echo -en "bowtie2 --threads $THREADS -x $MISCANTHUS_MetaG -1 $TRIMMED/miscan/$sample.fastq.pe1.gz -2 $TRIMMED/miscan/$sample.fastq.pe2.gz -U $TRIMMED/miscan/$sample.fastq.se12.gz -S $SAMFILE >$FLAGSTATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sh
	fi

	#################################### Step 5. Sort the reads and get the alignment numbers #################################### 
	echo "echo \"\$(date) Completed Alignment to the metagenomic assembly.\""          >>$SAMPLE_SCRIPTS/$sample.sh
	echo "samtools view --threads $THREADS -bt $CONTIGS_FILE.fai $SAMFILE &> $BAMFILE"         >>$SAMPLE_SCRIPTS/$sample.sh
	echo "samtools sort $BAMFILE -o $BAMS/$sample.$PLANT.sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sh
	echo "samtools index -@ $THREADS $BAMS/$sample.$PLANT.sorted.bam "                              >>$SAMPLE_SCRIPTS/$sample.sh
	echo "echo \"\$(date) Completed sorting MAG aligned reads.\""          >>$SAMPLE_SCRIPTS/$sample.sh
	
	# # #################################### Step 6. Get the read counts along the contigs #####################################
	echo "bedtools coverage -hist -a $MAG_BED -b $BAMS/$sample.$PLANT.sorted.bam > $COUNTS/$sample.map.hist" >> $SAMPLE_SCRIPTS/$sample.sh

	echo "python $GLBRC/scripts/get_coverage_for_genes.py -i <(echo $COUNTS/$sample.map.hist) > $COUNTS/$sample.mags.coverage" >> $SAMPLE_SCRIPTS/$sample.sh
	# echo -en "samtools idxstats $BAMFILE >$IDXSTATS/$sample.mags.tsv\n\n" >> $SAMPLE_SCRIPTS/$sample.sh
	

	# #Filtering EC file took the number of genes from 524.470 to 54,630
	ANNO=$GLBRC/data/annotations/prokkaAnnotation
	echo "python $GLBRC/scripts/genes.to.kronaTable.py -i $ANNO/PROKKA.$MAG_SOURCE.ec -m $GLBRC/data/annotations/kegg/ec.to.pwy -H $GLBRC/data/annotations/kegg/pwy.hierarchy -n $sample -l <(grep \"minpath 1\" $ANNO/PROKKA.$MAG_SOURCE.kegg.minpath) -c $COUNTS/$sample.mags.coverage -o $COUNTS/$sample.mags.krona.kegg.minpath.tab" >> $SAMPLE_SCRIPTS/$sample.sh

	# cat $SAMPLE_SCRIPTS/$sample.sh
	cat /home/GLBRCORG/dooley.shanek/GLBRC/scripts/hpc_scripts/SubHeader.sub >$HPC_SCRIPTS/sampleSubs/$sample.sub
	echo "executable = $SAMPLE_SCRIPTS/$sample.sh" >>$HPC_SCRIPTS/sampleSubs/$sample.sub
	echo "queue" >>$HPC_SCRIPTS/sampleSubs/$sample.sub
	condor_submit $HPC_SCRIPTS/sampleSubs/$sample.sub

	cat $SAMPLE_SCRIPTS/$sample.sh
	
	
	
done







































# export last=$SAMPLE_SCRIPTS/$sample.sb
# echo "\n" $last $SCRATCH_DIR/run_logs/$sample.$READTYPE.err


#################################### Alternative way to do the count files ####################################
#module load bedtools
#python coverage-bed-reference.py $CONTIGS_FILE > $CONTIGS_FILE.bed

#for x in mapping-data/*sorted; do echo "bamToBed -i $x > $x.bed"; done > bamtobed.sh
#cat bamtobed.sh | $PAR_PATH/parallel
#for x in mapping-data/*bed; do echo "coverageBed -a $CONTIGS_FILE.bed -b $x -d > $x.bed2"; done > coveragebed.sh
#cat coveragebed.sh | $PAR_PATH/parallel
#for x in mapping-data/*bed2; do echo "python bedcoverage-to-coverage.py $x > $x.counts"; done > bedcoveragefinal.sh
#cat bedcoveragefinal.sh | $PAR_PATH/parallel

#################################### Step 4. Remove the host related reads #################################### 
	#4a) bowtie2 mapping against host sequence
	# if [[ $sample == *"G5"* ]]; then
 # 		SAMFILE=$SAMS/$sample.SWGRASS.sam
 # 		BAMFILE=$BAMS/$sample.SWGRASS.bam
 # 		echo -en "bowtie2 --threads $THREADS -x $SWITCHGRASS -1 $TRIMMED/$sample.fastq.pe1.gz -2 $TRIMMED/$sample.fastq.pe2.gz -U $TRIMMED/$sample.fastq.se12.gz -S $SAMS/$sample.SWGRASS.sam >$HCLEANING_STATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sb
 # 	else
 # 		SAMFILE=$SAMS/$sample.MISCANTS.sam
 # 		BAMFILE=$BAMS/$sample.MISCANTS.bam
 # 		echo -en "bowtie2 --threads $THREADS -x $MISCANTHUS -1 $TRIMMED/$sample.fastq.pe1.gz -2 $TRIMMED/$sample.fastq.pe2.gz -U $TRIMMED/$sample.fastq.se12.gz -S $SAMS/$sample.MISCANTS.sam >$HCLEANING_STATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	# fi
	# echo "samtools view -bS $SAMFILE > $BAMFILE" >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "echo \"\$(date) Completed Alignment to plant assembly.\""          >>$SAMPLE_SCRIPTS/$sample.sb

	#4b) filter to get reads that are unmapped to the plant assemblies
	# echo "samtools view -b -f 4 -F 256 $BAMFILE > $BAMS/$sample.unmapped.bam"  >>$SAMPLE_SCRIPTS/$sample.sb

	#4c) split paired-end reads into separated fastq files .._R1 .._R2
	# echo "samtools sort -n $BAMS/$sample.unmapped.bam -o $BAMS/$sample.unmapped_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "bedtools bamtofastq -i $BAMS/$sample.unmapped_sorted.bam -fq $CLEANED_FASTQS/$sample.R1.fastq -fq2 $CLEANED_FASTQS/$sample.R2.fastq"  >>$SAMPLE_SCRIPTS/$sample.sb

	# echo "bowtie2 --threads $THREADS -x $FUNGAL -1 $CLEANED_FASTQS/$sample.R1.fastq -2 $CLEANED_FASTQS/$sample.R2.fastq -S $SAMS/$sample.fungal.sam 2>$FCLEANING_STATS/$sample.fungal.stat" >> $SAMPLE_SCRIPTS/$sample.sb
	# echo "samtools view -bS $SAMS/$sample.fungal.sam > $BAMS/$sample.fungal.bam" >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "samtools view -b -f 4 -F 256 $BAMS/$sample.fungal.bam > $BAMS/$sample.fungal.unmapped.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "samtools sort -n $BAMS/$sample.fungal.unmapped.bam -o $BAMS/$sample.fungal.unmapped_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "rm $BAMS/$sample.fungal.unmapped.bam" >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "bedtools bamtofastq -i $BAMS/$sample.fungal.unmapped_sorted.bam -fq $CLEANED_FASTQS/$sample.FR1.fastq -fq2 $CLEANED_FASTQS/$sample.FR2.fastq"  >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "echo \"\$(date) Completed Alignment to the fungal assemblies/\""          >>$SAMPLE_SCRIPTS/$sample.sb
	#################################### Step 5. Run MicrobeCensus to normalize samples by single copy gene #################################### 
	# echo -en "run_microbe_census.py -t $THREADS $CLEANED_FASTQS/$sample.FR1.fastq,$CLEANED_FASTQS/$sample.FR2.fastq $SINGLEGENES/$sample.txt\n\n" >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "echo \"Completed microbe census\"" >>$SAMPLE_SCRIPTS/$sample.sb

	################################### Step 6. Align the cleaned reads to the metagenomic assembly #################################### 
	# echo -en "bowtie2 --threads $THREADS -x $CONTIGS_FILE -1 $CLEANED_FASTQS/$sample.FR1.fastq -2 $CLEANED_FASTQS/$sample.FR2.fastq -S $SAMS/$sample.sam >$FLAGSTATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	# echo "echo \"\$(date) Completed Alignment to the metagenomic assembly.\""          >>$SAMPLE_SCRIPTS/$sample.sb
	
	#################################### Step 7. Compress the sam file to make a bam #################################### 
	# echo -en "samtools view --threads $THREADS -bt $CONTIGS_FILE.fai $SAMS/$sample.sam &>$BAMS/$sample.bam\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	# echo -en "\n\nsamtools view --threads $THREADS -bt $CONTIGS_FILE.fai $SAMS/$sample.metaG.sam &>$BAMS/$sample.metaG.bam  \n"
	# echo \"Done converting\"; \n 
	# samtools sort -n $BAMS/$sample.metaG.bam -o $BAMS/$sample.metaG_sorted.bam \n
	# samtools index -@ $THREADS $BAMS/$sample.metaG_sorted.bam \n\n" >> $SAMPLE_SCRIPTS/$sample.sb

	# Remove the sam file to save space, commented out because everything done on scratch which will auto delete after time runs out
	# echo -en "rm $SAMS/$sample.sam\n\n">> $SAMPLE_SCRIPTS/$sample.sb

	#################################### Step 8. Sort the reads #####################################
	# echo -en "samtools sort -o $BAMS/$sample.sorted.bam $BAMS/$sample.bam\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	#Remove the unsorted bam file to save space and because we have the sorted bam now, and no need to keep both
	# echo -en "rm $BAMS/$sample.bam\n\n" >> $SAMPLE_SCRIPTS/$sample.sb

	# #################################### Step 9. Index the reads #####################################
	# echo -en "samtools index -@ $THREADS $BAMS/$sample.sorted.bam\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
################################### Step 3.5. Align reads to the MAG assembly #################################### 
	# echo -en "bowtie2 --threads $THREADS -x $MAG_ASSEMBLY -1 $TRIMMED/switchgrass/$sample.fastq.pe1.gz -2 $TRIMMED/switchgrass/$sample.fastq.pe2.gz -S $SAMS/$sample.metaG_MAGs.sam 2>$MAGS_STATS/$sample.stat \n\n" >> $SAMPLE_SCRIPTS/$sample.sh
	# echo "echo \"\$(date) Completed Alignment to the metagenomic assembly.\""          >>$SAMPLE_SCRIPTS/$sample.sh
	# echo "samtools view --threads $THREADS -bt $CONTIGS_FILE.fai $SAMS/$sample.metaG_MAGs.sam &> $BAMS/$sample.metaG_MAGs.bam"         >>$SAMPLE_SCRIPTS/$sample.sh
	# echo "samtools sort $BAMS/$sample.metaG_MAGs.bam -o $BAMS/$sample.metaG_MAGs_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sh
	# echo "samtools index -@ $THREADS $BAMS/$sample.metaG_MAGs_sorted.bam "                              >>$SAMPLE_SCRIPTS/$sample.sh
	# echo "echo \"\$(date) Completed sorting MAG aligned reads.\""          >>$SAMPLE_SCRIPTS/$sample.sh
	
	# echo "Submitted $sample"
	# break
