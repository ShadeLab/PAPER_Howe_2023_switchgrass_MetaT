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
ASSEMBLY_TYPE="fullAssembly"
GLBRC=/mnt/research/ShadeLab/GLBRC
META_DIR=$GLBRC/mapping/$READTYPE/$ASSEMBLY_TYPE
SCRATCH_DIR=$SCRATCH/$READTYPE

#bowtie2 output dirs
BAMS=$SCRATCH_DIR/$ASSEMBLY_TYPE/bams
FUNGALBAMS=$BAMS/toFungal
HOSTBAMS=$BAMS/toHost
SAMS=$SCRATCH_DIR/$ASSEMBLY_TYPE/sams
FUNGALSAMS=$SAMS/toFungal

#fastqs dirs
CLEANED_FASTQS=$SCRATCH_DIR/cleaned_fastqs
PAIRED=$SCRATCH_DIR/paired
TRIMMED=$SCRATCH_DIR/trimmed
UNPAIRED=$GLBRC/mapping/$READTYPE/unpaired
FUNGAL_REMOVED=$SCRATCH_DIR/fungalCleaned

#Log files dir
HCLEANING_STATS=$META_DIR/hostRemovalFlagstats
FCLEANING_STATS=$META_DIR/fungalRemovalFlagstats

#Executables dir
SAMPLE_SCRIPTS=$SCRATCH_DIR/scripts/$READTYPE

#################################### House keeping run Variables ####################################
HEADER=$GLBRC/scripts/hpc_scripts/MappingHeader.sb #Header for individual samples to be run on the HPC
cd $UNPAIRED  #Starting Location

# HPC params
THREADS=20    ## of threads for trimming and alignment
MEM="100G"	  ## HPC RAM for data prep


SWITCHGRASS=$SCRATCH/bowtieDB/Pvirgatum_516_v5.0.fa
MISCANTHUS=$SCRATCH/bowtieDB/Msinensis_497_v7.0.fa
FUNGAL=$SCRATCH/bowtieDB/CombinedFungalAssembly.fa
KEGGFiles="$GLBRC/annotations/metagenomics-workshop/reference_db/kegg"

files=(*.gz)   #Get the sample files to process
nsamples=${#files[@]} #Number of samples to process
counter=0 #Counter to keep track of what sample we are processing

cd $SCRATCH_DIR
echo

#################################### For Each Sample Build HPC Script and Launch ####################################
for fastq in "${files[@]}"; do 
	counter=$((counter + 1))
	sample=${fastq/\.fastq\.gz/}
	
	#################################### Step 1. Make a batch script file for each sample and fill out sample specific info ####################################
	#Go to the scratch dir
	cat $HEADER >$SAMPLE_SCRIPTS/$sample.sb
	echo -en "cd $SCRATCH_DIR\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	echo "echo \"Processing sample: $sample\""  >> $SAMPLE_SCRIPTS/$sample.sb
	echo "echo \"Staring at: \$(date)\"" >> $SAMPLE_SCRIPTS/$sample.sb
	echo -en "newgrp ShadeLab \n\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	echo -en "set -e \n" >> $SAMPLE_SCRIPTS/$sample.sb
	
	if [[ $sample == *"G5"* ]]; then
		SAMFILE=$HOSTSAMS/$sample.SWGRASS.sam
 		BAMFILE=$HOSTBAMS/$sample.SWGRASS.bam
		HNAME="SWGRASS"
		HOST=$SWITCHGRASS
 	else
		SAMFILE=$HOSTSAMS/$sample.MISCANTS.sam
 		BAMFILE=$HOSTBAMS/$sample.MISCANTS.bam
		HOST=$MISCANTHUS
		HNAME="MISCANTS"
	fi
	
	#################################### Step 2. Separate combined reads into separate files (PE1, PE2, SE) ####################################
	splitting=false
	if [ ! -f $PAIRED/$sample.pe1.fastq.gz ]; then
		echo "$sample is being split"
		echo -en "split-paired-reads.py --gzip -1 $PAIRED/$sample.pe1.fastq.gz -2 $PAIRED/$sample.pe2.fastq.gz $UNPAIRED/$fastq \n\n" >>$SAMPLE_SCRIPTS/$sample.sb
		echo "echo \"\$(date) Completed Splitting Reads\"" >>$SAMPLE_SCRIPTS/$sample.sb
		splitting=true
	fi
	
	################################ Step 3. Trim adapters and QC reads #################################### 
	if [ ! -f $TRIMMED/$sample.pe1.fastq.gz ]; then
		echo "echo \"\$(date) Trimming started.\""          >>$SAMPLE_SCRIPTS/$sample.sb
		echo -en "trimmomatic PE -phred33 -threads $THREADS $PAIRED/$sample.pe1.fastq.gz $PAIRED/$sample.pe2.fastq.gz $TRIMMED/$sample.pe1.fastq.gz $TRIMMED/$sample.se1.fastq.gz $TRIMMED/$sample.pe2.fastq.gz $TRIMMED/$sample.fastq.se2.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>$TRIMSTATS/$sample.E.log \n\n" >>$SAMPLE_SCRIPTS/$sample.sb
		echo "echo \"\$(date) Completed Trimming Reads\"" >>$SAMPLE_SCRIPTS/$sample.sb
	fi
	
	################################# Step 4. Remove the host related reads #################################### 
	# 4a) bowtie2 mapping against host assembly
	if [ ! -f $HOSTSAMS/$sample.$HNAME.sam ]; then
		echo -en "bowtie2 --threads $THREADS -x $HOST -1 $TRIMMED/$sample.pe1.fastq.gz -2 $TRIMMED/$sample.pe2.fastq.gz -S $HOSTSAMS/$sample.$HNAME.sam 2>$HCLEANING_STATS/$sample.stat \n"  >> $SAMPLE_SCRIPTS/$sample.sb
	fi
	
	if [ ! -f $BAMFILE ]; then
		echo "echo \"\$(date) sam -> bam.\""          >>$SAMPLE_SCRIPTS/$sample.sb
		echo "samtools view -bS $SAMFILE > $BAMFILE" >>$SAMPLE_SCRIPTS/$sample.sb
		echo "echo \"\$(date) Completed Alignment to plant assembly.\""          >>$SAMPLE_SCRIPTS/$sample.sb
	fi
	
	# 4b) filter to get reads that are unmapped to the plant assemblies
	if [ ! -f $HOSTBAMS/$sample.unmapped.bam ]; then
		echo "echo \"\$(date) Filter HOST reads.\""          >>$SAMPLE_SCRIPTS/$sample.sb
		echo "samtools view -b -f 12 -F 256 $BAMFILE > $HOSTBAMS/$sample.unmapped.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
	fi

	# 4c) split paired-end reads into separated fastq files .._R1 .._R2
	if [ ! -f $CLEANED_FASTQS/$sample.R1.fastq ]; then
		echo "echo \"\$(date) Starting HOST sort and bam ->fastq extraction.\"" >>$SAMPLE_SCRIPTS/$sample.sb
		echo "samtools sort --threads $THREADS -n $HOSTBAMS/$sample.unmapped.bam -o $HOSTBAMS/$sample.unmapped_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
		echo -en "bedtools bamtofastq -i $HOSTBAMS/$sample.unmapped_sorted.bam -fq $CLEANED_FASTQS/$sample.R1.fastq -fq2 $CLEANED_FASTQS/$sample.R2.fastq\n\n"  >>$SAMPLE_SCRIPTS/$sample.sb
	fi

	################################# Step 5. Remove the fungal related reads #################################### 
	# 5a) bowtie2 mapping against fungal assemblies
	if [ ! -f $FUNGALBAMS/$sample.fungal_sorted.bam ]; then
		echo "echo \"\$(date) Starting Fungal align and sort.\""          >>$SAMPLE_SCRIPTS/$sample.sb
		echo -en "bowtie2 --threads $THREADS -x $FUNGAL -1 $CLEANED_FASTQS/$sample.R1.fastq -2 $CLEANED_FASTQS/$sample.R2.fastq -S $FUNGALSAMS/$sample.fungal.sam 2>$FCLEANING_STATS/$sample.fungal.stat  \n\n" >> $SAMPLE_SCRIPTS/$sample.sb
		echo "samtools view -bS $FUNGALSAMS/$sample.fungal.sam > $FUNGALBAMS/$sample.fungal.bam" >>$SAMPLE_SCRIPTS/$sample.sb
		echo "samtools sort --threads $THREADS -n $FUNGALBAMS/$sample.fungal.bam -o $FUNGALBAMS/$sample.fungal_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
		echo "echo \"\$(date) Completed Alignment to the fungal assemblies/\""          >>$SAMPLE_SCRIPTS/$sample.sb
	fi
	
	# 5b) filter to get reads that are unmapped to the fungal assemblies
	if [ ! -f $FUNGALBAMS/$sample.fungal.unmapped.bam ]; then
		echo "echo \"\$(date) Filter fungal reads.\""          >>$SAMPLE_SCRIPTS/$sample.sb
		echo "samtools view -b -f 12 -F 256 $FUNGALBAMS/$sample.fungal_sorted.bam > $FUNGALBAMS/$sample.fungal.unmapped.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
	fi
	
	#5c) split paired-end reads into separated fastq files .._R1 .._R2
	if [ ! -f $CLEANED_FASTQS/$sample.FR1.fastq ]; then
		echo "echo \"\$(date) Starting Fungal extraction.\""          >>$SAMPLE_SCRIPTS/$sample.sb
		echo "samtools sort --threads $THREADS -n $FUNGALBAMS/$sample.fungal.unmapped.bam -o $FUNGALBAMS/$sample.fungal.unmapped_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
		echo -en "bedtools bamtofastq -i $FUNGALBAMS/$sample.fungal.unmapped_sorted.bam -fq $CLEANED_FASTQS/$sample.FR1.fastq -fq2 $CLEANED_FASTQS/$sample.FR2.fastq \n\n"  >>$SAMPLE_SCRIPTS/$sample.sb
	fi

	echo "echo \"Done \$(date)\"" >> $SAMPLE_SCRIPTS/$sample.sb
	
	HOURS=4
	if [ "$splitting" = true ] ; then
		HOURS=16
	fi
	
	# If it's the last sample to process, then add an email flag so I know when it is done
	if [ $counter = $nsamples ]; then
		echo "$SAMPLE_SCRIPTS/$sample.sb"
		sbatch --time=$HOURS:00:00 --mem=$MEM --cpus-per-task=$THREADS -e "$SCRATCH_DIR/run_logs/$sample.$READTYPE.err" -o "$SCRATCH_DIR/run_logs/$sample.$READTYPE.out" --mail-type=ALL --mail-user=<USER_EMAIL> "$SAMPLE_SCRIPTS/$sample.sb"
		echo "Done"
	else 
		echo -en "$counter. $sample "
		sbatch --time=$HOURS:00:00 --mem=$MEM --cpus-per-task=$THREADS -e "$SCRATCH_DIR/run_logs/$sample.$READTYPE.err" -o "$SCRATCH_DIR/run_logs/$sample.$READTYPE.out" "$SAMPLE_SCRIPTS/$sample.sb"
	fi
	
	
done

echo -en "\n$dCount samples completed.\n\n"

exit 0