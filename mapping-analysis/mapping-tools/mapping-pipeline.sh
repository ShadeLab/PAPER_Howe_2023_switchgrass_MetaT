#!/bin/bash -login
#PBS -l walltime=50:00:00,nodes=01:ppn=8,mem=80gb
#PBS -q main
#PBS -M adina.chuang@gmail.com
#PBS -m abe
#PBS -A ged

MY_PATH=/mnt/home/howead/metags-darte-round2-env
MAPPER_PATH=/mnt/home/howead/software/bowtie2-2.2.9
SAMTOOLS_PATH=/mnt/home/howead/software/samtools-1.3.1
THREADS=8
PAR_PATH=/mnt/home/howead/software/parallel-20120122/src
CONTIGS_FILE=/mnt/home/howead/metags-darte-round2-env/manure-assembly-all/manure.contigs.fa
CONTIGS_FILE_B=`basename $CONTIGS_FILE`

cd $MY_PATH
mkdir mapping-data
cd mapping-data
ln -s $CONTIGS_FILE .
$MAPPER_PATH/bowtie2-build $CONTIGS_FILE $CONTIGS_FILE_B
cd $MY_PATH
for x in *R1*gz; do echo $MAPPER_PATH/bowtie2 -x mapping-data/$CONTIGS_FILE_B -1 $x -2 ${x%*R1*}*R2*gz -S mapping-data/${x%*_R1*}.sam; done > map.sh
cat map.sh | $PAR_PATH/parallel
for x in mapping-data/*sam; do echo "$SAMTOOLS_PATH/samtools view -b -S $x -t mapping-data/$CONTIGS_FILE_B > $x.bam"; done > sam.sh
cat sam.sh | $PAR_PATH/parallel
for x in mapping-data/*bam; do echo "$SAMTOOLS_PATH/samtools sort $x -o $x.sorted"; done > samtools.sort.sh
cat samtools.sort.sh | $PAR_PATH/parallel
for x in mapping-data/*sorted; do echo "$SAMTOOLS_PATH/samtools index $x"; done > samtools-index.sh
cat samtools-index.sh | $PAR_PATH/parallel
for x in mapping-data/*sorted; do echo "$SAMTOOLS_PATH/samtools idxstats $x > $x.idxstats"; done > samtools-idx.sh
cat samtools-idx.sh | $PAR_PATH/parallel

module load bedtools
python coverage-bed-reference.py $CONTIGS_FILE > $CONTIGS_FILE.bed

for x in mapping-data/*sorted; do echo "bamToBed -i $x > $x.bed"; done > bamtobed.sh
cat bamtobed.sh | $PAR_PATH/parallel
for x in mapping-data/*bed; do echo "coverageBed -a $CONTIGS_FILE.bed -b $x -d > $x.bed2"; done > coveragebed.sh
cat coveragebed.sh | $PAR_PATH/parallel
for x in mapping-data/*bed2; do echo "python bedcoverage-to-coverage.py $x > $x.counts"; done > bedcoveragefinal.sh
cat bedcoveragefinal.sh | $PAR_PATH/parallel
