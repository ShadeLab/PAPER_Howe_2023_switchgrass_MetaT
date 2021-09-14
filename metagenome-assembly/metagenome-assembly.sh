#!/bin/bash

source ~/bin/anaconda3/bin/activate

cd /home/GLBRCORG/dooley.shanek/GLBRC/data/mapping/metaG/trimmed/switchgrass   

THREADS=20
MEM=100e9

echo "Starting the metagenome assembly at: ${date}"
#cat *.se1.gz &>se1s.fastq.gz
#cat *.se2.gz &>se2s.fastq.gz
#ls *.pe1.gz |tr '\n' ','>R1s.txt
#ls *.pe2.gz |tr '\n' ','>R2s.txt
megahit --kmin-1pass -1 $(cat R1s.txt) -2 $(cat R2s.txt) -r se1s.fastq.gz,se2s.fastq.gz --presets meta-large -o /home/GLBRCORG/dooley.shanek/GLBRC/data/mapping/metaG/trimmed/switchgrass/switchUnfiltered --num-cpu-threads 20 -m $MEM           

echo "Finished the metagenome assembly at: ${date}"
