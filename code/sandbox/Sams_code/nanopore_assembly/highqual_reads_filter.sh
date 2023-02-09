#!/bin/bash

#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=10:00:00

source /rds/general/user/shemming/home/anaconda3/etc/profile.d/conda.sh
conda activate nanopore3



###List of variables###
#Isolate name
	ISOLATE=C105_1
#Location of fastq rawdata files
	INPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-rawdata/live/minION/${ISOLATE}/fastq_pass
#Names of all the fastq files in the input directory
	TXT=$(cat /rds/general/user/shemming/home/Scripts/nanopore_assembly/${ISOLATE}_fastq_list.txt)
#Location to place filtered file
	OUTPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-results/live/Sam/${ISOLATE}
#Location of Lambda genome
	LAMBDA_GEN=/rds/general/user/shemming/home/Scripts/nanopore_assembly
#Ephermal location
	EPH=/rds/general/user/shemming/projects/fisher-aspergillus-results/ephemeral



###Loop to filter fastq files###
for i in ${TXT}
    do
    	porechop -i ${INPUT_DIR}/${i} -o ${EPH}/porechop_${i} -t 20
	zcat ${EPH}/porechop_${i} | \
	NanoLyse --reference ${LAMBDA_GEN}/lambda_genome.fasta | NanoFilt -q 10 -l 1000 | \
	gzip >> ${OUTPUT_DIR}/${ISOLATE}_high_qual_reads.fastq.gz
    done

###This section is to make plots###

source /rds/general/user/shemming/home/anaconda3/etc/profile.d/conda.sh
conda activate nanopore

NanoPlot -t 32 --fastq ${OUTPUT_DIR}/${ISOLATE}_high_qual_reads.fastq.gz\
 --loglength -o ${OUTPUT_DIR}/reports --prefix ${ISOLATE}_high_qual_ --plots dot --format pdf --huge
