#!/bin/bash

#PBS -l select=1:ncpus=32:mem=64gb
#PBS -l walltime=10:00:00

#Add these two lines to activate the correct conda env
source /rds/general/user/shemming/home/anaconda3/etc/profile.d/conda.sh
conda activate longread_assembly

#List variables
ISOLATE=C105_1
INPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-results/live/Sam/${ISOLATE}
OUTPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-results/live/Sam/${ISOLATE}

#Canu command
canu -d ${OUTPUT_DIR}/${ISOLATE}_canu_assembly -p ${ISOLATE}_canu_assembly \
	gridOptions="-lselect=1:ncpus=32:mem=64gb -lwalltime=10:00:00" \
	genomeSize=29m -nanopore ${INPUT_DIR}/${ISOLATE}_high_qual_reads.fastq.gz
