#!/bin/bash

#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=10:00:00


#Add these two lines to activate the correct conda env
source /rds/general/user/shemming/home/anaconda3/etc/profile.d/conda.sh
conda activate nanopore


#List of variables
ISOLATE=C105_1
INPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-results/live/Sam/${ISOLATE}

#Command to run NanoPlot
NanoPlot -t 32 --fastq ${INPUT_DIR}/${ISOLATE}_high_qual_reads.fastq.gz\
 --loglength -o ${INPUT_DIR}/reports --prefix ${ISOLATE}_high_qual_ --plots dot --format pdf --huge
