#!/bin/bash

#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=10:00:00


#Add these two lines to activate the correct conda env
source /rds/general/user/shemming/home/anaconda3/etc/profile.d/conda.sh
conda activate longread_qc


#Variables list
ISOLATE=C105_1
INPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-results/live/Sam/${ISOLATE}/${ISOLATE}_canu_assembly
OUTPUT_DIR=/rds/general/user/shemming/projects/fisher-aspergillus-results/live/Sam/${ISOLATE}/reports

#quast command
quast -o ${OUTPUT_DIR}/quastplot ${INPUT_DIR}/${ISOLATE}_canu_assembly.contigs.fasta
