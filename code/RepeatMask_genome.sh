#!/bin/bash 
## Usage: RepeatMask_genome.sh
## Inputs: isolate_name, PATH/TO/new_assembly.fasta, PATH/TO/short_reads.fasta
## Outputs: 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

##################################
### Set up working environment ###
##################################

# Activate conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate RepeatMask_genome

###########################
### Softmask the genome ### 
###########################

### Identify novel repeats in genome assembly

# Run RepeatModeler
BuildDatabase -name V_bombi 'VABO_genome-reordered-renamed-no-Ncontig00000000-fix-final-issues-Ns-at-end-reordered-renamed.fa'
RepeatModeler -database V_bombi -pa 1 -LTRStruct > out.log

mkdir RepeatMasker_RM_lib
mkdir fungi
mkdir bacteria
mkdir viruses

# Repeat masking with RepeatMasker using repeats were found by RepeatModeler
RepeatMasker -lib V_bombi-families.fa -s -parallel 10 -dir RepeatMasker_RM_lib 'VABO_genome-reordered-renamed-no-Ncontig00000000-fix-final-issues-Ns-at-end-reordered-renamed.fa'

# Repeatmask the genome with fungal, viral, and bacterial sequences
RepeatMasker -species fungi -pa 8 -dir fungi -gff -e ncbi -s 'RepeatMasker_RM_lib/VABO_genome-reordered-renamed-no-Ncontig00000000-fix-final-issues-Ns-at-end-reordered-renamed.fa.masked'
RepeatMasker -species bacteria -pa 8 -dir bacteria -gff -e ncbi -s 'RepeatMasker_RM_lib/VABO_genome-reordered-renamed-no-Ncontig00000000-fix-final-issues-Ns-at-end-reordered-renamed.fa.masked'
RepeatMasker -species viruses -pa 8 -dir viruses -gff -e ncbi -s 'RepeatMasker_RM_lib/VABO_genome-reordered-renamed-no-Ncontig00000000-fix-final-issues-Ns-at-end-reordered-renamed.fa.masked'
