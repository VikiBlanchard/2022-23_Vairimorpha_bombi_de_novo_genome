#!/bin/bash 
## Usage: 3.1.run_canu_assembler.sh <isolate_name, pipeline_step, PATH/TO/fastq-unaligned.fastq 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

################################
### Generate a canu assembly ###
################################

# Activate conda environment to run the analysis in 
source ~/anaconda3/etc/profile.d/conda.sh
conda activate canu_assembly

# Make folder for output files
mkdir "results_${1}/${2}_${1}_assembly"
assembly_dir="results_${1}/${2}_${1}_assembly"
#assembly_dir="results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly"

# move into output folder
cd ${assembly_dir}

# Run canu assembler on reads from species non-specific, quality filtered fasta or fastq file
canu -d ./ -p "${2}_${1}" genomeSize=10m corOutCoverage=100 corMhapSensitivity=high corMinCoverage=0 minInputCoverage=0 -nanopore-raw ${3}

# move back to base directory 
cd ../..

#########################################
### Assess quality of genome assembly ###
#########################################

# Make directory to store report 
mkdir reports_${1}/${2}_${1}_assembly_quastplot

# Assess assembly quality with Quast
quast.py ${assembly_dir}/${2}_${1}.contigs.fasta -o reports_${1}/${2}_${1}_quastplot_report

# Deactivate conda environment 
conda deactivate
