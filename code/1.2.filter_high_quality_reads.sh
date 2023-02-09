#!/bin/bash 
## Usage: 1.filter_high_quality_reads.sh [options] isolate_name, PATH/TO/raw_reads, PATH/TO/lambda_genome.fa
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
conda activate read_quality_filter

# Store isolate name as variable
isolate=${1}
#isolate="Vairimorpha_bombi_8.1-3"

### Set path directory for data input files:
# Raw reads from fastq_pass nanopore output
raw_reads_dir=${2}
#raw_reads_dir='/media/vlb19/Expansion/Vairimorpha_bombi_sequences/Combined_fastq_passes'
# Lambda reference genome
lambda_genome_dir=${3}
#lambda_genome_dir="$( echo $(pwd))/data"

### Initialise output directories

# High quality reads
mkdir results_${isolate}/1.${isolate}_high_qual_reads
results_dir="results_${isolate}/1.${isolate}_high_qual_reads"
# Log files 
mkdir reports_${isolate}/logs
mkdir reports_${isolate}/logs/high_qual_read_filter_log
log_output_dir="reports_${isolate}/logs/high_qual_read_filter_log"
# Sequences without barcodes
mkdir ${results_dir}/no_barcodes


##############################################
### Remove contaminant barcoding sequences ###
##############################################

# Store number of input fastq files 
total_input_fastqs=`ls ${raw_reads_dir}/*.fastq* | wc -l`

# Parse barcodes from read data 
files_complete=0
for fastq_file in ${raw_reads_dir}/*.fastq*; do
    # Bin reads by presence of barcodes with porechop 
    porechop -i ${fastq_file} -b ${results_dir} --barcode_threshold 85 --untrimmed -t 20

    # Print progress report
    files_complete=$((files_complete+1))
    echo "$files_complete / $total_input_fastqs files checked for barcodes"
done

# Get sequence IDs for all the reads with barcodes 
for barcode_file in ${results_dir}/B*.fastq*; do 
    # Store list of sequence IDs in text file 
    zgrep '^[@]' ${barcode_file} >> ${results_dir}/barcode_sequences.txt 
done

# Delete generated temp fastq files
rm ${results_dir}/*.fastq*

# Remove barcode sequences from fastq input files
files_complete=0
for fastq_file in ${raw_reads_dir}/*.fastq*; do

    # Scan raw fastq for barcode IDs and omit matching sequences from new fastq file  
    zcat ${fastq_file} | python3 code/remove_reads_from_fastq.py ${results_dir}/barcode_sequences.txt | gzip - > ${results_dir}/no_barcodes/$(basename -- "${fastq_file}")

    # Print progress report 
    files_complete=$((files_complete+1))
    echo "$files_complete / $total_input_fastqs files searched"
done


#################################
### Filter high quality reads ###
#################################

files_complete=0
for fastq_file in ${results_dir}/no_barcodes/*.fastq*; do
    
    # Trim Nanopore adapters 
    porechop -i ${fastq_file} -o ${results_dir}/porechop_$(basename "${fastq_file}") -t 20
    
    # Remove reads mapping to the lambda phage control genome 
	zcat ${results_dir}/porechop_$(basename "${fastq_file}") | NanoLyse --reference ${lambda_genome_dir}/lambda_genome.fasta --logfile ${log_output_dir}/Nanolyse_log.txt | NanoFilt -q 5 -l 1000 | gzip >> ${results_dir}/${isolate}_high_qual_reads.fastq.gz

    # Print progress report
    files_complete=$((files_complete+1))
    echo "$files_complete / $total_input_fastqs files checked for barcodes"
    
done

# Unzip the zipped fastq
gunzip ${results_dir}/${isolate}_high_qual_reads.fastq.gz

#########################################################
### Make summary plots of filtered Nanopore sequences ###
#########################################################

NanoPlot -t 32 --fastq ${results_dir}/${isolate}_high_qual_reads.fastq \
--loglength -o "$( echo $(pwd))/reports_${isolate}"/NanoPlot_reports --prefix ${isolate}_high_qual_ \
--plots dot --format pdf --huge

# Deactivate conda environment 
conda deactivate