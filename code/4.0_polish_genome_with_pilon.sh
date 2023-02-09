#!/bin/bash 
## Usage: 7.0_polish_genome_with_pilon.sh 
## Inputs: isolate_name, PATH/TO/new_assembly.fasta, PATH/TO/short_reads.fasta
## Outputs: 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

#######################################################
### Set up conda environment to run the analysis in ###
#######################################################

# Activate the environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pilon_polishing

# Store isolate name as variable
isolate=${1}
isolate="Vairimorpha_bombi_8.1-3"

### Set path directory for data input files:
# Path to assembly
assembly_path=${2}
assembly_path='/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta'

# Path to short read sequences
short_reads=${3}
short_reads='/media/vlb19/Expansion1/Vairimorpha_bombi_sequences/Short_read_sequences'

# Make output directory
mkdir "results_${isolate}/4.${isolate}_Pilon_Polish/"

################################
### Polish genome with Pilon ###
################################

# Combine all short reads into one file
zcat ${short_reads}/*/*/*.fastq.gz | gzip >> ${short_reads}/V_bombi_combined_short_reads.fastq.gz
gunzip ${short_reads}/V_bombi_combined_short_reads.fastq.gz

# Move short reads into results directory 
mv "${short_reads}/V_bombi_combined_short_reads.fastq" "results_${isolate}/7.${isolate}_Pilon_Polish/V_bombi_combined_short_reads.fastq"

# Align short read sequences to assembly
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${assembly_path}" "results_${isolate}/7.${isolate}_Pilon_Polish/V_bombi_combined_short_reads.fastq" unpaired


# Run pilon on the genome
pilon --genome ${assembly_path} --bam '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3/V_bombi_combined_short_reads.fastq-mem.sorted.bam'   --outdir "results_${isolate}/7.${isolate}_Pilon_Polish" --output "Pilon_polished_assembly" --verbose --vcf --variant

#pilon --genome "./3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta" --bam './V_bombi_combined_short_reads.fastq-mem.sorted.bam' --output "Pilon_polished_assembly" --verbose

conda deactivate