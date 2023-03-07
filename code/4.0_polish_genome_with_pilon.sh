#!/bin/bash 
## Usage: 4.0_polish_genome_with_pilon.sh 
## Inputs: isolate_name, PATH/TO/new_assembly.fasta, PATH/TO/short_reads.fasta
## Outputs: 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

##########################
### Run polishing step ### 
##########################

# Store variables
assembly_dir="$1"
polishing_round="$2"
short_reads=$3

# Index the genome and perform short read mapping
bash ../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh ${assembly_dir} ${short_reads} unpaired

# Run Pilon using max memory available
pilon --genome ${assembly_dir} --unpaired $short_reads-mem.sorted.bam --output "Polish_loop_${polishing_round}.fasta" --verbose --vcf --variant --changes --fix all --chunksize 30000

############################
### Generate window plot ###
############################

# Generate pileup
samtools mpileup -f ${assembly_dir} -s $short_reads-mem.sorted.bam -o "Polish_loop_${polishing_round}.pileup"

# Generate windows dataframe
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_for_VCFs_mpileups_or_tabs3.pl' -r ${assembly_dir} -v "Polish_loop_${polishing_round}.fasta.vcf" -p "Polish_loop_${polishing_round}.pileup" -w 1000 -n 140 > 1000_windows-Polish_loop_${polishing_round}-normalised.tab

# Make windows plot
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_dataframe_to_R_figure3.pl' -w 1000_windows-Polish_loop_${polishing_round}-normalised.tab