#!/bin/bash 
## Usage: 7.0_polish_genome_with_pilon.sh 
## Inputs: isolate_name, PATH/TO/new_assembly.fasta, PATH/TO/short_reads.fasta
## Outputs: 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)


# Store variables
assembly_dir="$1"
polishing_round="$2"
short_reads='/media/vlb19/Expansion1/Vairimorpha_bombi_sequences/Short_read_sequences/V_bombi_combined_short_reads.fastq'


# Generate pileup
samtools mpileup -f '../3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered-renamed.fasta' -s "V_bombi_combined_short_reads.fastq-mem.sorted.bam" -o "V_bombi_window_plot.pileup"

# Generate windows dataframe
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_for_VCFs_mpileups_or_tabs3.pl' -r '../3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered-renamed.fasta' -v 'Pilon_polished_assembly.vcf' -p "V_bombi_window_plot.pileup" -w 1000 -n 140 > 1000-windows_dataframe-normalised.tab

perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_for_VCFs_mpileups_or_tabs3.pl' -r '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta' -v 'Pilon_polished_assembly.vcf' -p "V_bombi_window_plot.pileup" -w 1000 > 1000-windows_dataframe_not_normalised.tab

# Make windows plot
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_dataframe_to_R_figure3.pl' -w 1000-windows_dataframe-normalised.tab

perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_dataframe_to_R_figure3.pl' -w 1000-windows_dataframe_not_normalised.tab
