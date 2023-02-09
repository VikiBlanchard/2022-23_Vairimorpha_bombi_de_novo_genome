#!/bin/bash 
## Usage: 4.2.0.isolate_microsporidia_sequences.sh 
## Inputs: isolate_name, PATH/TO/new_assembly.fasta
## Outputs: 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

# Save path to new assembly into variable 
isolate=${1}
path_to_assembly=${2}
#"""if input is empty then:""" 
#path_to_assembly='results/canu_assemblies/Non_Bombus/test-Vb.contigs.fasta'

###############################################################
### Record assembly sequence IDs that contain BUSCO matches ###
###############################################################
step_4_results_dir="results_${isolate}/4.isolated_microsporidia_sequences_raw"

# Copy in busco results table
cp $step_4_results_dir/busco_analysis_results/run_microsporidia_odb10/full_table.tsv results_${isolate}/4.isolated_microsporidia_sequences_raw/busco_results_table.tsv

# Remove first two rows that disrupt the table format 
sed -i '1,2d;' results_${isolate}/4.isolated_microsporidia_sequences_raw/busco_results_table.tsv

# Pull all BUSCO-matched sequence IDs into one file
python3 code/4.2.1.get_BUSCO_matched_sequences.py "results_${isolate}/4.isolated_microsporidia_sequences_raw/" 

##################################################################
### Copy all busco-matched sequences from assembly to new file ###
##################################################################

# Get each sequence from the assembly into one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${path_to_assembly} > results_${isolate}/4.isolated_microsporidia_sequences_raw/non_bombus_assembly_single_line.fa 

# Remove the first empty line
tail -n +2 results_${isolate}/4.isolated_microsporidia_sequences/non_bombus_assembly_single_line.fa > results_${isolate}/4.isolated_microsporidia_sequences_raw/non_bombus_assembly_single_line_2.fa 

# Extract microsporidia sequences from assembly

for sequenceID in $(cat results_${isolate}/4.isolated_microsporidia_sequences_raw/busco_match_sequences.txt); do
    microsporidia_assembly_sequences=C grep -B 1 $sequenceID results_${isolate}/4.isolated_microsporidia_sequences_raw/non_bombus_assembly_single_line_2.fa >> results_${isolate}/4.isolated_microsporidia_sequences_raw/microsporidia_sequences.fa
done
