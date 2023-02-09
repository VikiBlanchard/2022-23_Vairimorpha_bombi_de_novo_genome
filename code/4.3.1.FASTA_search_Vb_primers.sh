#!/bin/bash 
## Usage: 4.3.0.check_Vb_primers.sh < isolate_name, PATH/TO/assembly.fasta
## Function: Search FASTA files for sequences matching Vairimorpha bombi primer sequences 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)


##############################
### Set-up directory paths ###
##############################

# Make directory for primer searches
mkdir results_${1}/4.isolated_microsporidia_sequences/Primer_Searches
output_dir="results_${1}/4.isolated_microsporidia_sequences/Primer_Searches"

##################
### Prep files ###
##################

# Read in data
input_assembly=$2

# Get each sequence into one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${input_assembly} > ${output_dir}/$(basename $input_assembly)_single_line.fa

# Remove the first empty line
tail -n +2 ${output_dir}/$(basename ${input_assembly})_single_line.fa > ${output_dir}/$(basename ${input_assembly})_single_line_2.fa

#######################################
### Save primer sequences as arrays ###
#######################################

declare -A V_bombi_primers
V_bombi_primers=(
    ["SSUrRNA-f1"]="CACCAGGTTGATTCTGCCT"
    ["SSUrRNA-r1c"]="GTTACCCGTCACTGCCTTG"
    ["Nbombi-SSU-Jf1"]="CCATGCATGTTTTTGAAGATTATTAT"
    ["Nbombi-SSU-Jr1"]="CATATATTTTTAAAATATGAAACAATAA"
    ["ITS-f2"]="GATATAAGTCGTAACATGGTTGCT"
    ["ITS-r2_seq"]="CATCGTTATGGTATCCTATTGATC"
)

declare -A V_apis_primers
V_apis_primers=(
    ["Napis-SSU-Jf1"]="CCATGCATGTCTTTGACGTACTATG"
    ["Napis-SSU-Jr1"]="GCTCACATACGTTTAAAATG"
)

#########################################
### Extract primer-matching sequences ###
#########################################

# Check V bombi primers 
for V_bombi_primer in "${V_bombi_primers[@]}"; do 
    V_bombi_primer_matches=C grep -B 1 $V_bombi_primer ${output_dir}/$(basename $input_assembly)_single_line_2.fa >> ${output_dir}/$(basename $input_assembly)_V_bombi_primer_matches.txt
done

# Check V. apis primers
for V_apis_primer in "${V_apis_primers[@]}"; do 
    V_apis_primer_matches=C grep -B 1 $V_apis_primer ${output_dir}/$(basename $input_assembly)_single_line_2.fa >> ${output_dir}/$(basename $input_assembly)_V_apis_primer_matches.txt
done