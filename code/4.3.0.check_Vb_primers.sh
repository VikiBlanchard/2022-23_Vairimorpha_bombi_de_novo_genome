#!/bin/bash 
## Usage: 4.3.0.check_Vb_primers.sh isolate_name, PATH/TO/new_assembly.fasta
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

###################
### Import data ###
###################

non_bombus_assembly=$2
isolate = ${1}
non_bombus_assembly="results/canu_assemblies/Non_Bombus/test-Vb.contigs.fasta"
microsporidia_sequences="results_${isolate}/4.isolated_microsporidia_sequences/microsporidia_sequences.fa"

###########################
### Find primer matches ###
###########################

# First assembly (without Bombus terrestris sequences)
bash code/4.3.1.FASTA_search_Vb_primers.sh ${isolate} $non_bombus_assembly 

# Microsporidia-filtered sequences
bash code/4.3.1.FASTA_search_Vb_primers.sh ${isolate} $microsporidia_sequences

#######################
### Total sequences ###
#######################

grep ">" $microsporidia_sequences | wc -l
grep ">" $non_bombus_assembly | wc -l

#######################################
### Compare V. bombi primer matches ###
#######################################

# First assembly (without Bombus terrestris sequences)
grep ">" "results_${isolate}/4.isolated_microsporidia_sequences/Primer_Searches/$(basename "$non_bombus_assembly")_V_bombi_primer_matches.txt" | wc -l

# Microsporidia-filtered sequences
grep ">" "results_${isolate}/4.isolated_microsporidia_sequences/Primer_Searches/$(basename "$microsporidia_sequences")_V_bombi_primer_matches.txt" | wc -l

#######################################
### Compare V. apis primer matches ###
#######################################

# First assembly (without Bombus terrestris sequences)
grep ">" "results_${isolate}/4.isolated_microsporidia_sequences/Primer_Searches/$(basename "$non_bombus_assembly")_V_apis_primer_matches.txt" | wc -l

# Microsporidia-filtered sequences
grep ">" "results_${isolate}/4.isolated_microsporidia_sequences/Primer_Searches/$(basename "$microsporidia_sequences")_V_apis_primer_matches.txt" | wc -l
