#!/bin/bash 
## Usage: 1.filter_high_quality_reads.sh [options] isolate_name, PATH/TO/raw_reads, PATH/TO/lambda_genome.fa
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

Cleaned, filtered, ginsi partial trim merged 
* Sequences were filtered (Prequal), aligned (G-INS-i), filtered (Divvier), trimmed (trimAL), and merged (partial sequences belonging to the same taxon)



# Convert fastq to fasta for prequal filtering 
seqtk seq -a "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq.gz" > "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fa"

# Filter non-homologous DNA
prequal "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fa"

# Run multiple alignments 
G-INS-i
mafft --globalpair --maxiterate 1000 input_file > output_file
or
ginsi input_file > output_file

# Predict homology
./divvier -mincol 4 -divvygap myfile.fas

# Trim alignments 
trimal

#Merge partial sequences belonging to the same taxon
