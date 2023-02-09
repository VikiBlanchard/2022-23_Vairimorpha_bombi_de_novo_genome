#!/usr/bin/env bash

# This script will create the coverage table required to obtain the cumulative graph

STRAIN1=${isolate} # Specify your species or strain name

REF1=/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/7.Vairimorpha_bombi_8.1-3_Pilon_Polish/Pilon_polished_assembly.fa.fasta # should be changed for your genome file path

TYPE1=contig # Specify your genome assembly type, such as contig, scaffold, chromosome, etc.

LEN1=`bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF1` # Size of assembled genome

LEN1=`'/home/vlb19/Documents/Coding/Downloaded_Repositories/bioawk/bioawk' -c fastx '{sum+=length($seq)}END{print sum}' $REF1` # Size of assembled genome

# Create the output file having a header line

echo "line,length,type,coverage" > length.csv

# Calculate cumulative sum and write result to the output file (HiFi data)

cat $REF1 | bioawk -c fastx -v line="$STRAIN1" '{print line","length($seq)","length($seq)}' | sort -k3rV -t "," | awk -F "," -v len="$LEN1" -v type="$TYPE1" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> length.csv

cat $REF1 | '/home/vlb19/Documents/Coding/Downloaded_Repositories/bioawk/bioawk' -c fastx -v line="$STRAIN1" '{print line","length($seq)","length($seq)}' | sort -k3rV -t "," | awk -F "," -v len="$LEN1" -v type="$TYPE1" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> length.csv

# Calculate cumulative sum and write result to the output file (CLR data)

STRAIN2=CLR_Dmel

REF2=/path/to/CLR_Dmel.contigs.fasta # should be changed your genome name

TYPE2=contig

LEN2=`bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF2`

cat $REF2 | bioawk -c fastx -v line="$STRAIN2" '{print line","length($seq)","length($seq)}' | sort -k3rV -t "," | awk -F "," -v len="$LEN2" -v type="$TYPE2" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> length.csv

# Calculate cumulative sum and write result to the output file (ONT data)

STRAIN3=ONT_Dmel

REF3=/path/to/ONT_Assembly.fasta # should be changed your genome name

TYPE3=contig

LEN3=`bioawk -c fastx '{sum+=length($seq)}END{print sum}' $REF3`

cat $REF3 | bioawk -c fastx -v line="$STRAIN3" '{print line","length($seq)","length($seq)}' | sort -k3rV -t "," | awk -F "," -v len="$LEN3" -v type="$TYPE3" 'OFS=","{ print $1,$2,type,(sum+0)/len; sum+=$3 }' >> length.csv