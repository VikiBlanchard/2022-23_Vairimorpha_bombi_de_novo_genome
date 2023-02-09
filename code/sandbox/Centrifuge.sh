#!/bin/bash 

"""Assembling genome from .fastq files generated during Oxford Nanopore sequencing."""

__appname__ = 'genome_assembler.sh'
__author__ = 'Viki Blanchard (viki.blanchard.2018@live.rhul.ac.uk)'
__version__ = '0.0.1'

######################################
### Code adapted from Sam Hemmings ###
######################################

# Activate conda environment to run the analysis in and install centrifuge
conda create Vairimorpha_bombi_genome
""" Dependencies currently installed in the MicrosporidiaGenomics env - need to reinstall dependencies to this env to use"""

# Make directory for log files 
mkdir Reports/logs

###########################
### Specify the isolate ###
###########################

# This script should be bashed with the isolate name as an argument after
#isolate=$1

# Or set the isolate variable manually
isolate="Vairimorpha_bombi"

#####################################
### 1.) Filter high quality reads ###
#####################################

# Set path directory for data input files 
InputDirectory="$( echo $(pwd))/Data/$isolate/fastq_pass"

# Set path directory for lambda reference genome
Lambda_Genome_Directory="$( echo $(pwd))/Data"

# Make directory for high quality reads
mkdir Results/${isolate}_assembly
ResultsDir="$( echo $(pwd))/Results/${isolate}_assembly"

# Store names of all fastq files in the input directory
ls $InputDirectory/*.fastq* > $InputDirectory/${isolate}_fastq_list.txt
txt=$(cat $InputDirectory/${isolate}_fastq_list.txt)

# Trim Nanopore adapters and remove reads mapping to the lambda phage control genome we added from the fastq files output one fastq file of the high quality reads
for file in ${txt}; do
    porechop -i ${file} -o ${ResultsDir}/porechop_$(basename "${file}") -t 20
	zcat ${ResultsDir}/porechop_$(basename "${file}") | \
	NanoLyse --reference ${Lambda_Genome_Directory}/lambda_genome.fasta --logfile "$( echo $(pwd))/Reports/logs/NanoLyselog.txt" | NanoFilt -q 5 -l 1000 | \
	gzip >> ${ResultsDir}/${isolate}_high_qual_reads.fastq.gz
done

#########################################################
### Make summary plots of filtered Nanopore sequences ###
#########################################################

NanoPlot -t 32 --fastq ${ResultsDir}/${isolate}_high_qual_reads.fastq.gz --loglength -o "$( echo $(pwd))/Reports"/NanoPlot_reports --prefix ${isolate}_high_qual_ --plots dot --format pdf --huge


##################################################
### Make summary plots of non-bombus sequences ###
##################################################

NanoPlot -t 32 --fastq "$( echo $(pwd))/Results/3.align_to_bumblebee/Vairimorpha_bombi_high_qual_reads.fastq-unaligned.fastq" --loglength -o "$( echo $(pwd))/Reports"/NanoPlot_reports --prefix ${isolate}_non-bombus_DNA_ --plots dot --format pdf --huge

#######################
### BLAST sequences ###
#######################

# Parse FASTQ into a FASTA file
gunzip -c ${ResultsDir}/${isolate}_high_qual_reads.fastq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${ResultsDir}/${isolate}.fa

# Make directory for BLAST hits 
mkdir -p ${ResultsDir}/${isolate}/blast

# BLAST multiple sequences on their remote server
blastn -db nt -query ${ResultsDir}/${isolate}.fa -out "BLAST_Query_Sample_8.1.csv" -remote 

#########################################################
### Search FASTA file for Vairimorpha bombi sequences ###
#########################################################

# Make directory for primer searches
mkdir ${ResultsDir}/Primer_Searches

# Get each sequence into one line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${ResultsDir}/${isolate}.fa > ${ResultsDir}/Primer_Searches/${isolate}_single_line.fa

# Remove the first empty line
tail -n +2 ${ResultsDir}/Primer_Searches/${isolate}_single_line.fa > ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa

# Extract the target sequences containing Vairimorpha bombi primer sequences (including their names) and save it into a corresponding file 

# Make directory for primer search results
mkdir ${ResultsDir}/Primer_Searches/Search_Results

## SSUrRNA-f1
LC_ALL=C grep -B 1 CACCAGGTTGATTCTGCCT ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_SSUrRNA-f1_matches.txt

## SSUrRNA-r1c
LC_ALL=C grep -B 1 GTTACCCGTCACTGCCTTG ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_SSUrRNA-r1c_matches.txt

## Nbombi-SSU-Jf1 
LC_ALL=C grep -B 1 CCATGCATGTTTTTGAAGATTATTAT ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_Nbombi-SSU-Jf1_matches.txt

## Nbombi-SSU-Jr1
LC_ALL=C grep -B 1 CATATATTTTTAAAATATGAAACAATAA ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_Nbombi-SSU-Jr1_matches.txt

## Napis-SSU-Jf1
LC_ALL=C grep -B 1 CCATGCATGTCTTTGACGTACTATG ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_Napis-SSU-Jf1_matches.txt

## Napis-SSU-Jr1
LC_ALL=C grep -B 1 GCTCACATACGTTTAAAATG ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_Napis-SSU-Jr1_matches.txt

## ITS-f2
LC_ALL=C grep -B 1 GATATAAGTCGTAACATGGTTGCT ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_TS-f2_matches.txt

## ITS-r2
LC_ALL=C grep -B 1 CATCGTTATGGTATCCTATTGATC ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_ITS-r2_matches.txt

######################################
### Remove contaminating sequences ###
######################################

bash remove_contamination.sh "$isolate"

#######################
### Assemble genome ###
#######################

# Make folder for output files
mkdir Results/${isolate}_assembly/canu_assembly
AssemblyDir="$( echo $(pwd))/Results/${isolate}_assembly/canu_assembly"

# Run canu assembler on reads from species non-specific, quality filtered fastq file

canu -d ${AssemblyDir}/${isolate}_canu_assembly -p ${isolate}_canu_assembly genomeSize=50m -nanopore-raw ${ResultsDir}/${isolate}_high_qual_reads.fastq.gz

#########################################
### Assess quality of genome assembly ###
#########################################

quast -o Results/quastplot ${AssemblyDir}/${isolate}_canu_assembly.contigs.fasta