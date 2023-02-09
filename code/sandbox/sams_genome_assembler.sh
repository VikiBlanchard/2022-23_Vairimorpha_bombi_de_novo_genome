#!/bin/bash 

"""Assemble genome from .fastq files generated during Oxford Nanopore sequencing."""

__appname__ = 'genome_assembler.sh'
__author__ = 'Viki Blanchard (viki.blanchard.2018@live.rhul.ac.uk)'
__version__ = '0.0.1'

################################################################
### Code from Sam Hemmings - assembling from nanopore output ###
################################################################

###Sam's steps#
# 1) In all scripts change ISOLATE=Name of isolate 
# 2) Make a .txt file called "${ISOLATE}_fastq_list.txt" containing all of the names of the fastq files that are in the your rawdata repository
# 3) run qsub highqual_reads.filter.sh and check nanoplot graphs made at the end 
# 4) run after changing any parameters you name (e.g. isolate name, genome size) qsub canu_assemby
# 5) run qsub_quast.sh to find graphs showing quality of assembly

###########################
### Specify the isolate ###
###########################

read -p "Please enter isolate name: " isolate
echo $isolate

#####################################
### 1.) Filter high quality reads ###
#####################################

# Set path directory for data input files 
InputDirectory="$( echo $(pwd))/Data/$isolate/fastq_pass"

# Make directory for high quality reads
mkdir Results/${isolate}_assembly

# Trim Nanopore adapters and remove reads from control genome
for file in "$InputDirectory"/*; do 
    echo "$(basename $file)"
    porechop -i $file -o $InputDirectory/"$(basename $file)" -t 20 #remove adapters from Oxford Nanopore reads
    zcat $InputDirectory/"$(basename $file)" #pull contents of all files together into one place
    NanoLyse --reference $InputDirectory/../../lambda_genome.fasta | NanoFilt -q 10 -l 1000 | \ #remove reads mapping to the lambda phage control genome we added from the fastq files 
	gzip >> Results/${isolate}_high_qual_reads.fastq.gz #output one fastq file of the high quality reads
done

######################################
### Remove contaminating sequences ###
######################################

bash remove_contamination.pl

#######################
### Assemble genome ###
#######################

# Make folder for output files
mkdir Results/${isolate}_assembly/canu_assembly

canu -d Results/${isolate}_canu_assembly -p ${isolate}_canu_assembly -nanopore ${input_dir}/${ISOLATE}_high_qual_reads.fastq.gz

#########################################
### Assess quality of genome assembly ###
#########################################

quast -o ${OUTPUT_DIR}/quastplot ${INPUT_DIR}/${ISOLATE}_canu_assembly.contigs.fasta