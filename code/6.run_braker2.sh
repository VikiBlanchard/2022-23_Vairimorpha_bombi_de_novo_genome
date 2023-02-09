#!/bin/bash 
## Usage: 7.run_breaker2.sh < isolate_name, PATH/TO/new_assembly.fasta
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)
source ~/anaconda3/etc/profile.d/conda.sh
conda activate genome_annotation

# Make output directory
mkdir results_${1}/6.annotation_${3}
cd results_${1}/6.annotation_${3}

####################
### Run Augustus ###
####################
# test run 1 - sucessful augustus annotation 
#augustus --species=encephalitozoon_cuniculi_GB results/canu_assemblies/Non_Bombus/test-Vb.contigs.fasta > results_${isolate}/6.annotation/${isolate}_annotation.gff
# test run 2 (troubleshooting braker2) - successful augustus annotation 
#augustus --species=encephalitozoon_cuniculi_GB test-Vb.contigs.fasta > test-Vb_augustus_annotation.gff

augustus --singlestrand=true --species=encephalitozoon_cuniculi_GB --stopCodonExcludedFromCDS=false '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1_and_8.2/3.no_bombus_Vairimorpha_bombi_8_assembly/3.no_bombus_Vairimorpha_bombi_8.contigs.fasta' --outfile='/results_Vairimorpha_bombi_8.1_and_8.2/6.annotation/no_bombus_Vairimorpha_bombi_8_assembly_augustus_annotation.gff'

################################
### Prepare protein database ###
################################

#wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
#tar xvf odb10_arthropoda_fasta.tar.gz
#cat arthropoda/Rawdata/* > proteins.fasta
#wget https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz
#tar xvf odb10_fungi_fasta.tar.gz
#cat fungi/Rawdata/* >> proteins.fasta

##mv protein.fasta 

###################
### Run braker2 ###
###################

tail -n +2 ${2} > ${isolate}.contigs.fasta
#sed -i 's/ //g' ${isolate}.contigs.fasta

# Get each sequence from the assembly into one line
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < test-Vb.contigs.fasta > non_bombus_assembly_single_line.fa 

# Microsporidia genome 
braker.pl --genome=${isolate}.contigs.fasta --prot_seq= '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/data/microsporidia_odb10/refseq_db.faa' --cores 4

# Fungal genome
braker.pl --genome=${isolate}.contigs.fasta --prot_seq='/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/data/fungi_odb10/refseq_db.faa' --cores 4

#######################
### Troubleshooting ###
#######################

#export AUGUSTUS_CONFIG_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/config
#export GENEMARK_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/gmes_linux_64
#export BAMTOOLS_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/bamtools
#export SAMTOOLS_PATH=/usr/bin/samtools

# Running through braker2 examples 
#braker.pl --genome genome.fa --bam RNAseq.bam --softmasking --cores 8
#braker.pl --genome genome.fa --prot_seq proteins.fa --prg gth --trainFromGth --softmasking --cores 4
#braker.pl --genome=test-Vb.contigs.fasta --esmode --softmasking --cores 4