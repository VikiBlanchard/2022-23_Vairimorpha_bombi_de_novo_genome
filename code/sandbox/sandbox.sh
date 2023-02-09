#!/bin/bash 

####################
### SANDBOX CODE ###
####################

#####################################################


/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/bash_scripts/FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion/reference_sequences/GCF_910591885.1iyBomTerr1.2/GCF_910591885.1_iyBomTerr1.2.genome.fa' Vairimorpha_bombi_high_qual_reads.fastq unpaired

(base) vlb19@vlb19-ASUSLaptop:/media/vlb19/Expansion/2022_Vairimorpha_Genome/Results/4.assemble_with_canu$ 

canu -d ./ -p test-Vb genomeSize=10m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 -nanopore-raw /media/vlb19/Expansion/2022_Vairimorpha_Genome/Results/3.align_to_bumblebee/Vairimorpha_bombi_high_qual_reads.fastq-unaligned.fastq


#######################
### BLAST sequences ###
#######################

# Parse FASTQ into a FASTA file
gunzip -c results/1.${isolate}_high_qual_reads.fastq.gz | \
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${ResultsDir}/${isolate}.fa

# Make directory for BLAST hits 
mkdir -p ${ResultsDir}/${isolate}/blast
mkdir results/V_bombi_blast

# BLAST multiple sequences on their remote server
blastn -query /home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results/canu_assemblies/Non_Bombus/test-Vb.contigs.fasta -db nt -out results/V_bombi_blast/"BLAST_Query_Sample_8.1.csv" -remote 


################################################################
### 2.) Remove contaminating Bombus terrestris DNA sequences ###
################################################################

# Align all reads to the Bombus terrestris genome 
bash /home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/bash_scripts/FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion/reference_sequences/GCF_910591885.1iyBomTerr1.2/GCF_910591885.1_iyBomTerr1.2.genome.fa' Vairimorpha_bombi_high_qual_reads.fastq unpaired

# Convert SAM file to fasta file
perl /media/vlb19/Expansion/2022_Vairimorpha_Genome/code/perl_scripts/SAM_unaligned_reads_to_fasta.pl """SAM_file""" """output_dir"""

# Convert SAM file to fastq file 
perl /media/vlb19/Expansion/2022_Vairimorpha_Genome/code/perl_scripts/SAM_unaligned_reads_to_fastq.pl """SAM_file""" """output_dir"""

######################
### canu assembler ###
######################

canu -d ./ -p test-Vb genomeSize=10m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 -nanopore-raw /media/vlb19/Expansion/2022_Vairimorpha_Genome/Results/3.align_to_bumblebee/Vairimorpha_bombi_high_qual_reads.fastq-unaligned.fastq 
# loc of sequences is hard coded as shell is unable to symlink onto separate drive

canu -d ${AssemblyDir}/${isolate}_canu_assembly -p ${isolate}_canu_assembly genomeSize=50m -nanopore-raw ${ResultsDir}/${isolate}_high_qual_reads.fastq.gz

#########################################
### Assess quality of genome assembly ###
#########################################

# Make directory to store report 
mkdir Reports/quastplot

# Assess assembly quality with Quast
quast.py ~/Documents/Coding/2022_Vairimorpha_Genome/Results/canu_assemblies/Non_Bombus/test-Vb.contigs.fasta -o Reports/quastplot

#quast -o Results/quastplot ${AssemblyDir}/${isolate}_canu_assembly.contigs.fasta

# Assess assembly quality with mummer
# go through manual, align new genome against itself using "no simplify" - shows duplications, inversions. Outputs delta file which can be plugged into mummerplot for visualisation

########################################
### Run BUSCO analysis on new genome ###
########################################

busco -m genome -i ~/Documents/Coding/2022_Vairimorpha_Genome/Results/canu_assemblies/Non_Bombus/test-Vb.contigs.fasta  -l microsporidia_odb10 -o "busco_Vairimorpha_bombi_non-bombus"

#####################################################
### Pull out sequences that contain BUSCO matches ###
#####################################################

# Make directory for BUSCO match files
mkdir Results/matched_busco_sequences

# Pull all BUSCO matches into one file
for busco_match in results/busco_results/busco_Vairimorpha_bombi_non-bombus/run_microsporidia_odb10/busco_sequences/*/*; do
	cat ${busco_match} >> Results/matched_busco_sequences/all_busco_matched_sequences.fa
done

# Make directory for BUSCO match files
mkdir results_{1}/4.matched_busco_sequences

## Plot busco summary chars 
python3 generate_plot.py -wd 



# Run BUSCO analysis on new assembly 
busco -m genome -i $path_to_assembly -l microsporidia_odb10 -o ../../results_${1}/4.isolate_microsporidia_sequences/"busco_analysis_results"


cp results/busco_results/busco_Vairimorpha_bombi_non-bombus/run_microsporidia_odb10/full_table.tsv results_${1}/4.isolate_microsporidia_sequences/busco_results_table.tsv

##################################
### Make BUSCO summary reports ###
##################################

# Plot busco summary charts 
generate_plot.py -wd ../../results_${1}/4."busco_analysis_results"

# Send summary plot to reports directory
mv ../../results_${1}/4.busco_analysis_results/"busco_figure.png" ../"busco_figure.png"

# Send summary plot R file to code directory
mv ../../results_${1}/4.busco_analysis_results/busco_figure.R ../../code/4.1.1.busco_figure.R

# Correct R plot output path
sed -i "s/\"..\/..\/results_*\/4\.busco_analysis\_results/my_output <- paste(\"..\/reports\_${1}/g" ../../code/4.1.1.busco_figure_${1}.R

# Move back to base working directory
cd ../..

# Deactivate the conda environment
conda deactivate 

#####################################################################
### Make text file of sequences that contain microsporidia BUSCOs ###
#####################################################################

# Import BUSCO results table 
#busco_results_df = pd.read_csv('results/matched_busco_sequences/full_table.tsv', sep='\t')  

#########################################################
### Search FASTA file for Vairimorpha bombi sequences ###
#########################################################

# Get each sequence into one line
"""Need to generalise paths"""

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${ResultsDir}/${isolate}.fa > ${ResultsDir}/Primer_Searches/${isolate}_single_line.fa


awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${ResultsDir}/${isolate}.fa > ${ResultsDir}/Primer_Searches/${isolate}_single_line.fa

# Remove the first empty line
tail -n +2 ${ResultsDir}/Primer_Searches/${isolate}_single_line.fa > ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa

#########################################
### Extract primer-matching sequences ###
#########################################

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

## ITS-f2
LC_ALL=C grep -B 1 GATATAAGTCGTAACATGGTTGCT ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_TS-f2_matches.txt

## ITS-r2
LC_ALL=C grep -B 1 CATCGTTATGGTATCCTATTGATC ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_ITS-r2_matches.txt

##################################
### Check for Vairimorpha apis ###
##################################

## Napis-SSU-Jf1
LC_ALL=C grep -B 1 CCATGCATGTCTTTGACGTACTATG ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_Napis-SSU-Jf1_matches.txt

## Napis-SSU-Jr1
LC_ALL=C grep -B 1 GCTCACATACGTTTAAAATG ${ResultsDir}/Primer_Searches/${isolate}_single_line_2.fa > ${ResultsDir}/Primer_Searches/Search_Results/${isolate}_Napis-SSU-Jr1_matches.txt

#####################################################################
# Metagenome analysis tutorial 
https://carpentries-incubator.github.io/metagenomics/aio/index.html
#####################################################################