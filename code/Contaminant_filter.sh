#!/bin/bash 
## Usage: Contaminant_filter
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)


##########################################
### Remove contaminating DNA sequences ###
##########################################

# Align all reads to the Bombus terrestris genome 
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${B_terrestris_genome}" "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq" unpaired

# Move the results to the relevant directory
mv "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq-mem.sam" "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_high_qual_reads.fastq-mem.sam"

# Convert SAM file to fasta file
perl code/2.2.SAM_unaligned_reads_to_fasta.pl "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_high_qual_reads.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fasta" 

# Convert SAM file to fastq file 
perl code/2.3.SAM_unaligned_reads_to_fastq.pl "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_high_qual_reads.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fastq" 

##################################
### Remove Ralstonia_pickettii ###
################################## 

# Align all reads to Ralstonia_wenshanesnsis_strain_56D2
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Ralstonia_pickettii/ncbi-genomes-2023-03-15/GCF_020341455.1_ASM2034145v1_genomic.fna' 'Vairimorpha_bombi_8.1-3_unaligned_high_qual_reads.fastq' unpaired > Ralstonia_pickettii_unaligned.fastq-mem.sam

# Convert SAM file to fasta file
perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl Vairimorpha_bombi_8.1-3_unaligned_high_qual_reads.fastq-mem.sam > "Ralstonia_pickettii_unaligned.fastq"

######################################
### Remove Ralstonia wenshanesnsis ###
###################################### 

# Align all reads to Ralstonia_wenshanesnsis_strain_56D2
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Ralstonia_wenshanesnsis_strain_56D2/ncbi-genomes-2023-03-15/GCF_021173085.1_ASM2117308v1_genomic.fna' "Ralstonia_pickettii_unaligned.fastq" unpaired

# Convert SAM file to fasta file
perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl Ralstonia_pickettii_unaligned.fastq-mem.sam > "Ralstonia_wenshanesnsis_strain_56D2_unaligned.fastq"


######################################
### Remove Snodgrassella_alvi_wkB2 ###
###################################### 

# Align all reads to Snodgrassella_alvi_wkB2
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Snodgrassella_alvi_wkB2/ncbi-genomes-2023-03-15/GCF_022870965.1_ASM2287096v1_genomic.fna' "Ralstonia_wenshanesnsis_strain_56D2_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Ralstonia_wenshanesnsis_strain_56D2_unaligned.fastq-mem.sam' > "Snodgrassella_alvi_wkB2_unaligned.fastq"

##################################
### Remove Ralstonia_insidiosa ###
################################## 

# Align all reads to  Ralstonia_insidiosa
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Ralstonia_insidiosa.fna' "Snodgrassella_alvi_wkB2_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Snodgrassella_alvi_wkB2_unaligned.fastq-mem.sam' > "Ralstonia_insidiosa_unaligned.fastq"

########################################
### Remove Ralstonia_mannitolilytica ###
######################################## 

# Align all reads to  Ralstonia_mannitolilytica
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Ralstonia_mannitolilytica.fna' "Ralstonia_insidiosa_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Ralstonia_insidiosa_unaligned.fastq-mem.sam' > "Ralstonia_mannitolilytica_unaligned.fastq"

########################################
### Remove Lactobacillus_panisapium ###
######################################## 

# Align all reads to  Ralstonia_mannitolilytica
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Lactobacillus_panisapium/ncbi_dataset/data/GCA_019469265.1/GCA_019469265.1_ASM1946926v1_genomic.fna' "Ralstonia_mannitolilytica_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Ralstonia_mannitolilytica_unaligned.fastq-mem.sam' > "Lactobacillus_panisapium_unaligned.fastq"

########################################
### Remove Gilliamella sp. ESL0441 ###
######################################## 

# Align all reads to  Gilliamella_sp._ESL0441
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Gilliamella_sp._ESL0441/ncbi-genomes-2023-03-15/GCF_019469185.1_ASM1946918v1_genomic.fna' "Lactobacillus_panisapium_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Lactobacillus_panisapium_unaligned.fastq-mem.sam' > "Gilliamella_sp._ESL0441_unaligned.fastq"


########################################
### Remove Gilliamella apicola ###
######################################## 

# Align all reads to  Gilliamella_apicola
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Gilliamella_sp._ESL0441/ncbi-genomes-2023-03-15/GCF_019469185.1_ASM1946918v1_genomic.fna' "Gilliamella_sp._ESL0441_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Gilliamella_sp._ESL0441_unaligned.fastq-mem.sam' > "Gilliamella_apicola_unaligned.fastq"

########################################
### Remove Cupriavidus gilardii ###
######################################## 

# Align all reads to  Cupriavidus_gilardii
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Cupriavidus_gilardii.fa' "Gilliamella_apicola_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Gilliamella_apicola_unaligned.fastq-mem.sam' > "Cupriavidus_gilardii_unaligned.fastq"

########################################
### Remove Cupriavidus basilensis ###
######################################## 

# Align all reads to  Cupriavidus_basilensis
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Cupriavidus_basilensis.fa' "Cupriavidus_gilardii_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Cupriavidus_gilardii_unaligned.fastq-mem.sam' > "Cupriavidus_basilensis_unaligned.fastq"

#########################################
### Remove Achromobacter xylosoxidans ###
#########################################

# Align all reads to  Achromobacter_xylosoxidans
bash ../../../../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh '/media/vlb19/Expansion1/reference_sequences/20230315_Contaminant_list/Achromobacter_xylosoxidans.fa' "Cupriavidus_basilensis_unaligned.fastq" unpaired

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Cupriavidus_basilensis_unaligned.fastq-mem.sam' > "Achromobacter_xylosoxidans_unaligned.fastq"

perl ../../../../../code/2.2.SAM_unaligned_reads_to_fasta.pl 'Cupriavidus_basilensis_unaligned.fastq-mem.sam' > "Achromobacter_xylosoxidans_unaligned.fasta"

#################################
### Save final filtered reads ###
#################################

cp "Achromobacter_xylosoxidans_unaligned.fasta" "../V.bombi_contaminant_filtered.fasta"