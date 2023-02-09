###########################################
### Remove contaminating  DNA sequences ###
###########################################

# Create output directory
mkdir results

# Align all reads to the Ralstonia_sp genome  
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh data/Ralstonia_sp_SRR11285236.fasta results/8.1-3.fasta unpaired

# Assemble reads without Ralstonia DNA 
conda activate canu_assembly
canu -d ./ -p "No_Ralstonia_assembly" genomeSize=10m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 minInputCoverage=0 -nanopore-raw results/8.1-3.fasta

#######################
### Assemble genome ### 
#######################
# Run canu assembler on raw high quality reads
bash code/3.run_canu_assembler.sh "${isolate}" "3.No_Ralstonia" '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Ralstonia_sp./R_unaligned.fasta'

# Look at GC content of every sequence in FASTA file 
seqkit fx2tab --name --only-id --gc "results_${isolate}/3.raw_${isolate}_assembly/3.raw_$isolate.contigs.fasta"

seqkit fx2tab --name --only-id --gc 'results_Vairimorpha_bombi_8.1-3/3.No_Ralstonia_Vairimorpha_bombi_8.1-3_assembly/3.No_Ralstonia_Vairimorpha_bombi_8.1-3.contigs.fasta' | cat > "results_Vairimorpha_bombi_8.1-3/${isolate}_content_per_tig.tsv"


























###########################################
### Remove contaminating  DNA sequences ###
###########################################

# Create output directory
mkdir results_${isolate}/2.align_reads_to_Bombus_terrestris 

# Align all reads to the Bombus terrestris genome 
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${B_terrestris_genome}" "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq" unpaired

# Move the results to the relevant directory
mv "results_${isolate}/1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq-mem.sam" "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_high_qual_reads.fastq-mem.sam"

# Convert SAM file to fasta file
perl code/2.2.SAM_unaligned_reads_to_fasta.pl "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_high_qual_reads.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fasta" 

# Convert SAM file to fastq file 
perl code/2.3.SAM_unaligned_reads_to_fastq.pl "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_high_qual_reads.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fastq" 

# Run canu assembler on raw high quality reads
bash code/3.run_canu_assembler.sh "${isolate}" "3.raw" "../1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq" 


########################################
### Remove contaminating S. alvi DNA ###
########################################

# Create output directory
mkdir results_${isolate}/2.align_reads_to_sella_alvi

# Align all reads to the sella genome 
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${sella_genome}" "results_${isolate}/2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fastq" unpaired

# Move the results to the relevant directory
mv "results_${isolate}/2.align_reads_to_Bombus_terrestris/Vairimorpha_bombi_8.3_unaligned_high_qual_reads.fastq-mem.sam" "results_${isolate}/2.align_reads_to_sella_alvi/${isolate}_sella_unaligned_reads.fastq-mem.sam" 

# Convert SAM file to fasta file
perl code/2.2.SAM_unaligned_reads_to_fasta.pl '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.3/2.align_reads_to_Snodgrasella_alvi/Vairimorpha_bombi_8.3_Snodgrasella_unaligned_reads.fastq-mem.sam' > '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.3/2.align_reads_to_Snodgrasella_alvi/Vairimorpha_bombi_8.3_Snodgrasella_unaligned_reads.fasta'

# Convert SAM file to fastq file 
perl code/2.3.SAM_unaligned_reads_to_fastq.pl "results_${isolate}/2.align_reads_to_sella_alvi/${isolate}_Snodgrassella_unaligned_reads.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_sella_alvi/${isolate}_Snodgrassella_unaligned_reads.fastq" 

##############################################
### Remove contaminating Ralstonia sp. DNA ###
##############################################

# Create output directory
mkdir results_${isolate}/2.align_reads_to_Ralstonia_sp.

# bash FASTQ to FASTA and back converter 
sed -n '1~4s/^@/>/p;2~4p' 'results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Snodgrasella_alvi/Vairimorpha_bombi_8.1_and_8.2_Snodgrasella_unaligned_reads.fastq' > 8_1and2.fa
sed -n '1~4s/^@/>/p;2~4p' 'results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Snodgrasella_alvi/Vairimorpha_bombi_8.3_Snodgrasella_unaligned_reads.fastq' > 8_3.fa

cat 8_1and2.fa 8_3.fa > 8.1-3.fa 

sed -n '1~4s/^>/@/p;2~4p' 8.1-3.fa >  8.1-3.fastq


# Align all reads to the Ralstonia genome 
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${Ralstonia_sp_genome}"  8.1-3.fastq unpaired

# Move the results to the relevant directory
mv "results_${isolate}/2.align_reads_to_Snodgrasella_alvi/Vairimorpha_bombi_8.1-3_Snodgrasella_unaligned_reads.fastq-mem.sam" "results_${isolate}/2.align_reads_to_Ralstonia_sp./R_unaligned.fastq-mem.sam" 

# Convert SAM file to fasta file
perl code/2.2.SAM_unaligned_reads_to_fasta.pl "results_${isolate}/2.align_reads_to_Ralstonia_sp./R_unaligned.fastq-mem.sam"  > "results_${isolate}/2.align_reads_to_Ralstonia_sp./R_unaligned.fasta" 

# Convert SAM file to fastq file 
perl code/2.3.SAM_unaligned_reads_to_fastq.pl "results_${isolate}/2.align_reads_to_Ralstonia_sp./R_unaligned.fastq-mem.sam"  > "results_${isolate}/2.align_reads_to_Ralstonia_sp./R_unaligned.fastq" 


###########################################
### Remove contaminating  DNA sequences ###
###########################################

# Create output directory
mkdir results_${isolate}/2.align_reads_to_Cupriavidus_necator_genome

# Align all reads to the Bombus terrestris genome 
bash code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${Cupriavidus_necator_genome}" '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Ralstonia_sp./R_unaligned.fastq' unpaired

# Move the results to the relevant directory
mv '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Ralstonia_sp./R_unaligned.fastq-mem.sam' "results_${isolate}/2.align_reads_to_Cupriavidus_necator_genome/${isolate}_no_Cupriavidus_necator.fastq-mem.sam"

# Convert SAM file to fasta file
perl code/2.2.SAM_unaligned_reads_to_fasta.pl "results_${isolate}/2.align_reads_to_Cupriavidus_necator_genome/${isolate}_no_Cupriavidus_necator.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_Cupriavidus_necator_genome/${isolate}_no_Cupriavidus_necator.fasta" 

# Convert SAM file to fastq file 
perl code/2.3.SAM_unaligned_reads_to_fastq.pl "results_${isolate}/2.align_reads_to_Cupriavidus_necator_genome/${isolate}_no_Cupriavidus_necator.fastq-mem.sam" > "results_${isolate}/2.align_reads_to_Cupriavidus_necator_genome/${isolate}_no_Cupriavidus_necator.fastq" 


# Convert SAM file to fastq file 
perl code/2.3.SAM_unaligned_reads_to_fastq.pl '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Ralstonia_sp./8.1-3.fastq-mem.sam' > '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Ralstonia_sp./8.1-3.fasta'

# Run canu assembler 
bash code/3.run_canu_assembler.sh "${isolate}" "3.No_Ralstonia_sp." '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Ralstonia_sp./8.1-3.fastq'

bash code/3.run_canu_assembler.sh "${isolate}" "3.raw" "../1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq" 

#################################################################
### Generate an assembly from reads without B. terrestris DNA ###
#################################################################

# Run canu assembler on raw high quality reads
# bash code/3.run_canu_assembler.sh "${isolate}" "3.raw" "../1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq" 

# Run canu assembler on Bombus-filtered reads 
bash code/3.run_canu_assembler.sh "${isolate}" "3.no_bombus" "../2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fastq"

# Save path to this assembly to variable 
path_to_no_bombus_assembly="$( echo $(pwd))/results_${isolate}/3.no_bombus_${isolate}_assembly/3.no_bombus_${isolate}.contigs.fasta"


######################################################################################
### Generate an assembly from reads without B. terrestris or Snograssella alvi DNA ###
######################################################################################

# Run canu assembler on Bombus and Snodgrassella-filtered reads 
bash code/3.run_canu_assembler.sh "${isolate}" "3.no_bombus_or_sella" '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_Snodgrasella_alvi/Vairimorpha_bombi_8.3_Snodgrasella_unaligned_reads.fasta'

# Save path to this assembly to variable 
path_to_no_sella_assembly="$( echo $(pwd))/results_${isolate}/3.no_bombus_or_sella_${isolate}_assembly/3.no_bombus_or_sella_${isolate}.contigs.fasta"