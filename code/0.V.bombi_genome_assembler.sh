#!/bin/bash 
## Usage: genome_assembler.sh isolate_name path_to_raw_reads path_to_lambda_reference_genome_directory path_to_bombus_genome
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

# bash FASTQ to FASTA and back converter 
#sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta
#sed -n '1~4s/^>/@/p;2~4p' INFILE.fasta > OUTFILE.fastq

##########################
### Set-up environment ###
##########################

# Store isolate name
isolate="Vairimorpha_bombi_8.1-3"


# Initialise output directories
bash code/1.1.initial_set-up.sh ${isolate}

# Path to raw nanopore reads
raw_nanopore_reads='/media/vlb19/Expansion/Vairimorpha_bombi_sequences/Combined_fastq_passes'

# Path to Lambda reference genome directory
Lambda_Genome_Directory="$( echo $(pwd))/data"

# Path to reference genomes
B_terrestris_genome='/media/vlb19/Expansion/reference_sequences/Bombus/Bombus_terrestris/GCF_910591885.1iyBomTerr1.2/GCF_910591885.1_iyBomTerr1.2.genome.fa'

#################################
### Filter high quality reads ###
#################################

# Trim adapters, remove sequences from the reference genome, remove contaminant barcode DNA 
bash code/1.2.filter_high_quality_reads.sh $isolate $raw_nanopore_reads $Lambda_Genome_Directory


###############################
### BLAST genomic sequences ###
###############################

/home/vlb19/Documents/Coding/Downloaded_Repositories/ncbi-blast-2.13.0+-x64-linux/ncbi-blast-2.13.0+/bin/blastn -db nt -query '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/2.align_reads_to_sella_alvi/Vairimorpha_bombi_8.3_sella_unaligned_reads.fasta' -out results_Vairimorpha_bombi_8.1-3/BLAST_results.out -remote

export PATH=$PATH:/home/vlb19/Documents/Coding/Downloaded_Repositories/ncbi-blast-2.13.0+-x64-linux/ncbi-blast-2.13.0+/bin
export BLASTDB=$HOME/blastdb

perl /home/vlb19/Documents/Coding/Downloaded_Repositories/ncbi-blast-2.13.0+-x64-linux/ncbi-blast-2.13.0+/bin/update_blastdb.pl --source ncbi --decompress


##########################################
### Remove contaminating DNA sequences ###
##########################################
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


#######################
### Assemble genome ### 
#######################

# Run canu assembler on reads that did not align to Bombus terrestris
bash code/3.run_canu_assembler.sh "${isolate}" "3.No_B_terrestris" "../2.align_reads_to_Bombus_terrestris/${isolate}_unaligned_high_qual_reads.fasta" 

# Run canu assembler on raw high quality reads
#bash code/3.run_canu_assembler.sh "${isolate}" "3.raw" "../1.${isolate}_high_qual_reads/${isolate}_high_qual_reads.fastq" 

# Look at GC content of every sequence in FASTA file 
seqkit fx2tab --name --only-id --gc "results_${isolate}/3.raw_${isolate}_assembly/3.raw_$isolate.contigs.fasta"

seqkit fx2tab --name --only-id --gc 'results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta' | cat > "results_Vairimorpha_bombi_8.1-3/${isolate}_content_per_tig.tsv"

####################################
### Re-order and re-name contigs ###
####################################

# Re-order contigs 
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s '3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta' -j b -p fasta > '3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered.fasta'

# Make list of contigs
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s  '3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered.fasta' -h s -g n | cut -f 1 > 'contigs.tab' 

# Add new contig ID prefix 
sed -e 's/$/\tV_bombi_8.1-3_tig0/' -i 'contigs.tab' 

# Add line number new contig ID prefix 
perl '../../code/quicklinenumber.pl' 'contigs.tab' > contigs-new.tab

# Rename contigs
perl ../../code/FASTA-rename-according-to-file.pl 3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered.fasta contigs-new.tab  > 3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered-renamed.fasta


###############################################
### Prep short read sequences for polishing ###
###############################################
sam_short_reads='/media/vlb19/Expansion1/Vairimorpha_bombi_sequences/Short_read_sequences/Sam_short_read_sequences'
novogene_short_reads='/media/vlb19/Expansion1/Vairimorpha_bombi_sequences/Short_read_sequences/Novogene_seqs'
no_bombus_assembly='/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered-renamed.fasta'

# Combine all forward short reads from Sam into one file
zcat ${sam_short_reads}/*/*_R1.fastq.gz | gzip >> ${sam_short_reads}/Sam_combined_short_reads_R1.fastq.gz
# Combine all reverse short reads from Sam into one file
zcat ${sam_short_reads}/*/*_R2.fastq.gz | gzip >> ${sam_short_reads}/Sam_combined_short_reads_R2.fastq.gz

# Combine all novogene forward short reads into one file
zcat ${novogene_short_reads}/*/*_1.fq.gz | gzip >> ${novogene_short_reads}/Novogene_combined_short_reads_R1.fastq.gz
# Combine all novogene reverse short reads into one file
zcat ${novogene_short_reads}/*/*_2.fq.gz | gzip >> ${novogene_short_reads}/Novogene_combined_short_reads_R2.fastq.gz

# Trim adapters from short reads from Sam
cd '/home/vlb19/Documents/Coding/Downloaded_Repositories/Trimmomatic/dist/jar'

trimlog_file="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/sam_reads_trimmomatic_log.txt"
summary_file="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/sam_reads_stats_summary.txt"

short_read_f1=${sam_short_reads}/Sam_combined_short_reads_R1.fastq.gz
short_read_f2=${sam_short_reads}/Sam_combined_short_reads_R2.fastq.gz

trimmed_short_read_f1="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Sam_combined_short_reads_R1_trimmed.fastq.gz"
trimmed_short_read_f2="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Sam_combined_short_reads_R2_trimmed.fastq.gz"

untrimmed_short_read_f1="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Sam_combined_short_reads_R1_untrimmed.fastq.gz"
untrimmed_short_read_f2="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Sam_combined_short_reads_R2_untrimmed.fastq.gz"

java -jar trimmomatic-0.40-rc1.jar PE -threads 4 -phred33 -trimlog $trimlog_file -summary $summary_file $short_read_f1 $short_read_f2 $trimmed_short_read_f1 $untrimmed_short_read_f1 $trimmed_short_read_f2 $untrimmed_short_read_f2 ILLUMINACLIP:../../adapters/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36 


# Trim adapters from novogene files
trimlog_file="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/novogene_reads_trimmomatic_log.txt"
summary_file="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/novogene_reads_stats_summary.txt"

short_read_f1=${novogene_short_reads}/Novogene_combined_short_reads_R1.fastq.gz
short_read_f2=${novogene_short_reads}/Novogene_combined_short_reads_R2.fastq.gz

trimmed_short_read_f1="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Novogene_combined_short_reads_R1_trimmed.fastq.gz"
trimmed_short_read_f2="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Novogene_combined_short_reads_R2_trimmed.fastq.gz"

untrimmed_short_read_f1="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Novogene_combined_short_reads_R1_untrimmed.fastq.gz"
untrimmed_short_read_f2="/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Novogene_combined_short_reads_R2_untrimmed.fastq.gz"

java -jar trimmomatic-0.40-rc1.jar PE -threads 4 -phred33 -trimlog $trimlog_file -summary $summary_file $short_read_f1 $short_read_f2 $trimmed_short_read_f1 $untrimmed_short_read_f1 $trimmed_short_read_f2 $untrimmed_short_read_f2 ILLUMINACLIP:../../adapters/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36 

# Change back to results directory 
cd '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep'

# Unzip compressed fastq files 
gunzip ./*/*.fastq.gz

# Align all reads from Sam to the canu assembly
bash ../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${no_bombus_assembly}" "/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Sam_short_read_sequences/Sam_combined_short_reads_R1_trimmed.fastq" "/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Sam_short_read_sequences/Sam_combined_short_reads_R2_trimmed.fastq" 

# Convert SAM file to fastq file 
perl ../../code/2.3.SAM_unaligned_reads_to_fastq.pl "Sam_short_read_sequences/Sam_combined_short_reads_R1_trimmed.fastq-mem.sam" > "Unaligned_Sam_short_reads.fastq"

# Align all Novogene reads to the canu assembly
bash ../../code/2.1.0.FASTQ_to_BAM_using_BWA.sh "${no_bombus_assembly}" "/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Novogene_seqs/Novogene_combined_short_reads_R1_trimmed.fastq" "/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/Short_read_prep/Novogene_seqs/Novogene_combined_short_reads_R2_trimmed.fastq" 

# Convert SAM file to fastq file 
perl ../../code/2.3.SAM_aligned_reads_to_fastq.pl "Novogene_seqs/Novogene_combined_short_reads_R1_trimmed.fastq-mem.sam" > "Aligned_Novogene_short_reads.fastq"

# Combine all reads
cat ./*.fastq > all_short_reads.fastq

# Normalise read coverage to 120x
~/Documents/Coding/Downloaded_Repositories/bbmap/bbnorm.sh in=all_short_reads.fastq out=all_short_reads_normalised.fastq target=120 min=5


###############################################
### Polish genome with short read sequences ###
###############################################
cd "/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3"

# Change into appropriate directory 
cd Novogene_reads_filtered 

# Copy filtered file over 
cp ../../Short_read_prep/Aligned_Novogene_short_reads.fastq ./Aligned_Novogene_short_reads.fastq

# Run 3 rounds of polishing on the genome
bash ../../../code/4.0_polish_genome_with_pilon.sh ${no_bombus_assembly} 1 ../Short_read_prep/all_short_reads.fastq # round 1
bash ../../../code/4.0_polish_genome_with_pilon.sh Polish_loop_1.fasta.fasta 2 ../Short_read_prep/all_short_reads.fastq # round 2
bash ../../../code/4.0_polish_genome_with_pilon.sh Polish_loop_2.fasta.fasta 3 ../Short_read_prep/all_short_reads.fastq # round 3

# Change into appropriate directory 
cd cd ../Sam_reads_filtered/

# Copy filtered file over 
cp ../../Short_read_prep/Aligned_Sam_short_reads.fastq ./Aligned_Sam_short_reads.fastq

# Run 3 rounds of polishing on the genome
bash ../../../code/4.0_polish_genome_with_pilon.sh ${no_bombus_assembly} 1 Aligned_Sam_short_reads.fastq # round 1
bash ../../../code/4.0_polish_genome_with_pilon.sh Polish_loop_1.fasta.fasta 2 Aligned_Sam_short_reads.fastq # round 2
bash ../../../code/4.0_polish_genome_with_pilon.sh Polish_loop_1.fasta.fasta 3 Aligned_Sam_short_reads.fastq # round 3


###################################################
### Filter contigs with read depth less than 1x ###
################################################### 

# mkdir "results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage"

# Generate file with depth per read
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/SAM_how_many_reads_align.pl' -s '../Short_read_prep/all_short_reads.fastq-mem.sam'

# Generate file with tig names with less than 1x coverage 
perl "../../../code/read_coverage_with_contig_names.pl" "all_short_reads.fastq-mem.sam-reads-aligned-to-which-contigs3.tab" 1 > "contigs_with_greater_than_1x_coverage.txt"

# Filter reads
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s 'Polish_loop_3.fasta.fasta' -l 'contigs_with_greater_than_1x_coverage.txt' -p fasta > 'contigs_with_more_than_1x_coverage.fasta'

# Reorder contigs 
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s contigs_with_more_than_1x_coverage.fasta -j b -p fasta > contigs_with_more_than_1x_coverage.fasta-reordered.fasta 


##############################################################################
### Filter only contigs with matched sequences in V. apis and/or V.ceranae ###
##############################################################################

# Get list of contigs with matches in V. apis and/or V. ceranae 
sed -nr 's/.*Vairimorpha_bombi_8.1-3no1x\t(.*)/\1/p' '/home/vlb19/Documents/Coding/2021-2022_Microsporidia_Comparative_Genomics/synima_runs/Synima_subset/Vairimorpha_bombi_8.1-3no1x_vs_apis_and_ceranae/Repo_spec.txt.dagchainer.aligncoords' > matched_contig_IDs.txt
sed -i "s/Vairimorpha_bombi_8.1-3no1x.*/ /g" matched_contig_IDs.txt

# Get unique contig IDs with matches
cat matched_contig_IDs.txt  | sort | uniq > unique_matched_contig_IDs.txt

# Get list of contigs in assembly 
sed -nr 's/>(.*)/\1/p' '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.SynIma_files_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/Vairimorpha_bombi_8.1-3no1x.genome.fa' | sort > all_contigs.txt

# Get list of contigs with no matches 
join -v1 -v2 unique_matched_contig_IDs.txt all_contigs.txt > contigs_with_no_matches.txt

# Remove contigs with no matches 
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.SynIma_files_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/Vairimorpha_bombi_8.1-3no1x.genome.fa' -l 'contigs_with_no_matches.txt' -p fasta > 'Ortholog-match_filtered_assembly.fasta' 

###########################
### Softmask the genome ### 
###########################

# download a preformatted NCBI BLAST database
#perl update_blastdb.pl --decompress nt [*]

### Identify novel repeats in genome assembly
# format NCBI nt database for RepeatModeler
#BuildDatabase -name nt_database_RepeatModeler -engine ncbi -dir '/media/vlb19/Expansion/reference_sequences/NCBI_nt_database' 

# Run RepeatModeler
#RepeatModeler -database nt_database_RepeatModeler '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/7.Vairimorpha_bombi_8.1-3_Pilon_Polish/Pilon_polished_assembly.fa.fasta' -noisy -pa 10 -LTRStruct -poly -html
# -database: prefix name of the database that is used in the BuildDatabase function
# -pa: number of threads
# -LTRStruct: runs the LTR structural discovery pipeline for discovering LTR retrotransposons
# -noisy: verbose option 
# -poly: reports repeats that may be polymorphic into file.poly 
# -html: makes additional output file in html format 

# Repeatmask the genome with fungal sequences
#RepeatMasker -species Nosema_ceranae -pa 8 -dir results_${isolate}/RepeatMask -gff -e ncbi -s '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/3.No_Ralstonia_Vairimorpha_bombi_8.1-3_assembly/3.No_Ralstonia_Vairimorpha_bombi_8.1-3.contigs.fasta'

# Repeat masking with RepeatMasker using repeats were found by RepeatModeler
#RepeatMasker -lib CLR-families.fa -s -parallel 10 -xsmall -alignments CLR_scaffolst

# Repeatmask the genome with fungal, viral, and bacterial sequences
RepeatMasker -species fungi -pa 8 -dir ./ -gff -e ncbi -s 'Polish_loop_3.fasta.fasta'
RepeatMasker -species bacteria -pa 8 -dir ./ -gff -e ncbi -s 'Polish_loop_3.fasta.fasta.masked'
RepeatMasker -species viruses -pa 8 -dir ./ -gff -e ncbi -s 'Polish_loop_3.fasta.fasta.masked.masked'


#######################
### Annotate genome ###
#######################

# Run Braker2 on masked genome
braker.pl --genome="Polish_loop_3.fasta.fasta.masked.masked" --esmode --softmasking --cores 4 --AUGUSTUS_BIN_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/bin --AUGUSTUS_SCRIPTS_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts

# additional option "--fungus" is available


######################################
### Generate protein and CDS files ###
######################################
mkdir results_Vairimorpha_bombi_8.1-3/9.SynIma_files_Vairimorpha_bombi_8.1-3
cd "results_Vairimorpha_bombi_8.1-3/9.SynIma_files_Vairimorpha_bombi_8.1-3"

# Generate gff3 from gtf
cat 'braker.gtf' | /home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts/gtf2gff.pl --gff3 --out=braker.gff3 
# Remove braker sequence prediction notes
sed -i 's/_len=.*=no//g' braker.gff3 

# Remove Augustus and GeneMark additions
sed -i 's/AUGUSTUS/Vairimorpha_bombi_8_1/g' braker.gff3
sed -i 's/GeneMark.hmm3/Vairimorpha_bombi_8_1/g' braker.gff3

# Remove start and stop codons
awk -i inplace '!/start_codon/' braker.gff3
awk -i inplace '!/stop_codon/' braker.gff3

# Generate .CDS file using gffread  
gffread -g 'Polish_loop_3.fasta.fasta.masked.masked' -x ${isolate}.annotation.cds braker.gff3

# Generate protein.fa file using gffread  
gffread -g 'Polish_loop_3.fasta.fasta.masked.masked' -y ${isolate}.annotation.pep braker.gff3 

# Rename files for Synima
mv Polish_loop_3.fasta.fasta.masked.masked ${isolate}.genome.fa
mv braker.gff3 ${isolate}.annotation.gff3 # gff3


##################
### Run SynIma ###
##################

# Move into working directory

# Initiate repo_spec document 
printf "// \n" > Repo_spec.txt

# Append all species names in to repo_spec with correct requirements
for dir in ./*; do 
    if [ "$(basename -- "$dir")" != "Repo_spec.txt" ]; then
        printf "Genome $(basename -- "$dir")\nAnnotation $(basename -- "$dir")\n//\n" >> Repo_spec.txt
    fi
done    

# Run SynIma from data directory
perl ../util/Create_full_repo_sequence_databases.pl -f mRNA -r ./Repo_spec.txt -v
perl ../util/Blast_grid_all_vs_all.pl -r ./Repo_spec.txt
perl ../util/Blast_all_vs_all_repo_to_OrthoMCL.pl -r ./Repo_spec.txt
#ALTERNATIVELY 1: ../util/Blast_all_vs_all_repo_to_RBH.pl -r ./Repo_spec.txt
#ALTERNATIVELY 2: ../util/Blast_all_vs_all_repo_to_Orthofinder.pl -r ./Repo_spec.txt
perl ../util/Orthologs_to_summary.pl -o all_orthomcl.out
perl ../util/DAGchainer_from_gene_clusters.pl -r ./Repo_spec.txt \
             -c GENE_CLUSTERS_SUMMARIES.OMCL/GENE_CLUSTERS_SUMMARIES.clusters
perl ../SynIma.pl -a Repo_spec.txt.dagchainer.aligncoords -b Repo_spec.txt.dagchainer.aligncoords.spans


##############################################
### Run BUSCO analysis on the new assembly ###
##############################################

bash code/4.1.0.run_busco_analysis.sh ${isolate} '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3/Pilon.1.2.4/Pilon_polished_assembly.fasta-reordered-renamed.fasta' "No_bombus"

