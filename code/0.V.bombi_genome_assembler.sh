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

Snodgrassella_genome='/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/data/Snodgrassella_alvi_wkB2/GCF_000600005.1_ASM60000v1_genomic.fa'
Ralstonia_sp_genome='/media/vlb19/Expansion/reference_sequences/Ralstonia_sp_56D2/SRR11285236/SRR11285236.fasta'
Cupriavidus_necator_genome='/media/vlb19/Expansion/reference_sequences/Cupriavidus_necator_H16/GCF_004798725.1_ASM479872v1/GCF_004798725.1_ASM479872v1_genome.fa'

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
### Polish genome with short read sequences ###
###############################################
short_reads='/media/vlb19/Expansion1/Vairimorpha_bombi_sequences/Short_read_sequences'

# Polish genome 
bash code/4.0_polish_genome_with_pilon.sh ${isolate} "results_${isolate}/3.No_B_terrestris_${isolate}/3.No_B_terrestris_${isolate}.contigs.fasta" ${short_reads}

### Generate window plot 
# Generate pileup
samtools mpileup -f  ../3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered-renamed.fasta -s "V_bombi_combined_short_reads.fastq-mem.sorted.bam" -o "V_bombi_window_plot.pileup"

# Generate windows dataframe
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_for_VCFs_mpileups_or_tabs3.pl' -r '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta-reordered.fasta' -v 'Pilon_polished_assembly.vcf' -p "V_bombi_window_plot.pileup" -w 1000 -n 166.25 > 1000-windows_dataframe-normalised.tab

perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_for_VCFs_mpileups_or_tabs3.pl' -r '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/3.No_B_terrestris_Vairimorpha_bombi_8.1-3_assembly/3.No_B_terrestris_Vairimorpha_bombi_8.1-3.contigs.fasta' -v 'Pilon_polished_assembly.vcf' -p "V_bombi_window_plot.pileup" -w 1000 > 1000-windows_dataframe_not_normalised.tab

# Make windows plot
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_dataframe_to_R_figure3.pl' -w 1000-windows_dataframe.tab

perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/Windows_dataframe_to_R_figure3.pl' -w 1000-windows_dataframe_not_normalised.tab

###################################################
### Filter contigs with read depth less than 1x ###
################################################### 

mkdir "results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage"

# Generate file with depth per read
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/SAM_how_many_reads_align.pl' -s '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3/V_bombi_combined_short_reads.fastq-mem.sam'

# Generate file with tig names with less than 1x coverage 
perl "code/read_coverage_with_contig_names.pl" "results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3/V_bombi_combined_short_reads.fastq-mem.sam-reads-aligned-to-which-contigs3.tab" 1 > "results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_greater_than_1x_coverage.txt"

# Filter reads
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s 'results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3/Pilon_polished_assembly.fasta-reordered.fasta' -l 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_greater_than_1x_coverage.txt' -p fasta > 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_more_than_1x_coverage.fasta'

# Reorder contigs 
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_more_than_1x_coverage.fasta -j b -p fasta > results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_more_than_1x_coverage.fasta-reordered.fasta 

# Add line number to contigs 
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_more_than_1x_coverage.fasta-reordered.fasta  -h s -g n | cut -f 1 > 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs.tab' 

sed -e 's/$/\tV_bombi_8.1-3_tig0/' -i 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs.tab'  

perl 'code/quicklinenumber.pl' 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs.tab'   > 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs-new.tab' 

# Rename contigs
perl 'code/FASTA-rename-according-to-file.pl' results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_more_than_1x_coverage.fasta-reordered.fasta results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs-new.tab > 'results_Vairimorpha_bombi_8.1-3/8.Tigs_with_more_than_1x_coverage/contigs_with_more_than_1x_coverage.fasta-reordered-renamed.fasta'

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

mkdir results_${isolate}/5.RepeatMasked_${isolate}
mkdir results_${isolate}/9.RepeatMasked_${isolate}_1x_coverage_filtered
mkdir results_${isolate}/13.RepeatMasked_${isolate}_ortholog_matches_only
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
cd 5.RepeatMasked_${isolate}
RepeatMasker -species fungi -pa 8 -dir ./ -gff -e ncbi -s 'Pilon_polished_assembly.fasta-reordered-renamed.fasta'
RepeatMasker -species bacteria -pa 8 -dir ./ -gff -e ncbi -s 'Pilon_polished_assembly.fasta-reordered-renamed.fasta.masked'
RepeatMasker -species viruses -pa 8 -dir ./ -gff -e ncbi -s 'Pilon_polished_assembly.fasta-reordered-renamed.fasta.masked.masked'

cd 9.RepeatMasked_${isolate}_1x_coverage_filtered
RepeatMasker -species fungi -pa 8 -dir ./ -gff -e ncbi -s 'contigs_with_more_than_1x_coverage.fasta-reordered-renamed.fasta'
RepeatMasker -species bacteria -pa 8 -dir ./ -gff -e ncbi -s 'contigs_with_more_than_1x_coverage.fasta-reordered-renamed.fasta.masked'
RepeatMasker -species viruses -pa 8 -dir ./ -gff -e ncbi -s 'contigs_with_more_than_1x_coverage.fasta-reordered-renamed.fasta.masked.masked'

cd 13.RepeatMasked_${isolate}_ortholog_matches_only 
RepeatMasker -species fungi -pa 8 -dir ./ -gff -e ncbi -s 'Ortholog-match_filtered_assembly.fasta'
RepeatMasker -species bacteria -pa 8 -dir ./ -gff -e ncbi -s 'Ortholog-match_filtered_assembly.fasta.masked'
RepeatMasker -species viruses -pa 8 -dir ./ -gff -e ncbi -s 'Ortholog-match_filtered_assembly.fasta.masked.masked'

cd 16.Polish_step_2
RepeatMasker -species fungi -pa 8 -dir ./ -gff -e ncbi -s 'Pilon_polished_assembly2.fasta-reordered-renamed.fasta'
RepeatMasker -species bacteria -pa 8 -dir ./ -gff -e ncbi -s 'Pilon_polished_assembly2.fasta-reordered-renamed.fasta.masked'
RepeatMasker -species viruses -pa 8 -dir ./ -gff -e ncbi -s 'Pilon_polished_assembly2.fasta-reordered-renamed.fasta.masked.masked'

#######################
### Annotate genome ###
#######################

mkdir results_${isolate}/6.Braker2_annotation_${isolate}
cp  "results_${isolate}/5.RepeatMasked_${isolate}/Pilon_polished_assembly.fasta-reordered-renamed.fasta.masked" "results_${isolate}/6.Braker2_annotation_${isolate}/${isolate}_masked.fasta.masked"
cd results_${isolate}/6.Braker2_annotation_${isolate}

# Run Braker2 on masked genome
braker.pl --genome="Pilon_polished_assembly.fasta-reordered-renamed.fasta.masked.masked" --esmode --softmasking --cores 4 --AUGUSTUS_BIN_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/bin --AUGUSTUS_SCRIPTS_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts

braker.pl --genome="Pilon_polished_assembly2.fasta-reordered-renamed.fasta.masked.masked" --esmode --softmasking --cores 4 --AUGUSTUS_BIN_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/bin --AUGUSTUS_SCRIPTS_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts

mkdir results_${isolate}/10.Braker2_annotation_${isolate}_1x_coverage_filtered
cp 'results_Vairimorpha_bombi_8.1-3/9.RepeatMasked_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/contigs_with_more_than_1x_coverage.fasta-reordered-renamed.fasta.masked' "results_${isolate}/10.Braker2_annotation_${isolate}_1x_coverage_filtered/contigs_with_more_than_1x_coverage.fasta.masked"
cd results_${isolate}/10.Braker2_annotation_${isolate}_1x_coverage_filtered

braker.pl --genome='./contigs_with_more_than_1x_coverage.fasta.masked' --esmode --softmasking --cores 4 --AUGUSTUS_BIN_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/bin --AUGUSTUS_SCRIPTS_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts

mkdir results_${isolate}/14.Braker2_annotation_${isolate}Ortholog-match_filtered
cp 'results_Vairimorpha_bombi_8.1-3/13.RepeatMasked_Vairimorpha_bombi_8.1-3_ortholog_matches_only/Ortholog-match_filtered_assembly.fasta.masked' "results_${isolate}/14.Braker2_annotation_${isolate}Ortholog-match_filtered/Ortholog-match_filtered_assembly.fasta.masked"
cd results_${isolate}/14.Braker2_annotation_${isolate}Ortholog-match_filtered

braker.pl --genome='./Ortholog-match_filtered_assembly.fasta.masked' --esmode --softmasking --cores 4 --AUGUSTUS_BIN_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/bin --AUGUSTUS_SCRIPTS_PATH=/home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts
# additional option "--fungus" is available

# Run Braker2 on unmasked genome
#bash code/6.run_braker2.sh ${isolate} '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1_and_8.2/3.no_bombus_Vairimorpha_bombi_8_assembly/3.no_bombus_Vairimorpha_bombi_8.contigs.fasta' "no_bombus"

cd ../..

######################################
### Generate protein and CDS files ###
######################################
mkdir results_Vairimorpha_bombi_8.1-3/9.SynIma_files_Vairimorpha_bombi_8.1-3
cd "results_Vairimorpha_bombi_8.1-3/9.SynIma_files_Vairimorpha_bombi_8.1-3"

# Generate gff3 from gtf
cat '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/6.Braker2_annotation_Vairimorpha_bombi_8.1-3/braker/braker.gtf' | /home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts/gtf2gff.pl --gff3 --out=braker.gff3 
# Remove braker sequence prediction notes
sed -i 's/_len=.*=no//g' braker.gff3 

# Remove Augustus and GeneMark additions
sed -i 's/AUGUSTUS/Vairimorpha_bombi_8_1/g' braker.gff3
sed -i 's/GeneMark.hmm3/Vairimorpha_bombi_8_1/g' braker.gff3

# Remove start and stop codons
awk -i inplace '!/start_codon/' braker.gff3
awk -i inplace '!/stop_codon/' braker.gff3

# Generate .CDS file using gffread  
gffread -g '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/6.Braker2_annotation_Vairimorpha_bombi_8.1-3/Vairimorpha_bombi_8.1-3_masked.fasta.masked' -x Vairimorpha_CDS.fasta braker.gff3
#gffread -g Pilon_polished_assembly2.fasta-reordered-renamed.fasta.masked.masked -x Vairimorpha_bombi_8.1-3_polish2.annotation.cds braker.gff3

# Generate protein.fa file using gffread  
gffread -g '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/6.Braker2_annotation_Vairimorpha_bombi_8.1-3/Vairimorpha_bombi_8.1-3_masked.fasta.masked' -y Vairimorpha_proteins.fasta braker.gff3 
#gffread -g '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1_and_8.2/3.no_bombus_Vairimorpha_bombi_8_assembly/3.no_bombus_Vairimorpha_bombi_8.contigs.fasta' Vairimorpha_proteins.fasta -y '/results_Vairimorpha_bombi_8.1_and_8.2/6.annotation/no_bombus_Vairimorpha_bombi_8_assembly_augustus_annotation.gff'
#gffread -g Pilon_polished_assembly2.fasta-reordered-renamed.fasta.masked.masked -y Vairimorpha_bombi_8.1-3_polish2.annotation.pep braker.gff3

# Rename files for Synima
cp '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/6.Braker2_annotation_Vairimorpha_bombi_8.1-3/Vairimorpha_bombi_8.1-3_masked.fasta.masked' ${isolate}.genome.fa # genome
mv braker.gff3 ${isolate}.annotation.gff3 # gff3
mv Vairimorpha_CDS.fasta ${isolate}.annotation.cds # CDS
mv Vairimorpha_proteins.fasta ${isolate}.annotation.pep # pep

cd ../.. 

###################################################################
### Generate SynIma files for assembly with 1x coverage removed ###
###################################################################

mkdir results_Vairimorpha_bombi_8.1-3/10.SynIma_files_Vairimorpha_bombi_8.1-3_1x_coverage_filtered
cd "results_Vairimorpha_bombi_8.1-3/10.SynIma_files_Vairimorpha_bombi_8.1-3_1x_coverage_filtered"

# Generate gff3 from gtf
cat '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.Braker2_annotation_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/braker/braker.gtf' | /home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts/gtf2gff.pl --gff3 --out=braker.gff3 
# Remove braker sequence prediction notes
sed -i 's/_len=.*=no//g' braker.gff3 

# Remove Augustus and GeneMark additions
sed -i 's/AUGUSTUS/Vairimorpha_bombi_8_1/g' braker.gff3
sed -i 's/GeneMark.hmm3/Vairimorpha_bombi_8_1/g' braker.gff3

# Remove start and stop codons
awk -i inplace '!/start_codon/' braker.gff3
awk -i inplace '!/stop_codon/' braker.gff3

# Generate .CDS file using gffread  
gffread -g '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.Braker2_annotation_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/contigs_with_more_than_1x_coverage.fasta.masked' -x Vairimorpha_CDS.fasta braker.gff3

# Generate protein.fa file using gffread  
gffread -g '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.Braker2_annotation_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/contigs_with_more_than_1x_coverage.fasta.masked' -y Vairimorpha_proteins.fasta braker.gff3 
#gffread -g '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1_and_8.2/3.no_bombus_Vairimorpha_bombi_8_assembly/3.no_bombus_Vairimorpha_bombi_8.contigs.fasta' Vairimorpha_proteins.fasta -y '/results_Vairimorpha_bombi_8.1_and_8.2/6.annotation/no_bombus_Vairimorpha_bombi_8_assembly_augustus_annotation.gff'

# Rename files for Synima
cp '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/10.Braker2_annotation_Vairimorpha_bombi_8.1-3_1x_coverage_filtered/contigs_with_more_than_1x_coverage.fasta.masked' ${isolate}.genome.fa # genome
mv braker.gff3 ${isolate}.annotation.gff3 # gff3
mv Vairimorpha_CDS.fasta ${isolate}.annotation.cds # CDS
mv Vairimorpha_proteins.fasta ${isolate}.annotation.pep # pep


cp Pilon_polished_assembly2.fasta-reordered-renamed.fasta.masked.masked Vairimorpha_bombi_8.1-3_polish2.genome.fa
cp braker.gff3 Vairimorpha_bombi_8.1-3_polish2.annotation.gff3
cd ../..


####################################################################
### Generate SynIma files for assembly with only contigs matched ###
####################################################################

mkdir results_Vairimorpha_bombi_8.1-3/15.SynIma_files_Vairimorpha_bombi_8.1-3_ortholog_matches_only
cd "results_Vairimorpha_bombi_8.1-3/15.SynIma_files_Vairimorpha_bombi_8.1-3_ortholog_matches_only"

# Generate gff3 from gtf
cat '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/14.Braker2_annotation_Vairimorpha_bombi_8.1-3Ortholog-match_filtered/braker/braker.gtf' | /home/vlb19/Documents/Coding/Downloaded_Repositories/Augustus/scripts/gtf2gff.pl --gff3 --out=braker.gff3 
# Remove braker sequence prediction notes
sed -i 's/_len=.*=no//g' braker.gff3 

# Remove Augustus and GeneMark additions
sed -i 's/AUGUSTUS/Vairimorpha_bombi_8_1/g' braker.gff3
sed -i 's/GeneMark.hmm3/Vairimorpha_bombi_8_1/g' braker.gff3

# Remove start and stop codons
awk -i inplace '!/start_codon/' braker.gff3
awk -i inplace '!/stop_codon/' braker.gff3

# Generate .CDS file using gffread  
gffread -g '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/14.Braker2_annotation_Vairimorpha_bombi_8.1-3Ortholog-match_filtered/Ortholog-match_filtered_assembly.fasta.masked' -x Vairimorpha_CDS.fasta braker.gff3

# Generate protein.fa file using gffread  
gffread -g '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/14.Braker2_annotation_Vairimorpha_bombi_8.1-3Ortholog-match_filtered/Ortholog-match_filtered_assembly.fasta.masked' -y Vairimorpha_proteins.fasta braker.gff3 
#gffread -g '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1_and_8.2/3.no_bombus_Vairimorpha_bombi_8_assembly/3.no_bombus_Vairimorpha_bombi_8.contigs.fasta' Vairimorpha_proteins.fasta -y '/results_Vairimorpha_bombi_8.1_and_8.2/6.annotation/no_bombus_Vairimorpha_bombi_8_assembly_augustus_annotation.gff'

# Rename files for Synima
cp '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/14.Braker2_annotation_Vairimorpha_bombi_8.1-3Ortholog-match_filtered/Ortholog-match_filtered_assembly.fasta.masked' ${isolate}.genome.fa # genome
mv braker.gff3 ${isolate}.annotation.gff3 # gff3
mv Vairimorpha_CDS.fasta ${isolate}.annotation.cds # CDS
mv Vairimorpha_proteins.fasta ${isolate}.annotation.pep # pep

cd ../..
RepeatMasker -species fungi bacteria viruses -pa 8 -dir results_${isolate}/-gff -e ncbi -s "results_Vairimorpha_bombi_8.1-3/12.Ortholog-match_filtered_assembly/"
unctional_annotation_${isolate}
cd results_${isolate}/Functional_annotation_${isolate}

# Set up orthologer project 
bash '/home/vlb19/Documents/Coding/Downloaded_Repositories/orthologer_2.8.3/ORTHOLOGER-2.8.3/bin/setup.sh' empty
./setup_project_vlb19.sh

# Check the pipeline setup 
./common.sh
./sbin/check_pipeline_odb.sh

# Reformat fasta files 
'/home/vlb19/Documents/Coding/Downloaded_Repositories/orthologer_2.8.3/ORTHOLOGER-2.8.3/bin/build_import_file.sh'

# Import fasta files given by fasta.txt into existing project 
./sbin/build_import_file.sh '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/7.Vairimorpha_bombi_8.1-3_Pilon_Polish/Pilon_polished_assembly.fa.fasta'

#After checking your setup is OK, run all steps
./orthologer.sh -r all

# download OrthoDB data for Microsporidia 
./sbin/mapping.sh DOWNLOAD 6029

# create user 'V_bomi_genome'
/home/vlb19/Documents/Coding/Downloaded_Repositories/orthologer_2.8.3/ORTHOLOGER-2.8.3/bin/mapping.sh CREATE V_bomi_genome GO



##############################################
### Run BUSCO analysis on the new assembly ###
##############################################

bash code/4.1.0.run_busco_analysis.sh ${isolate} '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/4.Pilon_polished_Vairimorpha_bombi_8.1-3/Pilon.1.2.4/Pilon_polished_assembly.fasta-reordered-renamed.fasta' "No_bombus"

