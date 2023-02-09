#!/bin/bash 
## Usage: genome_assembler.sh isolate_name path_to_raw_reads path_to_lambda_reference_genome_directory path_to_bombus_genome
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

#####################################################
### Isolate microsporidia DNA from the new genome ###  
#####################################################

# Run BUSCO analysis on the newly assembled genome 
bash code/4.1.0.run_busco_analysis.sh ${isolate} $path_to_no_bombus_assembly
#bash code/4.1.0.run_busco_analysis.sh ${isolate} "/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results/canu_assemblies/Non_Bombus/"

# Convert high quality filtered fastq file to fasta file 
#sed -n '1~4s/^@/>/p;2~4p' results_Vairimorpha_bombi_8_1/1.Vairimorpha_bombi_8_1_high_qual_reads/Vairimorpha_bombi_8_1_high_qual_reads.fastq > results_Vairimorpha_bombi_8_1/1.Vairimorpha_bombi_8_1_high_qual_reads/Vairimorpha_bombi_8_1_high_qual_reads.fasta
# Run BUSCO analysis on raw reads
bash code/4.1.0.run_busco_analysis.sh ${isolate} /home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8_1/1.Vairimorpha_bombi_8_1_high_qual_reads/Vairimorpha_bombi_8_1_high_qual_reads.fasta

# Extract microsporidia-matched sequences to fasta file 
bash code/4.2.0.isolate_microsporidia_sequences.sh  ${isolate} $path_to_no_bombus_assembly
bash code/4.2.0.isolate_microsporidia_sequences.sh  ${isolate} /home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8_1/1.Vairimorpha_bombi_8_1_high_qual_reads/Vairimorpha_bombi_8_1_high_qual_reads.fasta

# Check there are no other sequences that contain Vairimorpha bombi primer regions 
bash code/4.3.0.check_Vb_primers.sh ${isolate} $path_to_no_bombus_assembly 
