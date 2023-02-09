---
title: "README.md"
author: "Victoria Blanchard"
date: "11/08/2022"
output: pdf_document
---

Data should be organised into subfolders named after the isolate

## Dependencies 
Before you run the code, please ensure you have these installed and added to your path:

   - anaconda3 
   - kraken2 and Kraken2-build
      - OpenMP
      - wget
      - rsync
   - phyloseq  

## Pipeline and associated files 

To just run through the entire pipeline from the downloaded Git repository: 
   **0.V.bombi_genome_assembler.sh**
      - **Inputs:** Isolate name, Path to folder containing the raw reads for that isolate
      - **Outputs:** Fasta file with assembled genome, quality reports, all intermediate files generated in the pipeline 

Otherwise, the code can be run file by file following the pipeline:
   
1. **Filter high quality reads**
   
   1.1.initial_set-up.sh
      Assuming the Git repository is cloned and immediately run, set up directory structure and set isolate name
      - **Input:** isolate name
      - **Output:** report and results directories for the isolate

   
   1.2.filter_high_quality_reads.sh
      Remove Oxford Nanopore adapters, reads from the lambda phage control genome added in the Ligation Library Preparation, and all reads with a q score lower than 10. 
      - **Input:** 
      - **Outputs:** One fastq file containing all reads, summary plots of the raw reads in the results directory

2. **Remove contaminating *Bombus terrestris* DNA sequences**
   Align the reads to the genome of the largest known contaminant - *Bombus terrestris*, then extract unaligned reads from the resulting SAM file into FASTA and FASTQ files. 

   2.1.0.FASTQ_to_BAM_using_BWA.sh
      - **Inputs:** Reference sequence, Reads in FASTQ, Reads in paired FASTQ or "unpaired"
      - **Output:** alignment results to SAM file
      - **Uses subscript:** 2.1.1.SAM_how_many_reads_align.pl

   2.2.SAM_unaligned_reads_to_fasta.pl
      - **Input:** Sam file from 2.1
      - **Output:** FASTA file with unaligned reads
   
   2.3.SAM_unaligned_reads_to_fastq.pl
      - **Input:** Sam file from 2.1
      - **Output:** FASTA file with unaligned reads

3. **Generate an assembly from reads without *B. terrestris* DNA**
   Run Canu to assemble the genome on the remaining reads. This will contain some contamination that will be filtered further later. Quality assess the assembly using quast.

   3.run_canu_assembler.sh
      - **Input:** 
      - **Output:**

4. **Isolate microsporidia DNA from the new assembly**
    and generate a FASTA file for matched and unmatched reads respectively. Search the unmatched reads for *V. bombi* primer sequences. 
   
   4.1.0.run_busco_analysis.sh
      Search newly assembled genome for microsporidia BUSCOs
      - **Input:** 
      - **Outputs:** 
         - summary plot code: 4.1.1.busco_figure_isolate.R
         - 

   4.2.0.isolate_microsporidia_sequences.sh
      Extract microsporidia sequences from the first assembly into a new fasta file
      - **Input:** 
      - **Output:**
      - **Uses subscript:** 4.2.1.get_BUSCO_matched_sequences.py

   4.3.0.check_Vb_primers.sh
      Ensure that all sequences from the first assembly that match *V. bombi*-specific primers are included in the microsporidia-specific sequences FASTA file.      
      - **Input:** 
      - **Output:**
      - **Uses subscript:** 4.3.1.FASTA_search_Vb_primers.sh
   
5. **Generate an assembly from microsporidia-specific reads**
   Generate a new assembly and thoroughly quality test it. 

   3.run_canu_assembler.sh
       Re-run canu script with this new input FASTA file to assemble a genome from reads we know are from the microsporidia, and therefore presumably only *V. bombi* DNA
      - **Input:** 
      - **Output:**

   5.1.genome_quality_assessment.sh
      - **Input:** 
      - **Output:**

6. **Annotate assembly**

   6.run_braker2.sh
      - **Input:** 
      - **Output:**
