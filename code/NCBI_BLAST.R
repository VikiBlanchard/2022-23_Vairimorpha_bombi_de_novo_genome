#!/usr/bin/Rscript
#' @title NCBI BLAST searches
#' @description Process BLAST results for my genome
#' @param 
#' @author viki.webster.2018@live.rhul.ac.uk

##########################
### Set up environment ###
##########################
# Clear environment 
rm(list=ls())

# Set working directory
library ("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))

# Load packages
library("ggplot2")
#install.packages("tidyverse")
library("tidyverse")
library(tidyr)
library(RColorBrewer)
#install.packages("viridis")
library("viridis")
#install.packages("readr")
library("readr")
#install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')
library ("rBLAST")
#install.packages("dplyr")
library ("dplyr")

#########################################
### Make fasta file into a data frame ###
#########################################

# Load in fasta file
V_bombi <- readDNAStringSet('/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/3.No_Ralstonia_Vairimorpha_bombi_8.1-3_assembly/3.No_Ralstonia_Vairimorpha_bombi_8.1-3.contigs.fasta')
# Store sequence names 
seq_names <- names(V_bombi)
# Store sequences
sequences <- paste(V_bombi)
# Generate data frame 
V_bombi_df <- data.frame(seq_names, sequences)

###########################################
### Manually BLAST individual sequences ###
###########################################

# Copy and paste each sequence into the NCBI webapp
V_bombi_df[47,1]
Single_Sequence <- paste( paste( '>', as.character(V_bombi_df[1,1])), as.character(V_bombi_df[1,2]))
clipr::write_clip(Single_Sequence)    

# Or save each sequence to its own fasta file to blast
for (sequence_number in 1:nrow(V_bombi_df)){
  writeLines(c(paste( '>', as.character(V_bombi_df[sequence_number,1]), sep =""),as.character(V_bombi_df[sequence_number,2])), paste("Sequence_", sequence_number, ".txt", sep = ""))}


#############################
### Read in BLAST results ###
#############################

# Change into BLAST result directory 
setwd("20221116_BLAST_Results")

# Read in the first line from every result file
  # store .csv filenames as row names
Overall_BLAST_hits <- do.call( rbind, sapply(
  list.files(".","*.csv"), read.csv, nrows=1, header = TRUE, simplify = FALSE, stringsAsFactors = TRUE
))

Overall_BLAST_hits$names <- rownames(df)

# Change back to wd
setwd("..")

# Read in GC content results
GC_per_contig <- read_delim(
  '/home/vlb19/Documents/Coding/2022_Vairimorpha_Genome/results_Vairimorpha_bombi_8.1-3/GC_content_per_tig.tsv', 
  show_col_types = FALSE, 
  col_names = FALSE)
                          
#########################
### Summarise results ###
#########################

# Look at which species have been matched
levels(Overall_BLAST_hits$Scientific.Name)

########################################
### Find contamination by GC content ###
########################################
separate(data = V_bombi_df, col = seq_names, into = c("seq_ID", "len", "read_count","class", "suggestRepeat","suggest_bubble", "suggestCircular", "trim"), sep = "\t")

# Filter GC content by 
GC_content_over_50 <- filter(GC_per_contig, X2 >= 50)
View(GC_content_over_50)
# Filter GC content by 
GC_content_under_30 <- filter(GC_per_contig, X2 <= 30)

# Search for probable contaminant sequences
V_bombi_df %>% filter_all(any_vars(. %in% c("tig00000042")))
which(V_bombi_df == "tig00000042", arr.ind = TRUE)
df %>% filter_all(any_vars(. %in% c('G')))
V_bombi_df %>% filter_all(any_vars(. %in% c(GC_content_over_50$X1)))

df %>% filter_all(any_vars(. %in% c('value1', 'value2', ...)))
V_bombi_df[which(V_bombi_df$seq_names == "tig00000042")]

#######################
### Plot BLAST hits ### 
#######################

# Make a normal barplot
ggplot(Overall_BLAST_hits, aes(x=reorder(Scientific.Name, Per..ident ))) +
    geom_bar(aes(fill = Per..ident, group = Per..ident)) + 
    #scale_fill_hue() +
    scale_color_viridis(option = "D") + 
    theme_minimal() +
    theme(text = element_text(size = 20)) + 
    labs(fill = "Percent identity")  + 
    ylab("Number of BLAST matches") + 
    xlab("BLAST Match Species") +
    coord_flip()
  #fill = Per..ident,


########################################
### Add species onto problem contigs ###
########################################

if (row.names(Overall_BLAST_hits) == GC_per_contig[1,]){
      print("yes")
}






Overall_BLAST_hits %>% 
  mutate(Scientific.Name = as.character(Scientific.Name)) %>%
  pivot_longer(-Scientific.Name) %>% 
  ggplot(aes(x=Scientific.Name, fill= Per..ident, group = Per..ident)) + 
  geom_col()

library(tidyverse)
dat %>% 
  pivot_longer(-Clades) %>% 
  ggplot(aes(x=Clades, y=value, fill=name)) + 
  geom_col()

# Make a lollipop plot 
ggplot(Overall_BLAST_hits, aes (x=Overall_BLAST_hits$Scientific.Name)) + geom_point() + geom_segment(x=Scientific.Name, xend=x, y=0 yend=y)

# Make a lollipop plot
ggplot(Overall_BLAST_hits, aes(x=Scientific.Name, y=Scientific.Name$Freq)) +
  geom_segment( aes(x=Scientific.Name, xend=Scientific.Name, y=0, yend=50), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

# Plot
ggplot(Overall_BLAST_hits, aes(x=Scientific.Name, y=c(1:180))) +
  geom_point() + 
  geom_segment( aes(x=Scientific.name, xend=x, y=0, yend=y))

       #coord_flip() + scale_color_brewer(palette = "PuOr") + theme(text = element_text(size = 20)) + labs(fill = "Percent identity")  + ylab("Number of BLAST matches") + xlab("BLAST Match Species")





####################
### Sandbox area ###
####################

# Install and load packages to use BLAST
#install.packages("Biostrings")
library ("Biostrings")
#install.packages("dplyr")
library(devtools)
library(githubinstall)
#install_github("gschofl/blastr")
#install.packages('blastr', repos = 'https://github.com/gschofl/blastr')
#library ("blastr")


# List available databases 
biomartr::listNCBIDatabases(db = "all")
# show all NCBI nr files
biomartr::listNCBIDatabases(db = "nr")
# show all NCBI nt files
biomartr::listNCBIDatabases(db = "nt")

# Make directory for nucleotide database
dir.create("nt")
# download the entire NCBI nt database
biomartr::download.database.all(db = "nt", path = "nt")


# Read in all .csv blast result files
for (tig in V_bombi_df[,1]){
  if (any(grepl(substr(tig, 1, 11),list.files())) == TRUE){
    do.call("<-", list(paste(substr(tig, 1, 11),"_blast_hits",sep = ""),
    read.csv(paste(substr(tig, 1, 11),".csv", sep = ""))))
  }}

################################################################
# Load metablastr package 
library ("metablastr")

# Test metablastr is working 
blast_test <- blast_nucleotide_to_nucleotide(
  query   = system.file("~/porechop_FAS94198_pass_756abde2_11.fastq.gz", package = 'metablastr'),
  subject = system.file('..Vairimorpha_bombi_8.1-3_Snodgrasella_unaligned_reads.fasta', package = 'metablastr'),
  output.path = tempdir(),
  db.import  = FALSE)

View(blast_test)


##########################
### BLAST my sequences ###
##########################

blastn –db nt –query ../'Vairimorpha_bombi_8.1-3_Snodgrasella_unaligned_reads.fasta' –out BLAST_results_all_ON_Reads.out  -remote

#blastn -query "test-Vb.contigs.fasta" -db nt -remote -task blastn-short -word_size 7 -evalue 500 -perc_identity 95 -entrez_query "Fungi [organism]" -outfmt 6 -out blast_result_Hsapiens.table -max_target_seqs 10 -max_hsps 5
#blastn –db nt –query "test-Vb.contigs.fasta" –out results.out -remote
#blastn –db nt –query test-Vb.contigs.fasta –out results.out
#zcat  | awk '{if (/^>/) { print ">" $2} else { print $_}}' > nt.fa
#makeblastdb(file.path(./, "nt.fa"), dbtype = "nucl")
#NucleotideDatabase <- blast (.)


# Convert my FASTQ raw read files into FASTA files for BLASTing 
#install.packages("remotes")
#remotes::install_github("robertdouglasmorrison/DuffyTools")
library(DuffyTools)
fastqToFasta("~/porechop_FAS94198_pass_756abde2_11.fastq.gz", "porechop_FAS94198_pass_756abde2_11.fasta", Qscores=TRUE)

# BLAST sequences
setwd("E:/reference_sequences/nt")
blast(db = "E:/reference_sequences/nt", type = "blastn")
blast_help(type = "blastn")

seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
                                    package="rBLAST"))

bl <- blast(db="./16S_ribosomal_RNA/16S_ribosomal_RNA")
bl
cl <- predict(bl, seq[1,])
cl[1:5,]



# Blast my sequences
V_bombi_blast <- blast_nucleotide_to_nt_database(
  query   = system.file('seqs/qry_nn.fa', package = 'metablastr'),
  output.path = tempdir(),
  db.import  = FALSE)



system2(command = blastn, 
        args = c("-db", blast_db, 
                 "-query", input, 
                 "-outfmt", format, 
                 "-evalue", evalue, 
                 "-ungapped"))

blast_out <- system2(command = blastn, 
                     args = c("-db", blast_db, 
                              "-query", input, 
                              "-outfmt", format, 
                              "-evalue", evalue,
                              "-ungapped"),
                     wait = TRUE,
                     stdout = TRUE)

