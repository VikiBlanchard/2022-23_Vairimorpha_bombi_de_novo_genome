#!/usr/bin/env Rscript
#' @title annotation_parser.R
#' @description Code for parsing .annot files
#' @param 
#' @author Viki Webster

##########################
### Set up environment ###
##########################

### Clear global environment 
rm(list = ls())

### Set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))

### Load packages 
require("readr")

###################
### Import data ###
###################

# Import annotation file 
annotation_file <- read_tsv("/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/14.Analysis/Enrichment_analyses/GO_terms/20230324_annotations_with_enzyme_names.annot", col_names = c("gene_name", "GO_term", "Description"))

# Import list of genes 
genes_of_interest <- read_tsv("/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/results_Vairimorpha_bombi_8.1-3/14.Analysis/Enrichment_analyses/SC_orthologs/V.bombi_SC_orthologs_microsporidia", col_names=FALSE)

# Get only results from genes of interest 
annotated_genes_of_interest <- annotation_file[ annotation_file$gene_name %in% genes_of_interest$gene_name, ]

# Export annotated genes of interest 
write_tsv(annotated_genes_of_interest, "SC_genes.annot")

