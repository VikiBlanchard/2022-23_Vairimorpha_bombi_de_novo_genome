#!/bin/bash 
## Usage: annotation_file_parser.sh
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.webster.2018@live.rhul.ac.uk)

# Get list of V.bombi orthologs across microsporidia 
sed "s/^.*V_bombi/V_bombi/g" shared_sc_orthogroups > V.bombi_SC_orthologs_microsporidia
sed -i "s/ (G022).*//g" V.bombi_SC_orthologs_microsporidia

# Run annotation parser 


# Run enrichment analysis
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/RunGoEnrichTests.pl' --annot_file_1 'SC_genes.annot' --annot_file_2 '20230324_annotations_with_enzyme_names.annot' --obo_file '20230403_goslim_metagenomic_obo.txt' > SC_ortholog_enrichment