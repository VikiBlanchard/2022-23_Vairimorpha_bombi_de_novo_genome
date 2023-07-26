#!/bin/bash 
## Usage: Tree_file_maker.sh
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.webster.2018@live.rhul.ac.uk)

# Add correctly formatted species names to each sequence in all the pep files
for pep_file_name in */*.pep; do 
    sed -e "/>/s/$/pep_organism_name="$(basename -- ${pep_file_name} .annotation.pep)"/" -i $pep_file_name
done 



##################################################################
### Find which orthologs genes are the same across the species ###
##################################################################

# Save location of all_orthomcl.out file 
all_orthomcl="./all_orthomcl.out"

# Save single and multi-copy orthogroups that are shared by all microsporidia 
grep "36 taxa" $all_orthomcl > all_shared_orthogroups

# Save single copy orthogroups that are shared by all microsporidia 
grep "(36 genes,36 taxa)" $all_orthomcl > shared_sc_orthogroups
# add space to orthogroupID for splitting later
sed -i "s/(/ (/g" shared_sc_orthogroups

# Make directory for ortholog files 
mkdir SC_orthogroups 
cd SC_orthogroups

# Split orthomcl file into 1 file per orthogroup 
awk '{ print > $1 ; close($1) }' ../shared_sc_orthogroups

# Format all orthogroup files in directory
sed -i "s/.*:	 //g" * # Remove orthogroup match information 
sed -i "s/([^)]*)/\\n/g" * # Put gene IDs on separate lines and remove orthoID from end of each line
sed -i "s/.*|/>/g" * # Remove orthoID from the start of each line 
sed -i '/^$/d' * # Remove empty lines


##########################################
### Get pep sequences for each gene ID ###
##########################################

# Save location of files from Synima run 
synima_run_files='../'

# Make _peps file for each orthogroup containing the corresponding pep IDs and sequences
for orthogroup in *; do 
        cat ${orthogroup} | while read pep_id; do 
        # look for pep ID in file, get matching pep entry, and save it to a new file
        sed -n "/^${pep_id}/,/^>/p" ${synima_run_files}/*/*.pep | sed \$d >> $(basename -- $orthogroup)_peps
        done
done 

# Move organism name to the start of each line in the file
sed -i '/pep_organism_name=/s/.*pep_organism_name=/>/g' *_peps
# Remove whitespace from the end of lines 
sed -i 's/ *$//' *_peps
# Replace dashes with underscores
sed -i '/>/s/-/_/' *_peps
sed -i 's/-/_/g' *_peps


########################################
### Align the orthologs using muscle ###
########################################

# Run a Muscle alignment for every ortholog ID
for orthogroup in *_peps; do
    muscle -in ${orthogroup} -out $(basename ${orthogroup}).afa
done

########################################
### Trim the alignments using trimAl ###
########################################

# Trim every muscle alignment using the parameters from Huang et al. 2021
for alignment in *.afa; do
    trimal -in $alignment -w 3 -gt 0.95 -st 0.01  -fasta -out $(basename $alignment).fa -htmlout $(basename $alignment).html  # generate nexus files for further processing and html files to view the results
done

# -w = This value establishes the number of columns at each side of a given position that trimAl has to consider when computing some scores, such as gap, conservation or consistency scores, for that position. When a window size is given, trimAl provides the average value of all columns considered.
# -gt = Sets the limit for an acceptable gap score (The gap score for a column of size n is the fraction of positions in the column without a gap). Columns below this score are deleted from the input alignment
# -st = Similarity score. Measures similarity value for each column from the alignment using the Mean Distance method (Blosum62 similarity matrix)


###################################################
### Concatenate files for phylogenetic analysis ###
###################################################

# Run concatenation code written with rfarrer to pull all matched orthologs together
perl '/home/vlb19/Documents/Coding/2022-23_Vairimorpha_bombi_de_novo_genome/code/folder_of_fastas_to_cat.pl' . > "All_matched_orthologs_across_species.fa"

# Run FASTA-parser code from rfarrer to convert the FASTA file into a NEXUS file for MrBayes
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s "All_matched_orthologs_across_species.fa" -p nexus -d protein > "All_matched_orthologs_across_species.nexus"

# Run FASTA-parser code from rfarrer to convert the FASTA file into a PHYLIP file for RAxML
perl '/home/vlb19/Documents/Coding/Downloaded_Repositories/2022_Farrer_Lab_Code/perl_scripts/FASTA-parser.pl' -s "All_matched_orthologs_across_species.fa" -p phylip -d protein > "All_matched_orthologs_across_species.phylip"

###################
### Run MrBayes ### 
###################

mb -i "All_matched_orthologs_across_species.nexus"
prset aamodelpr = mixed 
mcmc ngen = 1000000 samplefreq=10
# wait till finished then run: 
sump burninfrac=0.25
sumt burninfrac=0.25

#################
### Run RAxML ###
#################

raxmlHPC -p 7487643260 -f a -x 1784414621 -# 1000 -m PROTGAMMAILGF -n bootstrapped_raxml_tree -s All_matched_orthologs_across_species.phylip -T 2
# -f a runs the program in bootstrapping mode
# -p is random seed number to make your tree reproducible 
# -x is random number seed for bootstrapping to make your analysis reproducible
# -# sets the number of distinct starting trees 
# -s is the alignment file 
# -n is the output file 
# -T is the number of cores (optional - Rhys recommends no more than 3)
# -m sets the protein model (LG+I+F+G selected from ProtTest)


# open figtree 
figtree 




##########################
### Useless extra code ### 
##########################


# Remove bp count from pep sequences 
sed -i 's/ .*/ /' *.afa.fa

# Remove newline characters
sed -i -z 's/[\n\r]//g' *.afa.fa

# Add newlines to separate species
sed -i 's/>/\n>/g' *.afa.fa

# Sort pep files by species
for pep_file in *.afa.fa; do 
    sort -o $pep_file $pep_file
done

# Create function to merge all pep files into one file
function multijoin() {
    out=$1
    shift 1
    cat $1 | awk '{print $1}' > $out
    for f in $*; do join $out $f > tmp; mv tmp $out; done
}

# Merge all pep files into one file and remove empty first line
multijoin All_matched_orthologs_across_species.fa *afa.fa 

# Remove empty first line 
sed -i '1d' All_matched_orthologs_across_species.fa

# Remove extra spaces 
sed -i 's/ /;/' All_matched_orthologs_across_species.fa
sed -i 's/ //g' All_matched_orthologs_across_species.fa
sed -i 's/;/\n/g' All_matched_orthologs_across_species.fa

# Add newlines to species
sed -i 's/ /\n>/g' All_matched_orthologs_across_species.fa


grep \> All_matched_orthologs_across_species.fa | wc -l 
24
