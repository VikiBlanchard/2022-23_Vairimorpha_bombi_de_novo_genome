#!/bin/bash 

"""Visualise taxonomic composition of sequence generated by Oxford Nanopore run"""

__appname__ = 'kraken2_run.sh'
__author__ = 'Viki Blanchard (viki.blanchard.2018@live.rhul.ac.uk)'
__version__ = '0.0.1'


""" GET ISOLATE NAME """ 

#########################################
### Set up Krona in conda environment ### 
#########################################

# Activate conda environment to run the analysis in and install Krona

conda create --yes -n krona krona 
conda activate krona 

# Delete symbolic link which is not correct
rm -rf ~/anaconda3/envs/krona/taxonomy 

# Create directory in home where the krona database is stored 
mkdir -p ~/krona/taxonomy

# Make a new symbolic link to the database directory 
ln -s ~/krona/taxonomy ~anaconda3/envs/krona/taxonomy 

#############################################
### Build the taxonomy database for Krona ###
#############################################

# Build taxonomy database using Krona commands. If this fails we can download a pre-built database 
ktUpdateTaxonomy.sh ~/krona/taxonomy

# Unzip the file 
gzip -d taxonomy.tab.gz

# Move the unzipped file to the taxonomy directory 
mv taxonomy.tab ~/krona/taxonomy

#############################################################
### Visualise species composition from the Kraken2 report ###
#############################################################

# Navigate into Kraken 2 analysis directory
cd ~/analysis/kraken 

# Build a two-column file from the Kraken 2 results to import into Krona  
cat ${isolate}.kraken | cut -f 2,3 > ${isolate}.kraken.krona

# Import file into krona 
ktImportTaxonomy ${isolate}.kraken.krona

# Generate html report 
firefox taxonomy.krona.html

################################################################
### Visualise species composition from the Centrifuge report ###
################################################################

# Navigate into Centrifuge analysis directory
cd ~/analysis/centrifuge 

# Build a two-column file from the Centrifuge results to import into Krona  
cat ${isolate}-results.txt | cut -f 2,3 > ${isolate}-results.krona

# Import file into krona 
ktImportTaxonomy ${isolate}-results.krona

# Generate html report 
firefox taxonomy.krona.html