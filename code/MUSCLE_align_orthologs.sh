#!/bin/bash 
## Usage: Align orthologs 
##
## Options:
##   -h, --help    Display this message.
##   -n            Dry-run; only show what would be done.

### vlwebster (viki.blanchard.2018@live.rhul.ac.uk)

##################################
### Set up working environment ###
##################################

# Store isolate name as variable
isolate=${1}

# Set path directory for data input files
orthologs_dir=${2}

# Initialise output directories
mkdir results_${isolate}/$isolate_${3}



######################################################
### Make one file of all species for each ortholog ###
######################################################


######################################################
### Align each ortholog across species with MUSCLE ###
######################################################


# Concatenate all muscle alignments into a single file



#######################
### Trim alignments ###
#######################

trimal -in example1 -out output1 -htmlout output1.html -gt 1



################################################
### Generate tree from MUSCLE alignment only ###
################################################



############################################
### Generate tree from trimmed alignment ###
############################################