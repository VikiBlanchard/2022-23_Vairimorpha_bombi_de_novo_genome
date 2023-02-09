#!/usr/bin/env python3
import os
import sys
import pandas as pd 
from pandas import *

# Author: vlwebster viki.blanchard.2018@live.rhul.ac.uk
# Script: 4.2.get_BUSCO_matched_sequences.py
# Date: 20220813

""" Print BUSCO-matched sequences to text file """

#####################################################################
### Make text file of sequences that contain microsporidia BUSCOs ###
#####################################################################

# Save path to results directory 
results_dir = sys.argv[1]
#results_dir = "results_Vairimorpha_bombi_8_1/4.isolated_microsporidia_sequences/"

# Import BUSCO results table
busco_results_df = pd.read_csv(os.path.join(results_dir,"busco_results_table.tsv"), sep='\t')  

# Remove all missing BUSCOs 
busco_matches_df = busco_results_df[busco_results_df.Status != 'Missing']

# Remove all repeated sequence IDs 
unique_busco_matches = busco_matches_df.drop_duplicates(subset = ["Sequence"])

# Write BUSCO-matched sequence IDs into test file
with open (os.path.join(results_dir,'busco_match_sequences.txt'), 'w') as busco_match_sequences:
    for matched_sequence in unique_busco_matches['Sequence'].tolist():
        busco_match_sequences.write("%s\n" %matched_sequence)


