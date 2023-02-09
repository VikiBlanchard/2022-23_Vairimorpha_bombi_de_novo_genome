#!/usr/bin/env python3

# Author: Tao https://www.biostars.org/u/7761/ and viki.webster (viki.webster.2018@live.rhul.ac.uk)
# Script: remove_reads_from_fastq.py
# Source: https://www.biostars.org/p/199946/#254022  
# Date: 20220813

import fileinput
import sys
import os
reads_to_remove_path = sys.argv[1]

# Read in list of reads to remove
fp = open(reads_to_remove_path)
remove_list = []
while True:
    line = fp.readline()
    line = line.strip()
    if not line:
        break
    remove_list.append(line)

# Convert list to set
remove_list = set(remove_list)

# Set start switch to on 
switch = 1;
# Initialise line count for file
line_count = 0

# Open fastq file
for line_in in sys.stdin:
    
    line_count += 1
    # Remove any leading or trailing whitespace from  line
    line = line_in.strip()
    # Each sequence is 4 lines long, beginning with an ID starting @    
    # Check if line is an ID  
    if line_count%4 == 1:
        if line.startswith('@'):
            # If the sequence ID is in the list of reads to remove
            if line in remove_list:
                # Set switch to off 
                switch = 0
            else:
                # Set switch to on 
                switch = 1
        else:
            print ("error read Id", line)
            sys.exit(1)

    # If the switch is on, add line to the new file
    if switch == 1:
        print (line)
