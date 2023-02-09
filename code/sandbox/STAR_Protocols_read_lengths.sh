#!/usr/bin/env bash
# Create a new file and generate a header line
echo "platform,length" > length.csv
# Add each read length into the length.csv file.
bioawk -c fastx ’{print "PacBio_CLR," length($seq)}’ SRR11906525_WGS_of_drosophila_melanogaster_female_adult_subreads.fastq.gz >> length.csv
bioawk -c fastx ’{print "PacBio_HiFi," length($seq)}’ SRR12473480_Drosophila_PacBio_HiFi_UltraLow_subreads.fastq.gz >> length.csv
bioawk -c fastx ’{print "ONT," length($seq)}’ SRR13070625_1.fastq.gz >> length.csv