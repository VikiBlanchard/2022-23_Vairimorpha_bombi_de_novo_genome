#!/usr/bin/env bash
# KAT is a toolkit for addressing assembly completeness through k-mer counts (Mapleson et al.,
2017)
# More information about KAT in: https://github.com/TGAC/KAT
# You can use the short-read DNA sequencing data provided in the Key Resource Table (Accession
number: SRX8624462) to run the following script
# You need to provide the file path to the sequencing data or run this script in the same folder
where the sequencing data is saved
kat hist -o prefix -t 10 SRR12099722* 1> kat.output.txt
echo dme_size >> genome_size.txt
grep -i "Estimated" kat.output.txt >> genome_size.txt
# hist: a kat module for drawing histograms and estimating genome size
# -o: output prefix; you can specify ‘‘prefix’’ for your species or strain names
# -t: the number of threads that will be used to run the kat program
# You can replace SRR12099722* with your short-read DNA sequencing data# You can replace dme_-
size with the name of your species
# You can check the kat output by typing ‘‘cat genome_size.txt’’ in your terminal
> cat genome_size.txt
# Genome size can be estimated using the short-read DNA sequencing data
dme_size
Estimated genome size: 166.18 Mbp
Estimated heterozygous rate: 0.41%