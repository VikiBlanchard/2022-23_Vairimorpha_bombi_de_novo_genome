#!/usr/bin/perl -w
use Cwd;
my $path = Cwd::cwd();
print "$path\n";

use strict;
use Bio::SeqIO;
use Getopt::Std;
use Encode;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use MutationTools::read_FASTQ;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -s sequence.fastq > output\n
Output format:   -p Printing options (fasta, nexus, tba, color, fastq, none) [none]
                 -d If printing nexus (dna/protein), if printing TBA (source organism without dashes)  [dna]
		 -z Split fastq (color/base) into x files (n=ignore/all to output, 1=1 sequence per file etc.) [n]\n
Rename:          -a Append word to IDs (n=none) [n]
                 -b Replace IDs for Descriptions (yn) [n]
                 -c If opt_b, which part of Description is wanted (split by space)? [1]
		 -e Drop words from IDs (Separated by comma) E.g. supercontig,_,pilon,1. [n]\n
Reorder:         -n Make random subset of sequences/reads [0]
Modify sequence: -f Truncate x nt from 5' (left) [0]
                 -t Truncate x nt from 3' (right) [0]
	         -m Remove if smaller than this number of characters [1000000000000]
	         -i IDs seperated by comma (n=none): [n]
	         -r Reverse compliment on opt_i (y/n) [n]
		 -l Drop contigs in file []\n
Summarise:       -g Print summary for file (y/n) [y]
                 -h Print summary per entry (n=none, s=short (id and length) l=long (id, desc, seq, length etc.)) [n]\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_f, $opt_g, $opt_h, $opt_i, $opt_l, $opt_m, $opt_n, $opt_p, $opt_r, $opt_s, $opt_t, $opt_z);
getopt('abcdefghilmnprstz');
if(!defined $opt_a) { $opt_a = 'n'; }
if(!defined $opt_b) { $opt_b = 'n'; }
if(!defined $opt_c) { $opt_c = 1; }
if(!defined $opt_d) { $opt_d = 'dna'; }
if(!defined $opt_e) { $opt_e = 'n'; }
if(!defined $opt_f) { $opt_f = 0; }
if(!defined $opt_g) { $opt_g = 'y'; }
if(!defined $opt_h) { $opt_h = 'n'; }
if(!defined $opt_i) { $opt_i = 'n'; }
if(!defined $opt_m) { $opt_m = 1000000000000; }
if(!defined $opt_n) { $opt_n = 0; }
if(!defined $opt_p) { $opt_p = 'none'; }
if(!defined $opt_r) { $opt_r = 'n'; }
if(!defined $opt_t) { $opt_t = 0; }
if(!defined $opt_z) { $opt_z = 'n'; }
die $usage unless ($opt_s);
die "Cannot open $opt_s : $!" unless (-e $opt_s);
die "If defined -p, must be fasta, nexus, tba or none: $opt_p\n" if ($opt_p !~ m/fasta|nexus|none|tba|color|fastq/);
die "Setting -g must be y or n\n" if (($opt_g ne 'n') && ($opt_g ne 'y'));
die "Setting -h must be l, s or n\n" if (($opt_h ne 'n') && ($opt_h ne 'l') && $opt_h ne 's');

# Save sequences
my $fastq_struct = fastqfile::fastq_to_struct($opt_s);

# Random subsets? (if 0 returns full)
$fastq_struct = fastqfile::fastq_struct_to_random_subset($fastq_struct, $opt_n);

### Rename
# Append to ID
if($opt_a ne 'n') { $fastq_struct = fastqfile::fastq_struct_append_id($fastq_struct, $opt_a); }
# Replace ID with Description
if($opt_b ne 'n') { $fastq_struct = fastqfile::fastq_struct_replace_id_with_desc($fastq_struct, $opt_c); }
# Drop words from ID
if($opt_e ne 'n') { $fastq_struct = fastqfile::fastq_struct_remove_words_from_id($fastq_struct, $opt_e); }

### Modify
# Truncate and/or remove small sequences
if(($opt_f > 0) || ($opt_t > 0) || ($opt_m < 1000000000000)) { $fastq_struct = fastqfile::fastq_struct_truncate($fastq_struct, $opt_f, $opt_t, $opt_m); }
# Reverse compliment
if($opt_i ne 'n') { $fastq_struct = fastqfile::fastq_struct_reverse_compliment($fastq_struct, $opt_i, $opt_r); } 
# Remove reads/entries
if($opt_l) { $fastq_struct = fastqfile::fastq_struct_remove_entries($fastq_struct, $opt_l); }

### Summarise and print
if(($opt_g ne 'n') || ($opt_h ne 'n')) { $fastq_struct = fastqfile::fastq_struct_summary($fastq_struct, $opt_g, $opt_h); }
fastqfile::fastq_struct_print($fastq_struct, $opt_p, $opt_d, $opt_z);
