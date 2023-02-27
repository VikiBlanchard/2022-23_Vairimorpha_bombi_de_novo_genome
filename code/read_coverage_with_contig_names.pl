#!/usr/bin/perl -w 
use strict ; 

my $usage = "perl $0 <file> <dp_threshold> > file_out \n"; 

die $usage unless (@ARGV eq 2); 

my $count = 1; 
open my $fh, '<', $ARGV[0]; 
while (my $line = <$fh>) {
	chomp $line; 
    
    my @bits = split /\t/, $line;
	my ($contig, $read_tally, $transcript_aka_fragment_tally, $nt_tally, $Depth) = @bits;

    # Ignore contigs with coverage more than depth threshold 
    if ($Depth > 1) {
        next; 
    }
	

	print "$contig\n";
}
