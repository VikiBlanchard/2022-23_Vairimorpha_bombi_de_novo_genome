#!/usr/bin/perl -w 
use strict ; 

my $usage = "perl $0 <file> > file_out \n"; 

die $usage unless (@ARGV eq 1); 

my $count = 1; 
open my $fh, '<', $ARGV[0]; 
while (my $line = <$fh>) {
	chomp $line; 
	print "$line$count\n"; 
	$count++; 
}
