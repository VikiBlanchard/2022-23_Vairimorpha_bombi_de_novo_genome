#!/usr/bin/perl -w
use strict;

### rfarrer and vikiblanchard

# opening commands
my $usage = "perl $0 <sam file> > my.fastq\n";
die $usage unless(@ARGV eq 1);
my $sam = $ARGV[0];

# open sam and find unaligned reads
open my $fh, '<', $sam or die "Error: Cannot open $sam : $!";
while(my $line=<$fh>) {
	chomp $line;

	# ignore header
	next if($line =~ m/^\@/);

	my @bits = split /\t/, $line;
	my ($qname, $bitwise, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = @bits;
	
	# ignore aligned
	next if($rname ne '*');

	# print fastq
	print "\@$qname\n$seq\n+\n$qual\n";
}
