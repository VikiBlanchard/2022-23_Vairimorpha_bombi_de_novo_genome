#!/usr/bin/perl -w
use strict;

### rfarrer and vwebster

# opening commands
my $usage = "perl $0 <sam file> > my.fastq\n";
die $usage unless(@ARGV eq 1);
my $sam = $ARGV[0];

my $aligned_total = 0;
my $unaligned_total = 0;

my $aligned_total_nt = 0;
my $unaligned_total_nt = 0;

# open sam and find unaligned reads
open my $fh, '<', $sam or die "Error: Cannot open $sam : $!";
while(my $line=<$fh>) {
	chomp $line;

	# ignore header
	next if($line =~ m/^\@/);

	my @bits = split /\t/, $line;
	my ($qname, $bitwise, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = @bits;
	
	# ignore aligned
	if($rname ne '*') {
		$aligned_total++;
		$aligned_total_nt += length($seq);
		next;
	}

	$unaligned_total++;
	$unaligned_total_nt += length($seq);

	# print fastq
	print ">$qname\n$seq\n";
}

# reads
my $total = $unaligned_total + $aligned_total;
warn "total reads in SAM = $total\n";
warn "aligned reads = $aligned_total\n";
warn "unaligned reads = $unaligned_total\n";

my $percent_unaligned = sprintf "%0.2f", (($unaligned_total / $total) * 100);

warn "percent unaligned = $percent_unaligned\n\n\n";

# nt
my $total_nt = $unaligned_total_nt + $aligned_total_nt;
warn "aligned nt = $aligned_total_nt\n";
warn "unaligned nt = $unaligned_total_nt\n";

my $percent_unaligned_nt = sprintf "%0.2f", (($unaligned_total_nt / $total_nt) * 100);
warn "percent unaligned (nt) = $percent_unaligned_nt\n";
