#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use MutationTools::read_SAM;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -s <SAM>\n
Optional: -m Memory intensive mode to calculate transcript/fragment count (y/n) [n]\n";
our($opt_s, $opt_m);
getopt('sm');
die $usage unless ($opt_s);
die "Cannot open $opt_s : $!\n" unless (-e $opt_s);
if(!defined $opt_m) { $opt_m = 'n'; }

# Parse SAM file
my $SAM_summary = samlines::SAM_to_summary_hash($opt_s, $opt_m);

# Calculate DOC
my $depth_of_coverage = "Unknown";
my $depth_of_coverage_non_unique = "Unknown";
my $depth_of_no_coverage = "Unknown";
if(defined $$SAM_summary{'total_reference_sequence_length'}) { 
	$depth_of_coverage = sprintf("%.2f", ($$SAM_summary{'nt_aligned'} / $$SAM_summary{'total_reference_sequence_length'})); 
	$depth_of_coverage_non_unique = sprintf("%.2f", ($$SAM_summary{'nt_aligned_non_unique'} / $$SAM_summary{'total_reference_sequence_length'})); 
	$depth_of_no_coverage = sprintf("%.2f", ($$SAM_summary{'nt_not_aligned'} / $$SAM_summary{'total_reference_sequence_length'})); 
}

# Make summary
my $summary = "Total reads: $$SAM_summary{'total_reads'}\n";
$summary .= "Total nucleotides: $$SAM_summary{'total_nt'}\n\n";
$summary .= "Reads aligned: $$SAM_summary{'reads_aligned'}\n";
$summary .= "Nucleotides aligned: $$SAM_summary{'nt_aligned'}\n";
$summary .= "Depth of coverage: $depth_of_coverage\n\n";
$summary .= "Reads aligned (non-unique/secondary): $$SAM_summary{'reads_aligned_non_unique'}\n";
$summary .= "Nucleotides aligned (non-unique/secondary): $$SAM_summary{'nt_aligned_non_unique'}\n";
$summary .= "Depth of coverage (non-unique/secondary): $depth_of_coverage_non_unique\n\n";
$summary .= "Reads not aligned: $$SAM_summary{'reads_not_aligned'}\n";
$summary .= "Nucleotides not aligned: $$SAM_summary{'nt_not_aligned'}\n";
$summary .= "Depth of no coverage: $depth_of_no_coverage\n\n";

# Print summary
warn "$summary\n";
my $output = ($opt_s . '-reads-aligned-to-which-contigs3.tab');
open my $ofh, '>', $output or die "Cannot open $output: $!\n";
print $ofh "contig\tread_tally\ttranscript_aka_fragment_tally\tnt_tally\tDepth\n";
foreach my $contig(sort keys %{$$SAM_summary{'reads_aligned_to_contigs'}}) {
	my ($tally, $nt_tally, $frag_tally, $depth_of_contig_coverage) = (0,0,'not_calculated','not_calculated');
	if(defined $$SAM_summary{'reads_aligned_to_contigs'}{$contig}) { $tally = $$SAM_summary{'reads_aligned_to_contigs'}{$contig}; }
	if(defined $$SAM_summary{'nt_aligned_to_contigs'}{$contig}) { $nt_tally = $$SAM_summary{'nt_aligned_to_contigs'}{$contig}; }
	if(defined $$SAM_summary{'transcripts_aligned_to_contigs'}{$contig}) { $frag_tally = $$SAM_summary{'transcripts_aligned_to_contigs'}{$contig}; }
	if(defined $$SAM_summary{'reference_sequence_length'}{$contig}) { $depth_of_contig_coverage = ($nt_tally / $$SAM_summary{'reference_sequence_length'}{$contig}); }
	print $ofh "$contig\t$tally\t$frag_tally\t$nt_tally\t$depth_of_contig_coverage\n"; 
}
print $ofh "\n$summary\n";
close $ofh;
