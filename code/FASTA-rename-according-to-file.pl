#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

# Opening commands
my $usage = "Usage: perl $0 <fasta> <file of old to new names>\n";
die $usage if (@ARGV != 2);
my ($fasta, $file_of_ids) = @ARGV;
warn "Opening $ARGV[0]...\n";

# save file to memory
my $ids = &save_ids($file_of_ids);


# run through FASTA file and print with new ids

my $inseq = Bio::SeqIO->new('-file' => "<$fasta",'-format' => 'fasta' ) ;
while (my $seq_obj = $inseq->next_seq ) {
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my $desc = $seq_obj->description;
	
	$seq =~ s/(\S{60})/$1\n/g;

	# get new ids

	die "nothing found for id $id\n" if(!defined $$ids{$id});
	my $new_id = $$ids{$id};

	print ">$new_id\n$seq\n";

}

sub save_ids {
	my $file = $_[0];
	my %old_to_new;
	open my $fh, '<', $file;
	while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($old, $new) = @bits;

		$old_to_new{$old} = $new;
	}
	return \%old_to_new;
}