#!/home/vlb19/opt/perl/ -w
use strict;
use Bio::SeqIO;

### mfba004@live.rhul.ac.uk 

my $usage = "usage: perl $0 <directory_of_fasta_folders> \n";
die $usage unless (@ARGV eq 1);
my $dir = $ARGV[0];

# Go through folder of fasta files 
my @files = <$dir/*.afa.fa>;
my %species_to_sequence;
my %species_to_ID_to_1;
foreach my $file (@files){
    # Open file and save sequence
    my $inseq = Bio::SeqIO->new('-file' => "<$file",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
        #warn "$id \n";
        # Check if ID exists for everything and save sequence for everything too 
        die "I've already seen $id $file \n" if (defined $species_to_ID_to_1{$id}{$file});
        $species_to_ID_to_1{$id}{$file}=1;

        $species_to_sequence{$id}.=$seq;

    }
}

# print sequences 
foreach my $species(keys %species_to_sequence){
    my $sequence = $species_to_sequence{$species};
    print ">$species\n$sequence\n";
}
