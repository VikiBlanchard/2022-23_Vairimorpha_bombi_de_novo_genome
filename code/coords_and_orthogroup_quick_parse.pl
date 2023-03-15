#!/usr/bin/perl -w
use strict;
use Data::Dumper;

# usage
my $usage = "perl $0 <coords> <.clusters> > out_summary\n";
die $usage unless(@ARGV eq 2);
my ($coords, $clusters) = @ARGV;

# save coords
my ($contig_to_gene_count, $gene_to_contig) = &save_coord_info($coords);

# save clusters
my $contig_to_results = &save_orthogroups($clusters, $gene_to_contig);

# summary per contig

#header
print "contig\tgene_count";
foreach my $key_name(sort keys %{$contig_to_results}) {
    print "\t$key_name";
}
print "\n";

#data
foreach my $contig(sort keys %{$contig_to_gene_count}) {
    print "$contig";

    # gene count
    my $gene_count = $$contig_to_gene_count{$contig};
    print "\t$gene_count";

    foreach my $key_name(sort keys %{$contig_to_results}) {
        my $value = 0;

        if(defined $$contig_to_results{$key_name}{$contig}) {
            $value = $$contig_to_results{$key_name}{$contig};
        }
        print "\t$value";
    }
    print "\n";
}

sub save_orthogroups {
    my ($file, $gene_to_contig) = @_;

    my %results;

    my $interest = 'Vairimorpha_bombi_8.1-3no1x';

    my %info;
    open my $fh, '<', $file or die "Cannot open $file : $!";
    while(my $line=<$fh>) {
        chomp $line;
        next if($line eq '');

        my @bits = split /\t/, $line;
        my ($orthogroup, $species, $ignore1, $gene) = @bits;

        die "line = $line\n" if(!defined $orthogroup);

        # does orthogroup contain Vaira or other things
        if($species eq $interest) {
            $info{$orthogroup}{'contains_vaira'} = 1;
            $info{$orthogroup}{'genes'} .= "$gene\n";
        } else {
            $info{$orthogroup}{'contains_other'} = 1;
        }

        # total species
        $info{$orthogroup}{'unique_species'}{$species}++;


    }
    close $fh;

    # figure out the important stuff!
    ORTHOGROUP: foreach my $orthogroup(sort keys %info) {
        next ORTHOGROUP if(!defined $info{$orthogroup}{'contains_vaira'});
        next ORTHOGROUP if(!defined $info{$orthogroup}{'contains_other'});

        # We have found an orthogroup in both vaira and another species
        #print "orthogroup = $orthogroup\n";
        #print Dumper($info{$orthogroup});
        
        # make a summary of the numbers of species in each orthogroup
        my $number_of_species_in_orthogroup = scalar(keys(%{$info{$orthogroup}{'unique_species'}}));

        # max number of paralogs across all orthogroups on a given contig
        my $max_number = 0;
        SPECIES: foreach my $species_name(sort keys %{$info{$orthogroup}{'unique_species'}}) {
            my $number_of_genes = $info{$orthogroup}{'unique_species'}{$species_name};
            if($number_of_genes > $max_number) { $max_number = $number_of_genes; }
        }

        

        #print "number of species = $number_of_species_in_orthogroup\n";
        #die;

        my $genes = $info{$orthogroup}{'genes'};
        my @gene_array = split /\n/, $genes;

        # try and make some kind of summary about paralogs in vaira orthogroups
        my $paralog_count = scalar(@gene_array);

            

        # look these up and find the contig for them
        my $counter2 = 0;
        foreach my $gene(@gene_array) {
            my $contig = $$gene_to_contig{$gene};

            $results{'gene_in_orthogroups_per_contig'}{$contig}++;

            # save max number of paralogs across all orthogroups
            if(!defined $results{'max_paralogs_per_orthogroup'}{$contig}) {
                $results{'max_paralogs_per_orthogroup'}{$contig} = $max_number;
            }
            elsif($max_number > $results{'max_paralogs_per_orthogroup'}{$contig}) {
                $results{'max_paralogs_per_orthogroup'}{$contig} = $max_number;
            }
            else {}

            $results{'species_per_orthogroup_per_contig'}{$contig} .= "$number_of_species_in_orthogroup,";

            # save num orthogroups
            if($counter2 eq 0) {
                

                $results{'paralogs_of_vaira_per_contig'}{$contig} .= "$paralog_count,";

                $counter2 = 1;
            }
        }
    }

    return \%results;
}

sub save_coord_info {
    my $file = $_[0];
    my %contig_info;
    my %gene_to_contig;
    open my $fh, '<', $file or die "Cannot open $file : $!";
    while(my $line=<$fh>) {
        chomp $line;
        my @bits = split /\t/, $line;
        my ($contig, $gene) = @bits;
        $contig_info{$contig}++;
        $gene_to_contig{$gene} = $contig;
    }
    close $fh;
    return (\%contig_info, \%gene_to_contig);
}