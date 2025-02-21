#!/usr/bin/env perl

use 5.038;
use warnings FATAL => 'all';
use autodie ':default';
use Util qw(ref_to_json_file json_file_to_ref);
use bioinf qw(fasta2hash);
use File::Temp 'tempfile';

=purpose
	this script is going to analyze the sub-region host-blast output data
=cut

sub msa_from_hash ($fasta) {
	
}

my $sub_region_blast = json_file_to_ref('json/fungal.blastp.region.hits.json');
my $start_query = fasta2hash('DEG20.edited.yeast.only.fa');
#p $sub_region_blast;
foreach my $species (sort keys %{ $sub_region_blast }) {
	my $species_name = $species;
	$species_name =~ s/\.region\.blast$//;
	say $species_name;
	foreach my $query (sort keys %{ $start_query }) {
		say $query;
		my @regions = grep {/^$query/} sort keys %{ $sub_region_blast->{$species} };
		next if scalar @regions == 0;
		p @regions;
		die;
	}
}
