#!/usr/bin/env perl

use 5.038;
use warnings FATAL => 'all';
use autodie ':default';
use Spreadsheet::XLSX;
use Util qw(list_regex_files json_file_to_ref violin_plot workbook_to_hash);# ref_to_json_file);
use Getopt::Long 'GetOptions';
#use latex qw(write_2d_array_to_tex_tabular);

=purpose
	this script shows the variety of targetable regional lengths 
=cut

my (%group_species, %files);
my $figwidth = 11;
GetOptions(
#	'blast-json=s'  => \$files{blast_json},
#	'msa-program=s' => \$msa_program,
#	'gene-list|g=s' => \$files{gene_list_file},
#	'help'          => \$help,
	'figwidth=s'	 => \$figwidth,
	'hosts=s'       => \$files{hosts},
#	'output-stem'	 => \$stem,
	'pathogens|p=s' => \$files{pathogens},
) or die "error with GetOptions: $!";

sub read_tsv {
	my $tsv_file = shift;
	my $group    = shift;
	if (not defined $tsv_file) {
		die "$group got an undefined file.";
	}
	open my $fh, '<', $tsv_file;
	while (<$fh>) {
		next if /^#/; # skip comment lines
		chomp;
		my @line = split /\t/;
		if (scalar @line != 2) {
			p @line;
			die "There must be exactly 2 items in the tab-delimited file $tsv_file at line $., but I got the above.";
		}
		push @{ $group_species{$group} }, $line[0];
	}
	return 1;
}

read_tsv( $files{pathogens}, 'pathogens' ); # writes to %group_species

my %regions;
foreach my $species (@{ $group_species{pathogens} }) {
	say $species;
	foreach my $xlsx (list_regex_files($species,'xlsx')) {
		say "\t$xlsx";
		my $wb = workbook_to_hash({ filename => $xlsx });
		$wb = $wb->{'worksheet A'};
		foreach my $row (0..scalar @{ $wb->{length} } - 1) {
			next if $wb->{'length'}[$row] < 15;
			push @{ $regions{$species} }, $wb->{'length'}[$row];
		}
	}
}
violin_plot({
	filename => 'svg/length_violin_species',
	data		=> \%regions,
	logscale	=> ['y'],
	xlabel	=> 'Species',
	title		=> 'Targetable Segment Lengths Per Pathogen',
	ylabel	=> 'Targetable Segment Length (a.a.)',
	figwidth	=> $figwidth
});
