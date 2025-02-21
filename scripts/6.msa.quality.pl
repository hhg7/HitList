#!/usr/bin/env perl

use 5.038;
use warnings FATAL => 'all';
use autodie ':default';
use JSON 'encode_json';
use List::Util qw(min max);
use Util qw(json_file_to_ref colored_table);# ref_to_json_file);list_regex_files
use File::Temp 'tempfile';
use Capture::Tiny 'capture';
use Term::ANSIColor;

sub first_letter ($string) {
	if ($string =~ m/^([A-Za-z])/) {
		return $1;
	}
}
=sub arg_check { # to save space in 
	my ($args) = @_;
	my @reqd_args = ('reqd.args', 'optional.args', 'ref.types', 'sub.name');
	my @undef_args = {!defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'The arguments above are necessary for proper function and weren\'t defined.';
	}
	unless (
				(ref $args->{'req.dargs'} eq 'ARRAY') &&
				(ref $args->{'optional.args'} eq 'ARRAY') &&
				(ref $args->{'ref.types'} eq 'HASH') &&
				(ref $args->{'sub.type'} eq '')
	) {
		p $args;
		my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
		die "Something is wrong with the args submitted to $current_sub";
	}
	
}
=cut

#the JSON that I'm downloading was made the by following blastp command:
#blastp -query clustal.DEG20010729.input.fa -subject clustal.DEG20010729.input.fa -out clustal.DEG20010729.blastp.json -outfmt 15 -num_threads 5
my $msa_blastp = json_file_to_ref('../who.critical/fa/clustal.DEG20010729.blastp.json');
my %data;
foreach my $query ( @{ $msa_blastp->{BlastOutput2} }) {
	foreach my $hit_list (@{ $query->{report}{results}{bl2seq} }) {
		my $query_title = $hit_list->{query_title};
		foreach my $hit (@{ $hit_list->{hits} }) {
			my $hit_title = $hit->{description}[0]{title};
			%{ $data{$query_title}{$hit_title} } = %{ $hit->{hsps}[0] };
		}
	}
}
#p %data;
my @pathogens = qw(cryptococcus.neoformans.B.3501A	cryptococcus.neoformans.grubii.H99	cryptococcus.neoformans.JEC21	candida.auris	aspergillus.fumigatus	candida.albicans
nakaseomyces.glabratus	histoplasma.capsulatum	candida.tropicalis	candida.parapsilosis);
my @hosts = qw( homo.sapiens oryza.sativa	solanum.tuberosum	zea.mays	glycine.max);
my (@data, @row_labels, @col_labels);
foreach my $key1 (@hosts, 'DEG20010729', @pathogens) {
	my @row;
	foreach my $key2 (@hosts, 'DEG20010729', @pathogens) {
		push @row, $data{$key1}{$key2}{score};
	}
	push @data, [@row];
	if ($key1 =~ m/
						([^\.]+)		# genus
						\.				# separator
						(.+)			# species
					/xx) {
		push @row_labels, ucfirst first_letter($1) . ".$2";
		push @col_labels, ucfirst first_letter($1) . '.' . first_letter($2) . '.';
	} elsif ($key1 =~ m/^DEG\d+$/) {
		push @row_labels, 'yeast';
		push @col_labels, 'yeast';
	} else {
		die "can't get labels for $key1";
	}
}
colored_table({
	'col.labels'	=> \@col_labels,
	filename 		=> '/home/con/identify.target/who.critical/svg/msa.qc/DEG20010729_msa_quality.svg',
	data				=> \@data,
	'row.labels'	=> \@row_labels,
	cb_label			=> 'Score',
	title				=> 'Multiple Sequence Alignment Score'
});
