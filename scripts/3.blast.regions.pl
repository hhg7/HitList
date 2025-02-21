#!/usr/bin/env perl

use 5.038;
use warnings FATAL => 'all';
#use warnings::unused;
use autodie ':default';
use Util qw(json_file_to_ref ref_to_json_file list_regex_files execute);
#use latex qw(write_2d_array_to_tex_tabular);
use bioinf qw(fasta2hash hash2fasta_file);

=purpose
	this script looks through highlighted regions of possible target proteins and finds if the targeted regions are themselves likely to have BLAST hits.
=cut

unless (defined $ARGV[0]) {
	die 'a hosts file must be specified, e.g. hosts.tsv';
}

my %hits;
sub simplify_json { # and write to %hits
	my $json = shift;
	my $blast = json_file_to_ref($json);
	say "Finished loading $json";
	if (defined $blast->{BlastOutput2}) {
		$blast = $blast->{BlastOutput2}; # simplify hash
	}
	my $species = $json;
	$species =~ s/\.json$//;
	if ($species =~ m/\/([^\/]+)$/) {
		$species = $1;
	}
	$hits{$species}{report} = $blast->[0]{report};
	foreach my $report (@{ $blast }) {
		if (
				(defined $report->{report}{results}{search}{message})
				&&
				($report->{report}{results}{search}{message} eq 'No hits found')
			) {
			say STDERR "$report->{report}{results}{search}{query_title} has no hits found for $species";
			next;
		}
		my ($best_evalue, $best_h, $best_hsps) = ('inf', 'inf', 'inf');
		my $query = $report->{report}{results}{search}{query_title};
		if (not defined $query) {
			die "\$report->{report}{results}{search}{query_title} isn't defined.";
		}
=if ( # skip non-yeast genomes
			($query =~ m/DEG(\d+)/) &&
			(
				($1 < 20010001)
				||
				($1 > 20011110)
			)
		) {
			say "skipping $query";
			next;
		}
=cut
		$report = $report->{report}{results}{search}{hits}; # simplify hash
		while (my ($h, $hit) = each @{ $report }) {
			while (my ($hsps, $pair) = each @{ $hit->{hsps} } ) {
				next if $pair->{evalue} > $best_evalue;
				$best_evalue = $pair->{evalue};
				$best_h = $h;
				$best_hsps = $hsps;
			}
		}
=if (($species eq 'fusarium.sambucinum') && ($query eq 'DEG20010641')) {
			p $report;
			say "\$best_evalue = $best_evalue";
			say "\$best_h      = $best_h";
			say "\$best_hsps   = $best_hsps";
			die;
		}
=cut
		if ($best_evalue >= 0.1) {
			say "$species;$query";
			say $best_evalue;
			next;
		}
		if ($best_h == 'inf') {
			p $report;
			die "$query has nothing";
		}
		die if not defined $report->[$best_h]{description};
		$hits{$species}{$query}{description} = $report->[$best_h]{description};
		foreach my $stat (keys %{ $report->[$best_h]{hsps}[$best_hsps] }) {
			$hits{$species}{$query}{$stat} = $report->[$best_h]{hsps}[$best_hsps]{$stat};
		}
	}
}
my %db;
open my $fh, '<', $ARGV[0];
while (<$fh>) {
	next if /^#/;
	chomp;
	my @line = split /\t/;
	if (scalar @line != 2) {
		p @line;
		die "There must be exactly 2 items in the tab-delimited file $ARGV[0] at line $., but I got the above.";
	}
	unless (-f -r $line[1]) {
		die "$line[1] isn't a file or can't be read, which was read from $ARGV[0] at line $.";
	}
	$db{$line[0]} = $line[1];
}
close $fh;
my (%files, $help);

my $targetable_regions = json_file_to_ref('json/all.fungi.target.regions.json');
my @msa_files = list_regex_files('\.output\.aln$', 'fa');
#p $targetable_regions, array_max => 12;
my %targets;
foreach my $region (@{ $targetable_regions }) {
	my @fa_files = grep {/$region->[0]/} @msa_files;
	if (scalar @fa_files != 1) {
		p @fa_files;
		die "!= 1 fasta/aln file (above) for $region->[0], so a file to read is indeterminate.";
	}
	my $fasta = fasta2hash( $fa_files[0] );
	my $original_query = $fasta->{$region->[0]};
	my $n_aa_start = length $original_query;
	my @original_query = split '', $original_query;
	my $sub_query      = join ('', @original_query[$region->[9]..$region->[10]]);
	next if $sub_query =~ m/\*/; # skip matches that have stop codons inside
	my $aa_count = () = $sub_query =~ m/[A-Ya-y]/g; # all amino acid codes
	next if $aa_count < 15;
	my $key_name = join ('.', @{ $region }[0,9,10]);
	$targets{$key_name} = $sub_query;
}
p %targets;
foreach my $s (sort keys %db) { # remove the suffix, e.g. ."fa" for each file to get the stem
	$db{$s} =~ s/\.[^\.]+$//;
}
my $blast_regions_file = 'json/blast.regions.fa';
hash2fasta_file( \%targets, $blast_regions_file );
my @out;
foreach my $s (sort keys %db) {
	my $out = "json/$s.region.blast.json";
	push @out, $out;
	next if ((-f $out) && ((-s $out) > 21));
	my $blast = 'tblastn';
	$blast = 'blastp' if $db{$s} =~ m/protein/i;
	execute("$blast -query $blast_regions_file -out $out -outfmt 15 -num_threads 20 -max_target_seqs 4000 -db $db{$s}");
}

foreach my $json (@out) {
	say $json;
	simplify_json( $json );
}
mkdir 'json' unless -d 'json';
ref_to_json_file(\%hits, 'json/fungal.blastp.region.hits.json');
