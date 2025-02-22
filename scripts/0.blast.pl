#!/usr/bin/env perl

use strict;
use warnings FATAL => 'all';
use autodie ':default';
use FindBin 1.51 qw( $RealBin ); # use local versions so users don't have to install
use lib $RealBin;
use util qw(json_file_to_ref ref_to_json_file execute);
use Getopt::Long 'GetOptions';
use feature 'say';
use DDP;
my %hits;

sub tblastn_or_blastp {
	my $file	= shift;
	my $stem = $file;
	$stem =~ s/\.(?:fna|faa|fasta)$//;
	if (not grep {not -f "$stem.$_"} ('pdb','phd','phi','phr','pin','pjs','pog','pot','psq','ptf','pto')) {
		return 'blastp';
	}
	if (not grep {not -f "$stem.$_"} ('ndb')) {#,'nhr','nhd','njs','nsq','nhi','nto','nin','ntf','nog')) {
		return 'tblastn';
	}
	die "Couldn't find a blast type (tblastn/blastp) for $file; check that file's directory to see if \"makeblastdb\" was done.";
}

sub simplify_json {#($json) { # and write to %hits
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
		if ($best_evalue >= 0.1) {
			print "$species;$query\n";
			print "$best_evalue\n";
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
my (%group_species, %db);
sub read_tsv {
	my $tsv_file = shift;
	my $group    = shift;
	my %return;
	open my $fh, '<', $tsv_file;
	while (<$fh>) {
		next if /^#/;
		chomp;
		my @line = split /\t/;
		if (scalar @line != 2) {
#			p @line;
			die "There must be exactly 2 items in the tab-delimited file $tsv_file at line $., but I got the above.";
		}
		unless ((-f $line[1]) || (-l $line[1])) {
			die "$line[1] isn't a file, which was read from $tsv_file at line $.";
		}
		$db{$line[0]} = $line[1];
		push @{ $group_species{$group} }, $line[0];
	}
	return 1;
}

sub usage {
	my %h = (
#		'blast-json'=> 'output: a json summary of ',
		'query-list' => 'a fasta file, listing all of the genes ',
		hosts       => 'a tab-delimited, two-column file of host species and their proteome/genome files',
		pathogens   => 'a tab-delimited, two-column file of pathogen species and their proteome/genome files',
	);
	print "----------\n";
	print 'Execution: perl ' . __FILE__ . " --hosts hosts_file --pathogens pathogen_file\n";
	print "----------\n";
	p %h;
	print "----------\n";
}

my (%files, $help);
GetOptions(
#	'blast-json=s'  => \$files{blast_json},
	'query-list|q=s'	=> \$files{gene_list_file},
	'help'          	=> \$help,
	'hosts=s'       	=> \$files{hosts},
	'pathogens|p=s'	=> \$files{pathogens},
) or die 'error with GetOptions';

my @undef_files = grep {not defined $files{$_}} ('gene_list_file', 'hosts', 'pathogens');
if (scalar @undef_files > 0) {
	p @undef_files;
	say 'the above files are necessary';
	usage();
	die;
}
my $home = execute('echo $HOME', 'stdout');
my $blastp = execute('which blastp', 'stdout');
die 'no blastp' if $blastp eq '';
read_tsv( $files{hosts}, 'hosts' );				# writes to %db
read_tsv( $files{pathogens}, 'pathogens' );	# writes to %db
foreach my $s (sort keys %db) {
	$db{$s} =~ s/\.(?:bz2|gz)$//; # remove compression suffix if present
	$db{$s} =~ s/\.[^\.]+$//; # remove the suffix, e.g. ."fa" for each file to get the stem
}
my $outdir = 'blast.results';
mkdir $outdir unless -d $outdir;
my @out;
foreach my $s (sort keys %db) {
	my $out = "$outdir/$s.json";
	push @out, $out;
	next if ((-f $out) && ((-s $out) > 21));
	my $blast_type = tblastn_or_blastp($db{$s});
	execute("$blast_type -query $files{gene_list_file} -out $out -outfmt 15 -num_threads 20 -max_target_seqs 4000 -db $db{$s}");
}

foreach my $json (@out) {
	say $json;
	simplify_json( $json );
}
mkdir 'json' unless -d 'json';
ref_to_json_file(\%hits, 'json/blast.hits.json');
