#!/usr/bin/env perl

use 5.040;
use warnings FATAL => 'all';
use autodie ':default';
use Util 'execute';
use bioinf qw(hash2fasta_file fasta2hash);

# purpose
#	the output compressed .tar.bz2 of all test files is too large for GitHub.
#	So the size of the test directory must be reduced

my $makeblastdb = '/home/con/programs/ncbi-blast-2.16.0+/bin/makeblastdb';
die "Missing/not executable $makeblastdb" unless -f -x $makeblastdb;
my %fasta_files = (
	'b.graminis'	=> 'blastdb/blumeria.graminis/FungiDB-65_BgraminisTritici96224_AnnotatedProteins.fasta',
	'c.truncatum'	=> 'blastdb/colletotrichum.truncatum/GCF_014235925.1_CTRU02_protein.faa',
	'o.sativa'		=> 'blastdb/oryza.sativa/GCF_001433935.1_IRGSP-1.0_protein.faa',
);

my @bad_files = grep {not -f $_} values %fasta_files;
if (scalar @bad_files > 0) {
	p @bad_files;
	die 'the above files are missing';
}
@bad_files = grep {-s $_ == 0} values %fasta_files;
if (scalar @bad_files > 0) {
	p @bad_files;
	die 'the above files have 0 size';
}

foreach my $species ('b.graminis', 'c.truncatum', 'o.sativa') {
	say $species;
	my $blast_output_file = 'blast.results/total.' . ucfirst "$species.json";
	my %hits;
	open my $fh, '<', $blast_output_file;
	while (<$fh>) {
		next unless /^\s+"title":/;
		chomp;
		$_ =~ s/^\s+"title":\s+"//;
		$_ =~ s/"$//;
		$hits{$_} = 1;
	}
	close $fh;
	if (scalar keys %hits == 0) {
		die "$blast_output_file has no hits.";
	}
	my $fasta = fasta2hash( $fasta_files{$species} );
	printf("There are %u keys.\n", scalar keys %{ $fasta });
	foreach my $missing_hit (grep {not defined $hits{$_}} keys %{ $fasta }) {
		delete $fasta->{$missing_hit};
	}
	printf("There are %u keys.\n", scalar keys %{ $fasta });
	my $trimmed_fasta_file = "blastdb/$species/trimmed.source.fasta";
	mkdir "blastdb/$species" unless -d "blastdb/$species";
	hash2fasta_file( $fasta, $trimmed_fasta_file);
	my $stem = $trimmed_fasta_file;
	$stem =~ s/\.fasta$//;
	my $cmd = "$makeblastdb -in $trimmed_fasta_file -out $stem -dbtype 'prot' -hash_index";
	say $cmd;
	execute($cmd);
}

