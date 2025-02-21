#!/usr/bin/env perl

use FindBin 1.51 qw( $RealBin );
use lib $RealBin;
use strict;
use feature 'say';
use warnings FATAL	=> 'all';
use autodie ':default';
use util qw(msa_phylo_plot break_long_string_into_array fasta2hash json_file_to_ref);
use Getopt::Long 'GetOptions';
use List::Util 'min';
use Cwd 'getcwd';
use File::Temp 'tempfile';
use Term::ANSIColor;
use DDP;

sub usage {
	my %h = (
		'blast-json'=> 'a json summary of hits',
		'gene-key' 	=> 'a hash/dict JSON file of query titles',
		'query-list'	=> 'a fasta file, listing all of the genes ',
		'group'		=> 'the name of the group, e.g. "WHO Critical"',
		hosts      	=> 'a tab-delimited, two-column file of host species and their proteome/genome files',
		pathogens  	=> 'a tab-delimited, two-column file of pathogen species and their proteome/genome files',
	);
	print "----------\n";
	say 'Execution: perl ' . __FILE__ . ' --hosts hosts_file --pathogens pathogen_file';
	print "----------\n";
	p %h;
	print "----------\n";
}
my (%group_species, %db, %missing_files);
sub read_tsv { # wites to %db and %group_species
	my $tsv_file = shift;
	my $group    = shift;
	open my $fh, '<', $tsv_file;
	while (<$fh>) {
		next if /^#/;
		chomp;
		my @line = split /\t/;
		if (scalar @line != 2) {
			p @line;
			die "There must be exactly 2 items in the tab-delimited file $tsv_file at line $., but I got the above.";
		}
		unless (-f $line[1]) {
			$missing_files{$line[1]} = 1;
			next;
		}
		$db{$line[0]} = $line[1];
		push @{ $group_species{$group} }, $line[0];
	}
	return 1;
}

my (%files, $help, $group);
$files{query_species} = 'query species';
GetOptions(
	'blast-json=s' 		=> \$files{blast_json},
	'gene-key=s'			=> \$files{gene_key},
	'group=s'				=> \$group,
	'help'         		=> \$help,
	'hosts=s'      		=> \$files{hosts},
	'pathogens|p=s'		=> \$files{pathogens},
	'query-list|q=s'		=> \$files{query_list_file},
	'query-species|s=s'	=> \$files{query_species},
) or die 'error with GetOptions';

die 'group not defined' unless defined $group;
if ($help) {
	usage();
	die;
}
my $py = "$RealBin/x.my.align.py";
unless (-r -f $py) {
	die "$py doesn't exist.";
}
my @undef_files = grep {not defined $files{$_}} ('blast_json', 'query_list_file');#,'hosts');#, 'pathogens');
if (scalar @undef_files > 0) {
	p @undef_files;
	print "the above files are necessary\n";
	usage();
	die;
}
my @missing_files =  grep {not -f $files{$_}} ('blast_json', 'query_list_file');#,'hosts');#, 'pathogens');
if (scalar @missing_files > 0) { # get the error ASAP, so as not to waste time
	p @missing_files;
	print "the above files are necessary, but are either not files or don't exist\n";
	usage();
#	die;
}
if( defined $files{hosts}) {
	read_tsv( $files{hosts},         'hosts' ); # writes to %db, %group_species, & %missing_files
}
read_tsv( $files{pathogens}, 'pathogens' ); # writes to %db, %group_species, & %missing_files
if (scalar keys %missing_files > 0) {
	my @missing_files = sort keys %missing_files;
	p @missing_files;
	die 'missing files are listed above.';
}
my $all_blast = json_file_to_ref($files{blast_json});
my %all_fasta;
foreach my $species (keys %db) {
	$all_fasta{$species} = fasta2hash($db{$species});
}
foreach my $species (grep {not defined $db{$_}} keys %{ $all_blast }) {
	delete $all_blast->{$species};
}
my $gene_key;
if (grep {defined $files{$_}} 'gene_key') { # grep defined prevents autoinitialization
	$gene_key = json_file_to_ref( $files{gene_key} );
}
my $all_query_seq = fasta2hash( $files{query_list_file} ); # the list of query sequences
say '@ ' . __FILE__ . '& line' . __LINE__;
p $all_query_seq;
my %query; # get a list of all queries
foreach my $s (grep {defined $db{$_}} keys %{ $all_blast }) {
	say $s;
	foreach my $q (grep {$_ ne 'report'} keys %{ $all_blast->{$s} }) {
		$query{$q} = 1;
	}
}
foreach my $d ('svg', 'svg/msa', 'svg/phylo', 'tex', 'fa') {
	mkdir $d unless -d $d;
}
# make the input for clustalo
my $cwd = getcwd();
my @order;
foreach my $group ('pathogens', 'hosts') {
	if (defined $group_species{$group}) {
		push @order, @{ $group_species{$group} };
	}
}
say colored(['black on_yellow'], '@order:');
p @order;
foreach my $query (sort keys %query) {#'DEG20010160') {
	say "\$query = $query";
=foreach my $species (sort keys %{ $all_blast }) {
		say $species;
		if (defined $all_blast->{$species}{$query}) {
			say "\t$all_blast->{$species}{$query}{evalue}";
		} else {
			say "\tundef";
		}
	}
	die;
=cut
	my $msa_image  = "$cwd/svg/msa/$query" . '_msa.svg';
	my $tex_outfile = "$cwd/tex/$query.alignment.info.tabular.tex";
	my $out_fasta  = "$cwd/fa/clustal.$query.input.fa";
=if (
			(-f $msa_image)   &&
			(-f $tex_outfile) &&
			(-f $phyl_image)
		) {
		say "finished $query, no need to repeat";
		next; # don't waste computer time repeating work already done
	}
=cut
	my $title = $gene_key->{$query} // $query;
#	my $title = $deg2gene_name->{$query}{'Gene Name'} // $query;
#	my $key = $title;
#	$title .= " ($db->{$title}{Description}[0])" if defined $db->{$title}{Description}[0];
#	my $deg_number = $deg_id;
#	$deg_number =~ s/^DEG//;

	say '@ ' . __FILE__ . '& line' . __LINE__;
	my $files = msa_phylo_plot({
		'blast.json'          	=> $all_blast,	# {$species}{$query}
		db								=> \%all_fasta,
		debug						 	=> 0,
#		'msa.image.filename'  	=> $msa_image,
#		'msa.ylabel'          	=> 'Species',
		'out.fasta'           	=> $out_fasta,
#		'phylo.image.filename'	=> $phyl_image,
#		'phylo.image.width'   	=> 13,
		query                 	=> $query,
		'query.fasta'         	=> $all_query_seq->{$query},
		'query.species'		 	=> $files{query_species},
		'run.clustal'			 	=> 1,
		'species.order'       	=> [$query, @order],
		'species.only'        	=> 1,
		'tex.outfile'         	=> $tex_outfile,
		title                 	=> "$title ($group)",
	});
	say '@ ' . __FILE__ . '& line' . __LINE__;
}
