#!/usr/bin/env perl

use strict;
use feature 'say';
use File::Path 'rmtree';
use warnings FATAL => 'all';
use autodie ':default';
use Getopt::Long 'GetOptions';
use FindBin 1.51 '$RealBin';
use lib $RealBin;
use DDP;
use Devel::Confess 'color';
use Capture::Tiny 'capture';

sub execute {
	my $cmd = shift;
	my ($stdout, $stderr, $exit) = capture {
		system( $cmd );
	};
	if ($exit != 0) {
		say "exit = $exit";
		say "STDOUT = $stdout";
		say "STDERR = $stderr";
		die "$cmd failed"
	}
	return [$stdout, $stderr]
}
my %help = (
	'hosts'				=> 'The tab-delimited file of species name and source file for hosts',
	'output-svg'		=> 'Output scatterplot in SVG format',
	'pathogens'			=> 'The tab-delimited file of species name and source file for pathogens',
	'plot-title'		=> 'Title on plot for resulting scatterplot',
	'query-list'		=> 'The fasta file of protein/RNA sequences',
	'query-species'	=> 'Source species for gene/RNA-list',
	'test'				=> 'Run test to ensure that ProGen works; every other option is ignored'
);

my %args;

GetOptions(
	'help|h'         	=> \$args{help},
	'hosts=s'      	=> \$args{hosts},
	'output-svg=s'		=> \$args{'output-svg'},
	'pathogens|p=s'	=> \$args{pathogens},
	'plot-title=s'		=> \$args{'plot-title'},
	'query-list|g=s'	=> \$args{'query-list-file'},
	'query-species=s'	=> \$args{'query-species'},
	'test'				=> \$args{test}
) or die 'error with GetOptions';

if (defined $args{help}) {
	p %help;
	exit;
}
if (defined $args{test}) {
	$args{pathogens}			= "$RealBin/../test/pathogens.tsv";
	$args{hosts}		 		= "$RealBin/../test/hosts.tsv";
	$args{'query-list-file'}= "$RealBin/../test/test.progen.fa";
	$args{'output-svg'}	  	= 'scatterplot.test.svg';
	$args{'plot-title'}	  	= 'Title';
}
my @undef_options = grep {not defined $args{$_}} ('query-list-file', 'output-svg', 'hosts', 'pathogens', 'plot-title');

if (scalar @undef_options > 0) {
	p @undef_options;
	say STDERR 'the above options were not defined at the command line, and must be defined.';
	p %help;
	die;
}
foreach my $dir (grep {-d $_} ('blast.results','fa','svg','tex', 'xlsx', 'json')) {
	rmtree $dir;
}
execute("perl $RealBin/0.blast.pl --hosts $args{hosts} --pathogens $args{pathogens} --query-list $args{'query-list-file'}");
execute("perl $RealBin/1.make.input.pl --hosts $args{hosts} --pathogens $args{pathogens} --blast-json json/blast.hits.json --group Test --query-list $args{'query-list-file'}");
execute("perl $RealBin/2.find.regions.pl --hosts $args{hosts} --pathogens $args{pathogens}");
execute("perl $RealBin/9.plot.protein.targetable.length.pl --title $args{'plot-title'} --output-svg $args{'output-svg'} --target-json 'json/target.regions.json'");
