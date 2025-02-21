#!/usr/bin/env perl

use 5.038;
use warnings FATAL => 'all';
use autodie ':default';
use Util qw(json_file_to_ref);# ref_to_json_file);list_regex_files
use latex qw(write_2d_array_to_tex_tabular);

foreach my $d ('xlsx','tex') {
	mkdir $d unless -d $d;
}
my @header = ('DEG key','Region','Full Gene Name','length','# pathogen species','mean no. of pathogens', 'min(Similarity)','mean(Pathogen Similarity)','max(Similarity)','Start','End','Domains');
#----------------
sub make_table (
		$col,				# which column should be the sort key?
		$cols,			# which columns should appear in the final table?
		$tex_filename,	# the output tex filename
		$json
	) {
	my $xls = json_file_to_ref( $json );
	@{ $xls } = grep { $_->[3] >= 15} @{ $xls }; # eliminate very short regions
	foreach my $row (@{ $xls }) {
		$row->[1] .= ".$row->[9]-$row->[10]"; # re-name to show that I'm talking about regions, not genes
	}
	my @data = sort {$b->[$col] <=> $a->[$col]} @{ $xls }; # sort by decreasing length of the region
	foreach my $row (@data) {
		@{ $row } = @{ $row }[@{ $cols }]; # only use some rows
	}
	write_2d_array_to_tex_tabular({
	#	xlsx_filename  => 'xlsx/short.table.0.xlsx',
		data           => \@data,
		format			=> 1,
		tex_filename   => $tex_filename,
		header         => [@header[@{ $cols }]],
		'max.rows.tex' => 25,
		'max.width'    => 90
	});
}
=header
[
    [0]  "DEG key",
    [1]  "Gene Name",
    [2]  "Full Gene Name",
    [3]  "length",
    [4]  "# pathogen species",
    [5]  "mean # of pathogens",
    [6]  "min(Similarity)",
    [7]  "sum(Pathogen Similarity)",
    [8]  "max(Similarity)",
    [9]  "Start",
    [10] "End",
    [11] "Domains"
] (757B)
=cut
foreach my $d ('who.critical', 'top10.fungal') {
	make_table(
		3,		# sort by this key
		[1,2,3,5,7],
		"tex/short.table.$d.tex",
		"/home/con/identify.target/$d/json/all.fungi.target.regions.json"
	);
}
