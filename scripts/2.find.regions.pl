#!/usr/bin/env perl

use strict;
use feature 'say';
use FindBin 1.51 qw( $RealBin );
use lib $RealBin;
use warnings FATAL => 'all';
use autodie ':default';
use util qw(write_2d_array_to_tex_tabular fasta2hash list_regex_files execute json_file_to_ref ref_to_json_file plots_in_rows);
use File::Temp 'tempfile';
use List::Util qw(min max sum);
use Excel::Writer::XLSX;
use Cwd 'getcwd';
use Getopt::Long 'GetOptions';
use DDP;

# example "perl ../pl/2.find.regions.pl --hosts hosts.tsv --pathogens pathogens.tsv"

my ($help, %files);
GetOptions(
#	'blast-json=s'  => \$files{blast_json},
#	'msa-program=s' => \$msa_program,
#	'gene-list|g=s' => \$files{gene_list_file},
	'help'          => \$help,
	'hosts=s'       => \$files{hosts},
#	'output-stem'	 => \$stem,
	'pathogens|p=s' => \$files{pathogens},
	'query-species|q=s'	=> \$files{query_species},
) or die 'error with GetOptions';

unless ((defined $files{hosts}) && (defined $files{pathogens})) {
	die 'the correct syntax is: "perl ../pl/2.find.regions.pl --hosts hosts.tsv --pathogens pathogens.tsv"';
}
my %group_species;

sub read_tsv {
	my $tsv_file = shift;
	my $group    = shift;
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

read_tsv( $files{hosts},     'hosts' );	  # writes to %group_species
read_tsv( $files{pathogens}, 'pathogens' ); # writes to %group_species

my $database = json_file_to_ref( "$RealBin/sneath.phi.distance.json" );

sub similarity_score {
=purpose
	this subroutine takes an MSA of proteins as input and outputs a similarity score based on the database selected (e.g. Sneath)
	
	example:
	-------
	similarity_score({
		'1d.heatmap.file'		=> "svg/sneath.$p.colobar.heatmap.svg", # output
		fasta						=> "../who.critical/fa/clustal.DEG2001$p.output.aln", # input
		database					=> $database,
		groups					=> { # divide the species within the MSA to groups
			pathogens => ['A.fumigatus','C.albicans','C.auris','C.tropicalis','C.parapsilosis',
	'C.neoformans.B.3501A','C.neoformans.grubii.H99','C.neoformans.JEC21','H.capsulatum','N.glabratus'],
			hosts		 => ['O.sativa', 'S.tuberosum', 'Z.mays', 'G.max', 'H.sapiens']
		},
		diff						=> ['pathogens', 'hosts', "json/DEG2001$p.sneath.diff.json"], # get a difference graph for pathogens - hosts
		figwidth					=> 12,
		'multilineplot.file' => "$d/svg/sneath.$p.svg",
		plots_in_rows_image	=> "$d/svg/sneath_plots_in_rows_$p.svg",
		title					 	=> "$title Amino Acid Similarity",
		xlabel					=> 'MSA amino acid index',
		ylabel					=> 'Sneath Similarity',
	}
=cut
	my ($args) = @_;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	unless (ref $args eq 'HASH') {
		die "args must be given as a hash ref, e.g. \"$current_sub({ data => \@blah })\"";
	}
	my @reqd_args = (
		'fasta',	# fasta file/hash, where each key is a *string* of amino acids, NOT an array
		'database',	# e.g. sneath similarity. Should be hash or JSON file
	);
	my @undef_args = grep {!defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'The arguments above are necessary for proper function and weren\'t defined.';
	}
	my @defined_args = ( @reqd_args,
		'1d.heatmap.file',	# output file to 1D heatmap to view, with rows if multiple groups
#		'all',					# make an "all" group for everything; by default off
		'domains',				# Interpro domains, add as a 3rd plot
		'diff',					# get a difference graph for groups in "groups", which are specified here.
# example "diff": diff => ['pathogens', 'hosts'], which will produce pathogens - hosts; 3rd item is file to output data to
		'groups',				# split the data into groups. Any key named "all" will have all of the keys
		'figwidth',				# default 4.8
		'multilineplot.file',# output file to multiline_plot, e.g. "file.svg"
		'plots_in_rows_image',# plot plots in rows, e.g. "file.svg"
		'scatterplot.file',	# output file to scatterplot to view
		'title',					# if "scatterplot.file" is defined, use a title
		'xlabel', 'ylabel',	# strings
		'ylim',					# ylim => [0,1]
	);
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args, array_max => scalar @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args, array_max => scalar @defined_args;
		die "The above args are accepted for $current_sub";
	}
	my ($d, %unique_lengths, $database, $n_aa, %array_form);
	if (ref $args->{fasta} eq 'HASH') {
		$d = $args->{fasta};
	} elsif (-f -r $args->{fasta}) {
		$d = fasta2hash( $args->{fasta} );
	} else {
		die "$args->{fasta} is neither a file nor a hash";
	}
	my %lengths = map { $_ => length $d->{$_}} keys %{ $d };
	foreach my $key (keys %lengths) {
		$unique_lengths{$lengths{$key}}++;
	}
	if (scalar keys %unique_lengths != 1) {
		p %lengths;
		die 'all lengths must be identical';
	}
	if (ref $args->{database} eq 'HASH') {
		$database = delete $args->{database}; # don't read in the file again
	} elsif (-f -r $args->{database}) {
		$database = json_file_to_ref($args->{database});
	} else {
		die "$args->{database} broke pigeonholes.";
	}
	if ((defined $args->{diff}) && (not defined $args->{groups})) {
		p $args;
		die '"diff" was defined, but "groups" was not defined.';
	}
	if ((defined $args->{diff}) && (ref $args->{diff} ne 'ARRAY')) {
		p $args->{diff};
		die '"diff" must be an array reference.';
	}
	foreach my $key (keys %{ $d }) {
		@{ $array_form{$key} } = split '', $d->{$key};
		$n_aa = scalar @{ $array_form{$key} } - 1;
	}
	my (%group_similarity, %groups);
	if (defined $args->{'groups'}) {
		%groups = %{ $args->{groups} };
	}
#	$groups{all} = $groups{all} // 0;
#	if ($groups{all} > 0) {
#		$groups{all} = [keys %{ $d }];
#	}
	foreach my $group (keys %groups) {
		foreach my $i (0..$n_aa) {
			$group_similarity{$group}[$i] = 0; # initialize
			my $n_seq = 0;
			foreach my $pr1 (@{ $groups{$group} }) { # protein1 = $pr1
				next if not defined $array_form{$pr1}[$i];
				next if $array_form{$pr1}[$i] eq '-';
				if (not defined $database->{$array_form{$pr1}[$i]}) {
					p $args;
					die "\$database->{$array_form{$pr1}[$i]} isn't defined in $group @ position $i";
				}
				foreach my $pr2 (grep {$_ ne $pr1} @{ $groups{$group} }) {
					next if not defined $array_form{$pr2}[$i];
					next if $array_form{$pr2}[$i] eq '-';
					if (not defined $database->{$array_form{$pr1}[$i]}) {
						p $args;
						die "\$database->{$array_form{$pr1}[$i]} isn't defined in $group @ position $i";
					}
					$group_similarity{$group}[$i] += $database->{$array_form{$pr1}[$i]}{$array_form{$pr2}[$i]};
					$n_seq++;
				}
			}
			next if $n_seq <= 1; # cannot divide by 0
			$group_similarity{$group}[$i] /= ((scalar @{ $groups{$group}} ) * ((scalar @{ $groups{$group}}) - 1));
			if (($group_similarity{$group}[$i] > 1) || ($group_similarity{$group}[$i] < 0)) {
				p $groups{$group};
				p $group_similarity{$group}, array_max => scalar @{ $group_similarity{$group} } + 1;
				die "position $i has a value outside of [0,1], which shouldn't be possible";
			}
		}
	}
	if (defined $args->{diff}) {
		my @diff = map { $group_similarity{$args->{diff}[0]}[$_] - $group_similarity{$args->{diff}[1]}[$_] } 0..$n_aa;
		ref_to_json_file(\@diff, $args->{diff}[2]);
	}
# do not use any of these plotting options if there are other types of figures in rows, it makes the scripts too complicated
	if (defined $args->{'scatterplot.file'}) {
		my %args;
		foreach my $group (sort keys %group_similarity) {
			$args{data}{$group}{x} = [0..$n_aa];
			$args{data}{$group}{y} = [@{ $group_similarity{$group} }];
		}
		$args{filename} = $args->{'scatterplot.file'};
		foreach my $arg (grep {defined $args->{$_}} ('title')) {
			$args{$arg} = $args->{$arg};
		}
		scatterplot({ %args });
	}
	if (defined $args->{'multilineplot.file'}) {
		my %args;
		foreach my $group (sort keys %group_similarity) {
			$args{data}{$group}{x} = [0..$n_aa];
			$args{data}{$group}{y} = [@{ $group_similarity{$group} }];
		}
		foreach my $arg (grep {defined $args->{$_}} ('figheight', 'figwidth', 'title', 'xlabel', 'ylabel')) {
			$args{$arg} = $args->{$arg};
		}
		$args{filename} = $args->{'multilineplot.file'};
		$args{xlim} = [0,$n_aa];
		multiline_plot({ %args });
	}
	if (defined $args->{'1d.heatmap.file'}) {
		my %args;
		foreach my $group (sort keys %group_similarity) {
			$args{data}{$group}[0] = [0..$n_aa];
			$args{data}{$group}[1] = [@{ $group_similarity{$group} }];
		}
		foreach my $arg (grep {defined $args->{$_}} ('figheight', 'figwidth', 'title', 'xlabel')) {
			$args{$arg} = $args->{$arg};
		}
		$args{filename} = $args->{'1d.heatmap.file'};
		heatmaps_in_rows({ %args });
	}
	if (defined $args->{plots_in_rows_image}) {
		my %args;
		foreach my $group (sort keys %group_similarity) {
			$args{data}{$group}[0] = [0..$n_aa];
			$args{data}{$group}[1] = [@{ $group_similarity{$group} }];
		}
		foreach my $arg (grep {defined $args->{$_}} ('figheight', 'figwidth', 'title', 'xlabel', 'ylabel', 'ylim')) {
			$args{$arg} = $args->{$arg};
		}
		$args{filename} = $args->{plots_in_rows_image};
		$args{xlim} = [0,$n_aa];
		$args{ylim} = [0,1];
		plots_in_rows({ %args });
	}
	return \%group_similarity;
}
sub get_pathogen_only_regions {
=purpose
	take an MSA as input and report which pathogen regions are not covered by any host.
	
	arguments:
	
	Required:
	--------
	msa_file:
		MSA alignment in fasta format

	hosts:
		an array of strings, which is species

	pathogens:
		an array of strings, which is species

	Optional:
	--------
	1d heatmap:

	annotation:
		InterPro Annotation file, which will label regions

	gene_name:
		gene name for XLSX output
	msa_key:
		the key from which the MSA will be 
	
	regions.json:
		a JSON of target regions, which is the result of this subroutine
=cut
	my ($args) = @_;
	my @reqd_args = ('msa_file', 'hosts', 'pathogens');
	my @undef_args = grep { !defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined.';
	}
	if (scalar @{ $args->{pathogens} } == 0) {
		p $args;
		die 'Pathogens have no species listed, there is nothing to do';
	}
	my @defined_args = (@reqd_args,
#	'1d heatmap',
	'annotation',
	'deg_key',
	'figwidth',				# for the plots in rows figure (msa.similarity.image), default 6.4
	'gene_name',
	'json.output',			# write the line segments and InterPro names to a JSON for later plotting
	'long.name', 'msa.key',
	'plots_in_rows_image','output.image', 'pathogen.subset', 'regions.json', 'title');
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args;
		die 'The above args are accepted.'
	}
	if (grep {not defined $args->{$_}} keys %{ $args }) {
		p $args;
		die "undefined arguments were given.";
	}
	my $fasta = fasta2hash($args->{msa_file});
	if (not defined $fasta->{$args->{'msa.key'}}) {
		p $fasta;
		die "\"$args->{'msa.key'}\" isn't defined in $args->{msa_file}";
	}
	my @original_seq = split '', $fasta->{$args->{'msa.key'}};
	my @pathogen_species = grep { defined $fasta->{$_}} @{ $args->{pathogens} };
#	say colored(['cyan on_bright_red'], scalar @pathogen_species);
	if (scalar @pathogen_species == 0) {
		return []; #there's no point in analyzing this without any pathogens
	}
	my @host_species =     grep { defined $fasta->{$_}} @{ $args->{hosts} };
	foreach my $species_to_skip (keys %{ $fasta }) {
		if (
				(not grep {$_ eq $species_to_skip} @host_species)
				&&
				(not grep {$_ eq $species_to_skip} @pathogen_species)
			) {
			delete $fasta->{$species_to_skip};
		}
	}
	# do all sequences have the same length? bug check
	my %length = map { $_ => length $fasta->{$_}} keys %{ $fasta };
	my (%uniq_length, $length);
	foreach my $l (values %length) {
		$uniq_length{$l}++;
		$length = $l - 1;
	}
	if (scalar keys %uniq_length != 1) {
		p $args;
		say STDERR "from $args->{msa_file}, the following strings:";
		p %length;
		say STDERR 'which have the following lengths:';
		p %uniq_length;
		die 'all keys should have the same length';
	}
	undef %length;
	my %data; # split the strings into arrays
	foreach my $species (keys %{ $fasta } ) {
		@{ $data{$species} } = split '', $fasta->{$species};
	}
	# only species that are defined in the alignment file will show up as pathogens, don't worry if other pathogens are absent
#	say colored(['red on_black'], 'pathogens:');
#	p @pathogen_species;
	my @target_segments;
	foreach my $pos (0..$length) { # which positions have only pathogen
		my (@pathogens, @hosts);
		foreach my $host (@host_species) {
			push @hosts, $data{$host}[$pos];
			last if $data{$host}[$pos] ne '-';
		}
		next if grep {$_ ne '-'} @hosts; # skip this a. acid if any hosts have this amino acid aligned
		foreach my $pathogen (@pathogen_species) {
			if (not defined $data{$pathogen}[$pos]) {
				p $data{$pathogen};
				p $args;
				die "\$data{$pathogen}[$pos] isn't defined.";
			}
			next if $data{$pathogen}[$pos] eq '-';
			push @pathogens, $data{$pathogen}[$pos];
		}
		next if scalar @pathogens == 0;
		if (not defined $target_segments[0][0]) { # the first amino acid
			$target_segments[0][0] = $pos;
			$target_segments[0][1] = $pos;
			push @{$target_segments[0][2] }, scalar @pathogens; # how many pathogens are in this particular segment?
			next;
		}
		if ($pos == ($target_segments[-1][1] + 1)) {
			$target_segments[-1][1] = $pos;	# extend current segment
			push @{ $target_segments[-1][2] }, scalar @pathogens;
			next;
		}
		if ($pos > ($target_segments[-1][1]+1)) {
			push @target_segments, [$pos,$pos, [scalar @pathogens]]; # start new segment
		}
	}
	if (scalar @target_segments == 0) {
		return [];
	}
	while (my ($index, $segment) = each @target_segments) {
		if ($segment =~ m/\*/) {
			splice @target_segments, $index, 1;
			say "$segment has stop codon, removing.";
			next;
		}
		if (($segment->[1] - $segment->[0]) < 15) {
			splice @target_segments, $index, 1; # remove segments that are too short
			next;
		}
		# if there are < 15 amino acids in the segment, delete this segment
		if ((scalar grep { /[A-Ya-y]/ } @original_seq[$segment->[0]..$segment->[1]]) < 15) {
			say 'skipping segment, too few amino acids';
			splice @target_segments, $index, 1;
		}
	}
#	$args->{debug} = $args->{debug} // 0; # run by default
	my %args = (
		fasta		=> $fasta,
		database	=> $database,
		groups	=> {
			hosts			=> $args->{hosts},
			pathogens	=> $args->{pathogens},
		},
		xlabel	=> 'MSA amino acid index',
		ylabel	=> 'Sneath Similarity',
		ylim		=> [-0.05, 1.05],
	);
	foreach my $arg (grep {defined $args->{$_}} ('figwidth', 'plots_in_rows_image', 'title')) {
		$args{$arg} = $args->{$arg};
	}
	my $similarity = similarity_score({ %args });
=if (defined $args->{'1d heatmap'}) {
# https://stackoverflow.com/questions/45841786/creating-a-1d-heat-map-from-a-line-graph
		my $unlink = 0;
		my ($py, $temp_py) = tempfile(DIR => '/tmp', SUFFIX => '.py', UNLINK => $unlink);
		say $temp_py if $unlink == 0;
		say $py 'import matplotlib.pyplot as plt';
		say $py 'import numpy as np';
		say $py 'plt.rcParams["figure.figsize"] = 5,2';
		say $py 'x = np.arange(0, ' . scalar @similarity_score . ')';
		say $py 'y = np.array([' . join (',', @similarity_score) . '])';
		say $py 'fig, ax = plt.subplots()';
#		say $py 'extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]';
		say $py 'im = ax.imshow(y[np.newaxis,:], cmap="gist_rainbow", aspect="auto")';#, extent=extent)';
		say $py 'ax.set_yticks([])';
#		say $py 'ax.set_xlim(extent[0], extent[1])';
#		say $py 'ax2.plot(x,y)';
		
		say $py 'cbar = ax.figure.colorbar(im, ax = ax, label = "Similarity")';
		say $py 'plt.xlabel("Amino Acid Residue (MSA Coordinate)")';
		say $py 'fig.tight_layout(pad = 1.0)';
		say $py "plt.savefig('$args->{'1d heatmap'}')";
		close $py;
		execute("python3 $temp_py");
		say colored(['cyan on_bright_yellow'], "wrote $args->{'1d heatmap'}");
	}
=cut
	my (@target_regions, @domains);
# compare the similarities
	if (defined $args->{annotation}) {# && (defined $args->{'output.image'})) {
		my %args = (
			'msa.input'			=> $args->{msa_file},
			'msa.key'         => $args->{'msa.key'},
			'interpro.output' => $args->{annotation},
			'shaded.regions'  => \@target_segments,
			title             => "$args->{title}: $args->{'long.name'}"
		);
		foreach my $arg (grep {defined $args->{$_}} ('json.output', 'output.image')) {
			$args{$arg} = $args->{$arg};
		}
	}
#---------------
# how many host & pathogen species are present in the multiple alignment
#---------------
	my %species;
	open my $aln_fh, '<', $args->{msa_file};
	while (<$aln_fh>) {
		if (/^>(.+)/) {
			$species{$1} = 1;
		}
	}
	close $aln_fh;
#	my $n_hosts     = scalar grep {defined $species{$_}} @{ $args->{hosts}     };
#	my $n_pathogens = scalar grep {defined $species{$_}} @{ $args->{pathogens} };
	foreach my $target_region (@target_segments) {
		# make a row for the output array, which will become a row in excel
		next if ($target_region->[1] - $target_region->[0]) == 0;
		my $length = $target_region->[1] - $target_region->[0] + 1;
		my $mean = sum( @{ $similarity->{pathogens}}[$target_region->[0]..$target_region->[1]]) / $length;
		die "$mean > 1" if $mean > 1;
		push @target_regions, {
			deg_key						=> $args->{deg_key},
			gene_name					=> $args->{'gene_name'},															# gene name
			long_name 					=> $args->{'long.name'},															# long name
			'length'						=> $target_region->[1] - $target_region->[0],								# length
#			$n_pathogens,
			'mean # of pathogens'	=> sum( @{ $target_region->[2] }) / $length,# scalar( @{ $target_region->[2] }),# mean # of pathogens at positions
			min_similarity				=> min( @{ $similarity->{pathogens}}[$target_region->[0]..$target_region->[1]]),
			mean_similarity			=> $mean,
			max_similarity				=> max( @{ $similarity->{pathogens}}[$target_region->[0]..$target_region->[1]]),
			start							=> $target_region->[0],																# start
			end							=> $target_region->[1],																# end
		};
		if (grep {not defined $target_regions[-1]{$_}} keys %{ $target_regions[-1] }) {
			p $target_regions[-1];
			die 'keys were undefined';
		}
	}
	if (defined $args->{'regions.json'}) {
		ref_to_json_file(\@target_regions, $args->{'regions.json'});
	}
	return \@target_regions;
}
my $deg_annotation = json_file_to_ref("$ENV{HOME}/identify.target/deg/all.deg.info.json");
foreach my $d ('svg/regions', 'svg/similarity', 'json/regions', 'xlsx', 'json') {
	mkdir $d unless -d $d;
}
my $cwd = getcwd();
my @all_fungi;#, %individual_pathogens_xlsx);
#my $db = json_file_to_ref('yeast.gene.descr.other.aliases.json');
my @aln_files = list_regex_files('\.output\.aln$', "$cwd/fa");
if (scalar @aln_files == 0) {
	die "no alignment files were found in \"$cwd/fa\"";
}
foreach my $fa (@aln_files) {
	my $aln = $fa;
	my $name;
	if ($fa =~ m/clustal\.(.+)\.output.aln$/) {
		$name = $1;
	} else {
		die "$fa failed regex.";
	}
	my $key;
	if ($name =~ m/^DEG(\d+)$/) {
		$key = $1;
	} else {
		die "$name failed regex";
	}
	my $gene_name = $key;
	if (defined $deg_annotation->{$key}{'Gene name'}) {
		$gene_name = $deg_annotation->{$key}{'Gene name'} unless $deg_annotation->{$key}{'Gene name'} eq '-';
		say "got $gene_name";
#	} else {
#		p $deg_annotation->{$key};
#		die "Couldn't get a gene name for $key";
	}
	my $title = $key;
	if (defined $deg_annotation->{$key}{Description}) {
		$title = $deg_annotation->{$key}{Description};
	}
#	die "$aln doesn't have a defined gene name" if not defined $gene_name;
#--------
	my $deg_results = get_pathogen_only_regions({
#		'1d heatmap'	=> "svg/similarity/$fa" . '_all_1d_heatmap.svg',
#		annotation		=> "../eyeball/json/DEG$fa.json", # interpro annotation
		deg_key			=> $name,
		figwidth			=> 12,
		gene_name		=> $gene_name,
#		'json.output'	=> "$cwd/json/regions/interpro.domains.$deg_key.json",
		'long.name'		=> $title,
		msa_file			=> $aln,
		'msa.key'		=> $files{query_species},#$deg_key,
#		'plots_in_rows_image' => "svg/similarity/$fa" . '_all_pathogens.svg',
		hosts		  		=> $group_species{hosts},
		pathogens		=> $group_species{pathogens},
		'regions.json'	=> "$cwd/json/regions/$name.regions.json",
#		'output.image'	=> "$cwd/svg/regions/all_fungi_$deg_key.svg",
		title				=> $title,
	});
	foreach my $result (@{ $deg_results }) {
		push @all_fungi, $result;
	}
=foreach my $pathogen (@{ $group_species{pathogens} }) {
		my $pathogen_name = $pathogen;
		$pathogen_name =~ s/\./_/g;
		my $individual_pathogen = get_pathogen_only_regions({
#			'1d heatmap'	=> "svg/similarity/$fa" . "_$pathogen" . '_1d_heatmap.svg',
			annotation		=> "../eyeball/json/DEG$fa.json", # interpro annotation
#			figwidth			=> 12,
			gene_name		=> $gene_name,
			'long.name'		=> $title,
			msa_file			=> $aln,
			'msa.key'		=> $deg_key,
#			'msa.similarity.image' => "svg/similarity/$fa" . "_$pathogen_name.svg",
			hosts		  		=> $group_species{hosts},
			pathogens		=> [$pathogen],
#			'output.image'	=> "$cwd/svg/regions/$pathogen_name" . "_$deg_key.svg",
			'regions.json'	=> "$cwd/json/regions/$pathogen.json",
			title				=> "$deg_key: $gene_name",
		});
		foreach my $row (@{ $individual_pathogen }) {
			push @{ $individual_pathogens_xlsx{$pathogen} }, [$deg_key, @{ $row } ];
		}
	}
=cut
}
my @header = ('DEG key', 'Gene Name', 'Full Gene Name', 'length', 
#'# pathogen species',
'mean # of pathogens', 'min(Similarity)', 'mean(Pathogen Similarity)', 'max(Similarity)', 'Start', 'End');
my @results;
foreach my $region (@all_fungi) {
	push @results, [@{ $region }{'deg_key', 'gene_name', 'long_name', 'length', 'mean # of pathogens', 'min_similarity', 'mean_similarity', 'max_similarity', 'start', 'end'}];
}
if (scalar @results == 0) {
	p @results;
	die 'No results were found';
}
ref_to_json_file(\@all_fungi, 'json/target.regions.json');
write_2d_array_to_tex_tabular({
	xlsx_filename  => 'xlsx/target.regions.xlsx',
	data           => \@results,
	tex_filename   => 'tex/target.regions.tex',
	header         => \@header,
	'max.rows.tex' => 20,
	'max.width'    => 90,
});
=foreach my $pathogen (sort keys %individual_pathogens_xlsx) {
	write_2d_array_to_tex_tabular({
		xlsx_filename  => "xlsx/$pathogen.target.regions.xlsx",
		data           => $individual_pathogens_xlsx{$pathogen},
		tex_filename   => "tex/$pathogen.target.regions.tex",
		header         => \@header,
		'json.output'  => "json/$pathogen.target.regions.json",
		'max.rows.tex' => 20,
		'max.width'    => 90,
	});
}
