#!/usr/bin/env perl

use strict;
use feature 'say';
use File::Temp 'tempfile';
use Scalar::Util 'looks_like_number';
use FindBin 1.51 qw( $RealBin );
use lib $RealBin;
use DDP;
use warnings FATAL => 'all';
use autodie ':default';
use util 'json_file_to_ref';
use Getopt::Long 'GetOptions';
use Capture::Tiny 'capture';
use Cwd 'getcwd';
use Term::ANSIColor;

my (%args, $help);
GetOptions(
#	'blast-json=s'  => \$files{blast_json},
#	'msa-program=s' => \$msa_program,
#	'gene-list|g=s' => \$files{gene_list_file},
	'help'         => \$help,
	'target-json=s'=> \$args{target_json},
	'orthologs=s'	=> \$args{exclude}, # list of genes to remove
	'gene-labels=s'=> \$args{gene_labels},
	'output-svg=s'	=> \$args{output_svg},
	'group=s'		=> \$args{group},
	'title=s'		=> \$args{title}
) or die 'error with GetOptions';

my @undef_args = grep {!defined $args{$_}} ('target_json', 'output_svg');
if (scalar @undef_args > 0) {
	p @undef_args;
	die 'The above args need to be defined.';
}

sub text { # this is a helper function to other matplotlib subroutines
=example
	this is called *inside* other subroutines to shorten code
	"text" should be a 2D array wherein each element looks like:
	
	[x coord, y coord, "txt"]
	OR
	[x coord, y coord, text, xtext, ytext] if text location is specified for arrows
=cut
	my ($args) = @_;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	unless (ref $args eq 'HASH') {
		die "args must be given as a hash ref, e.g. \"$current_sub({ data => \@blah })\"";
	}
	my @reqd_args = (
		'text',	# a 2D array
		'fh',		# e.g. $py, $fh, which will be passed by the subroutine
		'object'	# e.g. "plt"
	);
	my @undef_args = grep { !defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined.';
	}
	my @defined_args = ( @reqd_args,
		'annotate',		# label using plt.annotate instead; 0 or 1; cf. text
		'adjust_text',	# use the adjust_text https://adjusttext.readthedocs.io/en/latest/Examples.html
		'xshift'	# if using "annotate" option; how far to shift right/left as a % of x axis
	);
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args, array_max => scalar @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args, array_max => scalar @defined_args;
		die 'The above args are accepted.'
	}
	$args->{adjust_text} = $args->{adjust_text} // 0;
	$args->{annotate}		= $args->{annotate}    // 0;
	if ($args->{adjust_text} > 0) { # can't handle bold font
# make the strings all 1 array
#		my $str = 'txt = [(' . join ('),(', @{ $args->{text} }) . ')]';
		say {$args->{fh}} 'data = [';
		foreach my $text (@{ $args->{text} }) {
			my @t = @{ $text }[0..2];
			$t[2] = "'$t[2]'";
			say {$args->{fh}} '	(' . join (',', @t) . '),';
		}
		say {$args->{fh}} ']';
		say {$args->{fh}} 'from adjustText import adjust_text';
		say {$args->{fh}} 'texts = [plt.text(x, y, l) for x, y, l in data]';
		say {$args->{fh}} 'adjust_text(texts)';
	} elsif ($args->{annotate} > 0) { # use lines to point from text label
		if ((defined $args->{xshift}) && (($args->{xshift} < 0) || ($args->{xshift} > 100))) {
			p $args;
			die 'xshift must be within [0,100]';
		}
		say {$args->{fh}} 'xmin, xmax, ymin, ymax = plt.axis()';
		$args->{xshift} = $args->{xshift} // 15;
		say {$args->{fh}} 'xshift = (xmax - xmin) * 0.01 * ' . $args->{xshift};
		foreach my $text (@{ $args->{text} }) {
			my $xytext_coords;
			if (scalar @{ $text } == 5) { # xy coords of point, then xy coords of text
				$xytext_coords = join (',', @{ $text }[3,4]);
			} else { #3f (scalar @{ $text } == 3) { # no defined shift location
				$xytext_coords = "$text->[0] + xshift, $text->[1]";#join (',', @{ $text }[0,1]);
#			} else {
#				p $text;
#				die "at text index $i (above), there can only be 3 or 5 elements.";
			}
			say {$args->{fh}} "$args->{object}.annotate(text = '$text->[2]', xy = [$text->[0],$text->[1]], xytext = [$xytext_coords], arrowprops = dict(arrowstyle = 'simple', facecolor = 'red', mutation_scale = 0.5))";
		}
	} else {
		foreach my $text (@{ $args->{text} }) {
			my @t = @{ $text };
			$t[2] = "'$t[2]'";
			say {$args->{fh}} "$args->{object}.text(" . join (',', @t) . ')';
		}
	}
}

sub scatterplot_2d_color {
=subroutine info

example:
---------
scatterplot_2d_color({
	filename	=> $svg_name,
	title		=> $title,
	data		=> {
		'Alignment length (a.a.)'	=> [@length],
		'Mean # of Pathogens'	=> [@mean_pathogens],
		'Mean Similarity'			=> [@mean_similarity]
	},
	color_key=>	'Mean Similarity',
	plot_params => {
		alpha => 0.5
	}
});
	----------------
	Req'd options:
	----------------
	data:
input data should look like this:
{
    align_len   [
        [0] 59,
        [1] 59,
            (...skipping 16 items...)
    ],
    evalue      [
        [0] 0.291346,
        [1] 0.33064,
            (...skipping 16 items...)
    ],
    score       [
        [0] 81,
        [1] 81,
            (...skipping 16 items...)
    ]
} (3K)
	the 3rd key (sorted case-insensitive alphabetically) will be the color key if "color_key" isn't specified
	the xlabel, ylabel, and color label are the key names themselves.  "xlabel", "ylabel", etc. aren't accepted
	another example:
		scatterplot_2d_color({
#		adjust_text	=> 1,
		annotate		=> 1, # use the arrow annotation
#		function => { # plot this function within the scatterplot
#			min	=> 0,
#			max	=> $xint,
#			func	=> "-(10.0/($xint*$xint))*(x+$xint)*(x-$xint)",
#			step	=> 1
#		},
		cb_min		=> 0,
		cb_max		=> 1,
		filename		=> $args->{svg_name},
		title			=> $args->{title},
		data			=> {
			'Alignment length (a.a.)'			=> [@length],
#			'Total Protein Length (a.a.)'		=> [@full_protein_length],
			'Mean Pathogen Column Occupancy'	=> [@mean_pathogens],
			'Mean Pathogen ϕ'						=> [@mean_similarity]
		},
		color_key	=>	'Mean Pathogen ϕ',#'Total Protein Length (a.a.)',
		plot_params => {
			alpha => 0.5
		},
		text			=> [@text],
		xmin			=> 0,
#		xlim			=> [0, 2290],
		ylim			=> [1,10.2]
	});
	where @text is [x, y, text, xtext, ytext] for arrowed annotation
=cut
	my ($args) = @_;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	unless (ref $args eq 'HASH') {
		die "args must be given as a hash ref, e.g. \"$current_sub({ data => \@blah })\"";
	}
	my @undef_args = grep { !defined $args->{$_}} ('data', 'filename');
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined.';
	}
	my @defined_args = (
	'annotate',		# use annotate function to use arrows for scatterplots
	'adjust_text',	# adjust the text with the "text" helper subroutine.  Set to 1 for "on"
	'cb_min','cb_max', # minimum/maximum color values
	'color_key', 'custom_legend', 'data', 'figwidth', 'filename',
	'function', # hash which takes min, max, func, and step (using "function" helper subroutine)
	'grid',
	'keys', # specify the order, otherwise alphabetical
	'line_segment',
#example
#	line_segment	=> [{
#		'x0'	=> $min,	'x1'	=> $max,
#		'y0'	=> $min,	'y1'	=> $max,
#	}],
	'logscale',			# specify which axes should be logscale, e.g. logscale => ['x', 'y']
	'plot_params',
	'text',				# an array of text items to add, each text addition should look like [x, y, text]
	'title',
	'xlabel',
	'ylabel',
	'xlim',		# array of min,max for each
	'xmin', 		# 1 number
	'xshift',	# 1 number, for annotate; how much % of the xaxis is shifted for default arrows
	'ylim',	# array of min,max for each
	);
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args;
		die "The above args are accepted by $current_sub"
	}
	my @keys;
	if (defined $args->{'keys'}) {
		@keys = @{ $args->{keys} };
	} else {
		@keys = sort {lc $a cmp lc $b} keys %{ $args->{data} };
	}
	if (scalar @keys != 3) {
		p $args->{data};
		die 'There must be exactly 3 data keys for a 2D scatterplot with color as a 3rd dimension';
	}
	if (defined $args->{color_key}) {
		while (my ($i, $key) = each @keys) {
			next unless $key eq $args->{color_key};
			splice @keys, $i, 1; # remove the color key from @keys
		}
	} else {
		$args->{color_key} = pop @keys;
	}
	if (not defined $args->{data}{$args->{color_key}}) {
		p $args->{data};
		die "\"$args->{color_key}\" isn't defined as a key in \"data\"";
	}
	my %lengths = map { $_ => scalar @{ $args->{data}{$_} }} (@keys, $args->{color_key});
	my %unique_lengths = map { $_ => 1 } values %lengths;
	if (scalar keys %unique_lengths != 1) {
		p %lengths;
		die 'There is some mismatch in entry data, all keys should have the same length';
	}
	my $unlink = 0;
	my ($fh, $tempfile) = tempfile (UNLINK => $unlink, SUFFIX => '.py', DIR => '/tmp');
	say $tempfile if $unlink == 0;
	say $fh "# originally written by $current_sub in " . __FILE__ if $unlink == 0;
	say $fh 'import matplotlib.pyplot as plt';
	say $fh 'x = [' . join (',', @{ $args->{data}{$keys[0]} }) . ']';
	say $fh 'y = [' . join (',', @{ $args->{data}{$keys[1]} }) . ']';
	say $fh 'z = [' . join (',', @{ $args->{data}{$args->{color_key}} }) . ']';
	$args->{xlabel} = $args->{xlabel} // $keys[0];
	$args->{ylabel} = $args->{ylabel} // $keys[1];
	say $fh "plt.xlabel('$args->{xlabel}')";
	say $fh "plt.ylabel('$args->{ylabel}')";
	say $fh "plt.title('$args->{title}')" if defined $args->{title};
	if (defined $args->{plot_params}) {
		my @params;
		while (my ($key, $value) = each %{ $args->{plot_params} }) {
			if (looks_like_number($value)) {
				push @params, "$key = $value";
			} else {
				push @params, "$key = '$value'";
			}
		}
		say $fh 'im = plt.scatter(x, y, c = z, cmap = "gist_rainbow", ' . join (',', @params) . ')';
	} else {
		say $fh 'im = plt.scatter(x, y, c = z, cmap = "gist_rainbow")';
	}
	say $fh "plt.colorbar(im, label = '$args->{color_key}')";
	if (defined $args->{logscale}) {
		if (ref $args->{logscale} ne 'ARRAY') {
			p $args->{logscale};
			die 'logscale should be an array of axes like ["x", "y"]' # this error message is more meaningful
		}
		foreach my $axis (@{ $args->{logscale} }) { # x, y 
			say $fh "plt.$axis" . 'scale("log")';
		}
	}
	if (defined $args->{line_segment}) {
		foreach my $ls (@{ $args->{line_segment} }) { # hash ref
			my @reqd_coords = ('x0', 'x1', 'y0', 'y1');
			my @undef_args = grep {!defined $ls->{$_}} @reqd_coords;
			if (scalar @undef_args > 0) {
				p @undef_args;
				say STDERR 'the above args are not defined.';
				say STDERR 'the below args are defined:';
				p $ls;
				p @reqd_coords;
				die 'the above coordinates must be defined.'
			}
			my $x = join (',', @{ $ls }{'x0', 'x1'});
			my $y = join (',', @{ $ls }{'y0', 'y1'});
			my $linestyle = $ls->{linestyle} // 'dashed';
			my $color = $ls->{color} // 'red';
			my $options = "linestyle = '$linestyle', color = '$color'";
			if (defined $ls->{label}) {
				$ls->{label} =~ s/'/\\'/g;
				$options .= ", label = '$ls->{label}'";
			}
			say $fh "plt.plot([$x], [$y], $options)";
		}
	}
	$args->{legend} = $args->{legend} // 0;
	if (defined $args->{custom_legend}) {
# https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
		say $fh 'legend_elements = [' . join (',', @{ $args->{custom_legend} }) . ']';
		say $fh 'plt.legend(handles = legend_elements)';
	} elsif ($args->{legend} > 0) { # 0 turns legend off
		say $fh 'plt.legend()';
	}
	say $fh "plt.clim(vmin = $args->{cb_min})" if defined $args->{cb_min};
	say $fh "plt.clim(vmax = $args->{cb_max})" if defined $args->{cb_max};
	foreach my $arg (grep { defined $args->{$_} } ('xlim', 'ylim')) {
		say $fh "plt.$arg(" . join (',', @{ $args->{$arg} }) . ')';
	}
	if (defined $args->{text}) {
		my %args = (
			object	=> 'plt',
			fh			=> $fh,
			text		=> $args->{text}
		);
		foreach my $arg (grep {defined $args->{$_}} ('adjust_text', 'annotate', 'xshift')) {
			$args{$arg} = $args->{$arg};
		}
		text({%args});
	}
	if (defined $args->{function}) {
		$args->{function}{fh} = $fh;
		function( $args->{function} );
	}
	if (defined $args->{xmin}) {
		say $fh "plt.xlim(left = $args->{xmin})";
	}
	say $fh "plt.savefig('$args->{filename}', bbox_inches='tight', pad_inches = 0.1, metadata={'Creator': 'made/written by Util.pm\\'s $current_sub " . getcwd() . '/' . __FILE__ . '\'})';
	close $fh;
	my ($stdout, $stderr, $exit) = capture {
		system( "python3 $tempfile" );
	};
	if ($exit != 0) {
		say "exit = $exit";
		say "STDOUT = $stdout";
		say "STDERR = $stderr";
		die "\"python3 $tempfile\" failed";
	}
	say 'wrote ', colored(['blue on_bright_yellow'], $args->{filename});
	return $args->{filename}
}
my $who = json_file_to_ref($args{target_json});
sub plot_lengths {#($data, $svg_name, $title, $text) {
	my $data			= shift;
	my $svg_name 	= shift;
	my $title		= shift;
	my $text			= shift;
=d{
    deg_key                 "DEG20010001",
    domains                 "B-block_TFIIIC (PFAM 35.0); disorder_prediction (MOBIDB_LITE 2.0) (evalue undef); disorder_prediction (MOBIDB_LITE 2.0) (evalue undef); Tau138_eWH (CDD 3.20) (evalue undef)",
    end                     553,
    gene_name               "TFC3",
    length                  551,
    long_name               "transcription factor TFIIIC subunit TFC3",
    max_similarity          0.0666666666666667,
    "mean # of pathogens"   1.79166666666667,
    mean_similarity         0.0132392914653784,
    min_similarity          0,
    start                   2
}
=cut
	my (@length, @mean_pathogens, @mean_similarity, @text, %deg);#, @full_protein_length);
	foreach my $row (@{ $data }) {
		my $key = $row->{deg_key};
#		next if $row->[3] < 10; # next if length < 10
		if ((defined $deg{$key}) && ($deg{$key}{'length'} < $row->{'length'})) {
			$deg{$key} = $row;
		} elsif (not defined $deg{$key}) {
			$deg{$key} = $row;
		}
	}
	my %gene_labels;
	if ((defined $text->[0]) && ($text->[0] eq 'всё')) {
		foreach my $r (@{ $data }) {
			$gene_labels{$r->{gene_name}} = 1;
		}
	} else {
		%gene_labels = map {$_ => 1} @{ $text };
	}
#	my $xint = 650;
	foreach my $key (sort keys %deg) {
		push @length, $deg{$key}{'length'};
		push @mean_pathogens, $deg{$key}{'mean # of pathogens'};
		push @mean_similarity, $deg{$key}{mean_similarity};
#		push @full_protein_length, $deg_length->{$key};
		unless (defined $gene_labels{$deg{$key}{gene_name}}) {
			next;
		}
#https://matplotlib.org/stable/gallery/text_labels_and_annotations/text_commands.html
		my $style = 'fontsize = 9';#weight="bold", style = "italic", 
#		if ($deg{$key}{gene_name} =~ m/^(?:YEF3|TRL1|FAS[12])$/) {
#			$style = 'style = "italic", fontsize = 9';
#		}
		$deg{$key}{gene_name} = ucfirst lc $deg{$key}{gene_name};
		push @text, [
								$deg{$key}{'length'},					# x
								$deg{$key}{'mean # of pathogens'},	# y
								$deg{$key}{gene_name},					# text
								$style
						];
			if ($deg{$key}{gene_name} eq 'FBA1') { # ad hoc
				$text[-1][1] -= 0.3; # add to y
			} elsif ($deg{$key}{gene_name} eq 'AUR1') {
				$text[-1][0] -= 90
			}
	}
	scatterplot_2d_color({
		cb_min		=> 0,
		cb_max		=> 1,
		filename		=> $svg_name,
		title			=> $title,
		data			=> {
			'Alignment length (a.a.)'			=> [@length],
#			'Total Protein Length (a.a.)'		=> [@full_protein_length],
			'Mean Pathogen Column Occupancy'	=> [@mean_pathogens],
			'Mean Pathogen ϕ'						=> [@mean_similarity]
		},
		color_key	=>	'Mean Pathogen ϕ',#'Total Protein Length (a.a.)',
		plot_params => {
			alpha => 0.5
		},
		text			=> [@text],
		xmin			=> 0,
#		xlim			=> [0, 2290],
#		ylim			=> [1,10.2]
	});
}

my $ortholog_db;
if (defined $args{exclude}) {
	$ortholog_db = json_file_to_ref($args{exclude});
}
my %ortholog_genes = map {$_ => 1} keys %{ $ortholog_db };
my @indices = reverse 0..$#$who;
foreach my $r (grep {defined $ortholog_genes{$who->[$_]{gene_name}}} @indices) {
	splice @{ $who }, $r, 1; # remove genes that have human orthologs
}
my @protein_labels;
if (defined $args{gene_labels}) {
	@protein_labels = split ',', $args{gene_labels};
}
say 'there are ' . scalar @{ $who } . ' genes without orthologs in WHO data';
plot_lengths(
	$who,
	$args{output_svg},
	$args{title},
	 \@protein_labels #['всё'],
);
