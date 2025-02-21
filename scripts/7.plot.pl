#!/usr/bin/env perl

use 5.038;
use warnings::unused;
use warnings FATAL => 'all';
use autodie ':default';
use List::Util qw(max min);
use Util qw(list_regex_files json_file_to_ref violin_plot scatterplot_2d_color scatterplot_2d_groups_color pvalue scatterplot);

sub first_letter ($string) {
	if ($string =~ m/([A-Za-z])/) {
		return uc $1;
	}
}

my %violin;
foreach my $j (list_regex_files('\.json$', '/home/con/identify.target/who.critical/blastp.results')) {
	next if $j =~ m/fungal.blastp.hits\.json$/;
	next if $j =~ m/sambucinum\.json$/;
	my $species;
	if ($j =~ m/
					([^\/\.]+)	# the genus, excluded characters right-slash and period \$1
					\.				# the separator between genus and species
					([^\/]+)		# the species name, which can have ".", saved as \$2
					\.json$		# file suffix
					/xx) {
		$species = ucfirst first_letter($1) . ".$2";
	} else {
		die "$j failed regex.";
	}
#	if ($j =~ m/\/([^\/]+)\.json$/) {
#		$species = $1;
#	} else {
#		die "$j failed regex.";
#	}
	say $j;
	my $blast = json_file_to_ref( $j );
	if (defined $blast->{BlastOutput2}) {
		$blast = $blast->{BlastOutput2}; # simplify hash
	}
	my %keys;
	foreach my $n ('0555', '1054', '0729', '0294', '0455') {
		$keys{"DEG2001$n"} = 1;
	}
	foreach my $report (@{ $blast }) {
		if (
				(defined $report->{report}{results}{search}{message})
				&&
				($report->{report}{results}{search}{message} eq 'No hits found')
			) {
#			say STDERR "$report->{report}{results}{search}{query_title} has no hits found for $j";
			next;
		}
		my $query = $report->{report}{results}{search}{query_title};
		next unless defined $keys{$query};# eq 'DEG20010555';
#		if (not defined $query) {
#			die "\$report->{report}{results}{search}{query_title} isn't defined.";
#		}
		$report = $report->{report}{results}{search}{hits}; # simplify hash
		foreach my $hit (@{ $report }) {
			foreach my $pair ( @{ $hit->{hsps} } ) {
				foreach my $metric (grep {defined $pair->{$_}} ('align_len', 'bit_score', 'gaps', 'evalue', 'identity', 'score')) {
					push @{ $violin{$query}{$metric}{$species} }, $pair->{$metric};
				}
			}
		}
	}
}
my $deg_info = json_file_to_ref( '../json/deg.annotation.json' );
my @species = map {ucfirst $_} reverse qw(g.max h.sapiens o.sativa solanum.tuberosum zea.mays a.fumigatus c.albicans c.auris c.parapsilosis c.tropicalis
c.neoformans.B.3501 n.glabratus c.neoformans.grubii.H99 c.neoformans.JEC21 h.capsulatum);
my @hosts = map {ucfirst $_} qw(g.max homo.sapiens o.sativa s.tuberosum z.mays);
my @pathogens = map {ucfirst $_} qw(a.fumigatus c.albicans c.auris c.parapsilosis c.tropicalis c.neoformans.B.3501 nakaseomyces.glabratus c.neoformans.grubii.H99 c.neoformans.JEC21 h.capsulatum);
my %p;
foreach my $query (sort keys %violin) {
	my (%scatterplot_data, %scatterplot_group, %extrema, %combined_scatterplot);
	my @axes = ('evalue', 'score', 'align_len');
	foreach my $axis (@axes) { # initialize the components
		$extrema{$axis}{min} = 'inf';
		$extrema{$axis}{max} = '-inf';
	}
	foreach my $metric (@axes) {
#	while (my ($iterator, $metric) = each @axes) {
		foreach my $host (@hosts) {
			my @val = grep {defined} @{ $violin{$query}{$metric}{$host} };
			push @{ $scatterplot_data{hosts}{$metric} }, @val;
			push @{ $scatterplot_group{hosts}{$host}{$metric} }, @val if scalar @val > 0;
			$extrema{$metric}{min} = min( $extrema{$metric}{min}, @val);
			$extrema{$metric}{max} = max( $extrema{$metric}{max}, @val);
			if ($metric =~ m/^(?:evalue|align_len)$/) {
				push @{ $combined_scatterplot{hosts}{$metric} }, @val if scalar @val > 0;
			}
		}
		foreach my $pathogen (@pathogens) {
			my @val = grep {defined} @{ $violin{$query}{$metric}{$pathogen} };
			push @{ $scatterplot_data{pathogens}{$metric} }, @val;
			push @{ $scatterplot_group{pathogens}{$pathogen}{$metric} }, @val if scalar @val > 0;
			$extrema{$metric}{min} = min( $extrema{$metric}{min}, @val);
			$extrema{$metric}{max} = max( $extrema{$metric}{max}, @val);
			if ($metric =~ m/^(?:evalue|align_len)$/) {
				push @{ $combined_scatterplot{pathogens}{$metric} }, @val if scalar @val > 0;
			}
		}
	}
	foreach my $group (keys %scatterplot_group) {
		foreach my $species (keys %{ $scatterplot_group{$group} }) {
			foreach my $val (@{ $scatterplot_group{$group}{$species}{evalue} }) {
				$val += 10**(-180);
			}
		}
	}
	my $gene_name = $deg_info->{$query}{'Gene Name'};
	scatterplot({
		filename	=> 'simple_scatterplot',
		data		=> \%combined_scatterplot,
		logscale	=> ['y'],
		title		=> "WHO Fungal Pathogens Alignment with $query ($gene_name)",
		xlabel	=> 'Alignment Length (amino acids)',
		ylabel	=> 'log(Expection Value+10⁻¹⁸⁰)',
	});
	die;
	say $extrema{evalue}{max} + 0.2*($extrema{evalue}{max} - $extrema{evalue}{min});
	foreach my $group ('hosts', 'pathogens') {
		my $leg_pos = 'outside center left';
		$leg_pos = 'outside center right' if $group eq 'pathogens';
		scatterplot_2d_groups_color({
			data				=> $scatterplot_group{$group},
			color_key		=> 'score',
			'color.min'		=> $extrema{score}{min},
			'color.max'		=> $extrema{score}{max},
			'x.axis.key'	=> 'align_len',
			'y.axis.key'	=> 'evalue',
			filename			=> join ('_', $query, 'group_scatterplot', $group),
#			title				=> "($gene_name): $group",
			xlabel			=> 'Alignment length (a.a.)',
			ylabel			=> 'log(Expection Value+10⁻¹⁸⁰)',
			'legend.position' => $leg_pos,
			logscale			=> ['y'],
			xmin				=> 0,
			xmax				=> $extrema{align_len}{max} + .05*($extrema{align_len}{max}-$extrema{align_len}{min}),
			ymin				=> 10**(-182),
			ymax				=> 10**7,
		});
#		scatterplot_2d_color({
#			color_key=> 'score',
#			data		=> $scatterplot_data{$group},
#			filename	=> join ('_', $query, 'scatterplot', $group),
#			text		=> $plot_text{$group},
#			title		=> "$query ($gene_name): $group"
#		});
	}
	die;
	foreach my $metric (sort keys %{ $violin{$query} }) {
		my $x_label = $metric;
		my (@host_vals, @pathogen_vals);
		foreach my $host (grep {defined $violin{$query}{$metric}{$_}} @hosts) {
			if (scalar @{ $violin{$query}{$metric}{$host} } == 0) {
				delete $violin{$query}{$metric}{$host};
				next;
			}
			push @host_vals, min(@{ $violin{$query}{$metric}{$host} });
		}
		foreach my $path (grep {defined $violin{$query}{$metric}{$_}} @pathogens) {
			if (scalar @{ $violin{$query}{$metric}{$path} } == 0) {
				delete $violin{$query}{$metric}{$path};
				next;
			}
			push @pathogen_vals, min(@{ $violin{$query}{$metric}{$path} });
		}
		say "\$p{$query}{$metric}";
		$p{$query}{$metric} = pvalue(\@host_vals, \@pathogen_vals);
		my @defined_species = grep {defined $violin{$query}{$metric}{$_}} @species;
		violin_plot({
			filename => $metric . "_violin_plot_$query",
			flip		=> 1,
			title		=> "$query ($gene_name): Alignment Distribution of $metric",
			data		=> $violin{$query}{$metric},
			'keys'	=> \@defined_species,
			ylabel	=> $x_label,
			xlabel	=> 'Species'
		});
	}
}
p %p;
