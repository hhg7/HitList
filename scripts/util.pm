package util;
use FindBin 1.51 qw( $RealBin );
use DDP;
use Capture::Tiny 'capture';
use List::Util 'min';
use JSON qw(encode_json decode_json);
our @EXPORT = qw(comments execute json_file_to_ref ref_to_json_file msa_phylo_plot break_long_string_into_array fasta2hash list_regex_files plots_in_rows seq_type write_2d_array_to_tex_tabular);
use feature 'say';
use Term::ANSIColor;
use File::Temp 'tempfile';
use Cwd 'getcwd';
use DDP;
use Exporter 'import';

sub widecell2tabular {#($text, $width = 'inf', $align = 'c') { # helper function; not exported
	my $text = shift;
	my $width = 'inf';
	my $align = 'c';
	return $text if $width == 'inf';
	return $text if (length($text) < $width);
	my @text = split '', $text;
	my $return = '\begin{tabular}{' . "$align}";
	for (my $l = 0; $l <  (scalar @text) / $width; $l++) {
		my $i1 = $width * $l;
		my $i2 = min($i1 + $width-1, scalar @text - 1);
		my $line = join ('', @text[$i1..$i2]);
		if ( # should a ending "-" be added?
				(defined $text[$i2+1])  &&
				($text[$i2]   =~ m/\S/) &&
				($text[$i2+1] =~ m/\S/)
			) {
			$line .= '-';
		}
		$return .= "$line\\\\\n";
	}
	$return .= '\end{tabular}';
	return $return;
}

sub fasta2hash {#($file, $key = undef) {
=purpose
	this takes a fasta file as input
	if $key is not defined, returns a hash reference of the fasta sequence of all keys
	if $key is defined, returns only a string the fasta sequence
=cut
	my $file = shift;
	my $key = undef;
	die "\"$file\" doesn't exist or isn't a readable file" unless -f -r $file;
	my ($description, %r);
	if ($file =~ m/\.gz$/) {
		open my $zcat, "zcat $file|";
		while (<$zcat>) {# read through the file to get the fasta sequence(s)
      	chomp;
		   if (/^>(.+)/) {
				$description = $1;
	     	} else {
				die "\$description isn't defined" unless defined $description;
				if (defined $key) { # I'm only saving a single key to save RAM
					$r{$description} .= $_ if $description eq $key;
				} else {
					$r{$description} .= $_;
				}
			}
			if ((defined $key) && (defined $r{$description}) && ($description ne $key)) {
				last; # stop reading the file once the selected sequence has been extracted
			}
		}
	} else {
		open my $fh, '<', $file;
      while (<$fh>) { # read through the file to get the fasta sequence(s)
         chomp;
         if (/^>(.+)/) {
            $description = $1;
         } else {
           die "\$description isn't defined" unless defined $description;
           if (defined $key) { # I'm only saving a single key to save RAM
					$r{$description} .= $_ if $description eq $key;
				} else {
					$r{$description} .= $_;
				}
	      }
			if ((defined $key) && (defined $r{$description}) && ($description ne $key)) {
				last; # stop reading the file once the selected sequence has been extracted
			}
	   }
	}
	if (scalar keys %r == 0) {
		my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
		die "couldn't find any keys for $file with $current_sub";
	}
	if (defined $key) {
		return $r{$key} if defined $r{$key};
		my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
		die "Couldn't find $key in $file using $current_sub";
	}
	return \%r
}
sub break_long_string_into_array {#($string, $width = 80) {
	my $string = shift;
	my $width = 80;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	die "string is not defined in $current_sub" unless defined $string;
	my @return;
	my $len = length($string)-1;
	my @seq = split '', $string;
	for (my $l = 0; $l < length($string)/$width; $l++) {
		my $i1 = $width*$l;
		my $i2 = min($i1+$width-1, $len);
		my $substring = join ('', @seq[$i1..$i2]);
		push @return, $substring;
	}
	return \@return
}

sub execute {#($cmd, $return = 'exit', $die = 1) {
	my $cmd = shift;
	my $return = 'exit';
	my $die = 1;
	if ($return !~ m/^(exit|stdout|stderr|all)$/i) {
		die "you gave \$return = \"$return\", while this subroutine only accepts ^(exit|stdout|stderr)\$";
	}
	my ($stdout, $stderr, $exit) = capture {
		system( $cmd )
	};
	if (($die == 1) && ($exit != 0)) {
		say STDERR "exit = $exit";
		say STDERR "STDOUT = $stdout";
		say STDERR "STDERR = $stderr";
		die "$cmd\n failed";
	}
	if ($return =~ m/^exit$/i) {
		return $exit
	} elsif ($return =~ m/^stderr$/i) {
		chomp $stderr;
		return $stderr
	} elsif ($return =~ m/^stdout$/i) {
		chomp $stdout;
		return $stdout
	} elsif ($return =~ m/^all$/i) {
		chomp $stdout;
		chomp $stderr;
		return {
			exit   => $exit, 
			stdout => $stdout, 
			stderr => $stderr
		}
	} else {
		die "$return broke pigeonholes"
	}
	return $stdout
}

sub json_file_to_ref {#($json_filename) {
	my $json_filename = shift;
	die "$json_filename doesn't exist or isn't a file" unless -f $json_filename;
	die "$json_filename has 0 size" if (-s $json_filename == 0);
#	say "Reading $json_filename" if defined $_[0];
	open my $fh, '<:raw', $json_filename; # Read it unmangled
	local $/;                     # Read whole file
	my $json = <$fh>;             # This is UTF-8
#	$json =~ s/NaN/"NaN"/g;
	return decode_json($json); # This produces decoded text
}
sub ref_to_json_file {#($ref, $json_filename) {
	my $ref = shift;
	my $json_filename = shift;
	my $ref_json_filename = ref $json_filename;
	unless ($ref_json_filename eq '') {
		die "$json_filename isn't a scalar/string";
	}
	open my $fh, '>:raw', $json_filename; # Write it unmangled
	say $fh encode_json($ref);
	say 'Wrote ' . colored(['blue on_red'], $json_filename);
	return $json_filename;
}

sub hash2fasta_file {
	my $hash = shift;
	my $filename = shift;
	my $order = shift;
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	unless (ref $hash eq 'HASH') {
		die "1st argument to $current_sub must be a hash";
	}
	unless (ref $filename eq '') {
		die "1st argument to $current_sub must be a scalar/string";
	}
	my @order;
	if (defined $order) {
		@order = @{ $order };
	} else {
		@order = sort keys %{ $hash };
	}
	open my $fh, '>', $filename;
	foreach my $key (@order) {
		say $fh ">$key";
		my $aa_or_nt_array = break_long_string_into_array($hash->{$key});
		say $fh join ("\n", @{ $aa_or_nt_array });
	}
	return $filename;
}

sub file2string {
	my $file = shift;
	open my $fh, '<', $file;
	return do { local $/; <$fh> };
}

sub msa_phylo_plot {
=purpose
	this script takes a sequence fasta, and a blast database, makes an alignment from them to visualize with my own python3 script and R
	the blast database is only to find the protein, this subroutine, and its helper above, find the protein sequence above, and make a multiple sequence alignment input
	
	-inputs:
	
	'active.site.aa' (optional):
	
		plot a vertical line going from top to bottom for these amino acids.  The argument is given as an anonymous hash:
		'active.site.aa'  => {
			His395 => 395,
			Asn924 => 924,
			Asn702 => 702
		}
	where the label is the key, and the value is the MSA-position
	
	blast.json (required):
		a hash summary of alignment results, organized by {species}{query}. Could be either the filename itself, or the hash data structure
	
	db (required):
		a hash; keys are species, the next key are protein names, and their values are FASTA sequences
	
	debug (optional):
		don't make output files.  The clustalo command is very time-consuming.  Don't waste computer time if I'm just debugging another part.  No external commands will be called, will only do the quick parts
	
	msa.image.filename (optional):
		the multiple sequence alignment output image in SVG format, which is made by my altered CIAlign script. The file will be created if the output filename is specified
	
	msa.xlabel (optional):
		xlabel for Multiple Sequence Alignment plot. The default xlabel is "Amino Acid Residue"
	
	msa.ylabel (optional):
		ylabel for Multiple Sequence Alignment plot. The default xlabel is "Species"
	
	out.fasta (optional):
		the multiple sequence alignment fasta, which will be output by clustalo.  Must end in ".input.fa"
	
	phylo.image.filename (optional):
		the output filename of the phylogenetic tree image
	
	query (required)
		the query that was used in the blast search
	
	query.fasta (required):
		the sequence of the query protein/nucleic acid
	
	species.only (optional):
		only use the species name in the plots, in case the whole names don't fit
	
	species.order (optional):
		present the species in this order for clarity (e.g. fungi apart from plants)
	
	tex.outfile (optional)
		output the data into a LaTeX table
	
	title (optional):
		the title for both plots
=cut
	my ($args) = @_;
	unless (ref $args eq 'HASH') {
		my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
		die "args must be given as a hash ref, e.g. \"$current_sub({ query => 'DEG2023423' })\"";
	}
	unless (ref $args->{db} eq 'HASH') {
		p $args->{db};
		die "the \"db\" argument must be a hash, but a " . ref $args->{db} . ' was entered.';
	}
	my @reqd_args = ('blast.json', 'db', 'out.fasta', 'query', 'query.fasta');
	my @undef_args = grep {!defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined';
	}
	my @defined_args = (@reqd_args, 'active.site.aa',
	'debug', 'msa.image.filename', 'msa.xlabel', 'msa.ylabel', 'out.fasta',
	'run.clustal', # do not run clustal if this is 0
	'phylo.image.filename', 'phylo.image.width',
	'query.species',	# substitute the query name for a species
	'species.only', 'species.order',
	'tex.outfile',
	'title');
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args;
		die 'The above args are accepted.'
	}
	if (
			(defined $args->{'out.fasta'})
			&&
			($args->{'out.fasta'} !~ m/\.input\.fa$/)
		) {
		 die "\"out.fasta\" must end in \".input.fa\"";
	}
=if ((defined $args->{'query.fasta'}) && (-f $args->{'query.fasta'})) { # it's required, I want to localize the data
		my $test_fasta = fasta2hash( $args->{'query.fasta'} );
		my $n_entries = scalar keys %{ $test_fasta };
		if ($n_entries != 1) {
			p $test_fasta;
			die "The test fasta must only have 1 entry, but there is an errant file above with $n_entries entries";
		}
	}
=cut
	my $home_dir = execute('echo $HOME', 'stdout');
	my $clustalo = "$RealBin/clustalo-1.2.4-Ubuntu-x86_64";
	my $py3 = "$RealBin/x.my.align.py";
	if ((defined $args->{debug}) && ($args->{debug} > 0)) {
		die "$clustalo isn't accessible, doesn't exist, or isn't executable" unless -f -x $clustalo;
		die "$py3 isn't readable, or doesn't exist." unless -f -r $py3;
	}
# read blastp alignment data to hash
	my ($all, $out, $filename, @order, @data_table, %reformat);
# format the names into scientific format
	foreach my $i (0..scalar @{ $args->{'species.order'} }-1) {
		my $newformat;
		if ((defined $args->{'query.species'}) && ($args->{'species.order'}[$i] eq $args->{query})) {
			$newformat = $args->{'query.species'};
		} else {
			$newformat = $args->{'species.order'}[$i];
		}
		$newformat =~ s/\./ /g; # remove dots for subspecies/strain
		if ($newformat =~ m/^([A-Z]+)\s+([a-z]+)(\s*)(.*)/) { # re-do the dot
			$newformat = "$1. $2$3$4";
		}
		$reformat{$args->{'species.order'}[$i]} = $newformat;
	}
# only read in the JSON if it's a string
	if ((ref $args->{'blast.json'}) eq '') { # $args->{'blast.json'} is a string
		$all = json_file_to_ref($args->{'blast.json'});
	} elsif ((ref $args->{'blast.json'}) eq 'HASH') { # args->{'blast.json'} is the data itself, don't use the disk again to re-read it repeatedly
		$all = $args->{'blast.json'};
	} else {
		my $ref = ref $args->{'blast.json'};
		die "'blast.json' should be either a hash or a string, but you entered a $ref";
	}
	if (defined $args->{'out.fasta'}) {
		$filename = $args->{'out.fasta'};
	} else {
		($out, $filename) = tempfile( DIR => '/tmp', SUFFIX => '.fa', UNLINK => 1);
	}
	open $out, '>', $args->{'out.fasta'};
	if (
#			(length $args->{'query.fasta'} <= 255) &&
			(seq_type($args->{'query.fasta'}) eq 'other')
		) { # it's a file.  This will break if the file is only amino acid/nucleotide characters
		print $out file2string($args->{'query.fasta'});
	} else {
#		say $out ">$args->{query}"; # title line
		if (defined $args->{'query.species'}) {
			say $out ">$args->{'query.species'}";
		}
#		say $out '>S.cerevisiae';
		say $out $args->{'query.fasta'}; # sequence itself
	}
	if (defined $args->{'species.order'}) {
		@order = @{ $args->{'species.order'} };
	} else {
		@order = sort { uc $a cmp uc $b } keys %{ $all };
	}
	p @order;
	$args->{'species.only'} = $args->{'species.only'} // 1; # by default, only the species is included
	my @stat = ('evalue', 'align_len', 'bit_score', 'identity', 'positive', 'score', 'gaps', '% identity', '% positive');
	my $total_sequences = 0; # I need to keep track of the # of sequences, clustal will die if only 1 sequence is shown
	foreach my $species (grep {$_ ne $args->{query}} @order) { # write species to the input fasta file, being written to $out; query is already written
		if (not defined $args->{db}{$species}) {
			die "\$args->{db}{$species} isn't defined @ " . __FILE__ . ' line ' . __LINE__;
			if (defined $args->{'tex.outfile'}) {
				push @data_table, [$species, ('-') x (scalar @stat + 2)];
			}
			next
		}
		my @proteins = grep {$_ eq $args->{query}} keys %{ $all->{$species} };
		if (scalar @proteins == 0) {
			say "$args->{query} not found for $species @ " . __FILE__ . ' line ' . __LINE__;
			push @data_table, [$species, ('-') x (scalar @stat + 2)];
			next
		}
		# for each query, get the fasta sequence for the protein that hit
		foreach my $protein (@proteins) {
			my $title;
			foreach my $d (@{ $all->{$species}{$protein}{description} }) {
				$title = $d->{title};
				last
			}
			my $seq_type = 'nt';
			$seq_type = 'aa' if $args->{db}{$species} =~ m/protein\.fa/;
			my $seq;
			if ($seq_type eq 'aa') { # get the protein sequence from the proteome file
#				say "Getting $title from \$args->{db}{$species}";
#				$seq = get_protein_sequence($args->{db}{$species}, $title);
				p $args->{db}{$species};
				foreach my $title (grep {$_ eq $title} keys %{ $args->{db}{$species} }) {
					$seq = break_long_string_into_array( $args->{db}{$species}{$title} );
				}
			} elsif ($seq_type eq 'nt') {
				@{ $seq } = $all->{$species}{$protein}{hseq};
				unshift @{ $seq }, ">$species $protein";
			}
			die "no sequence found for\n\n$title\n\nwithin\n\n$args->{db}{$species}" if scalar @{ $seq } == 0;
			if (grep {/\*/} @{ $seq }) {
				say STDERR "$species->$protein has a stop codon";
				next;
			}
			if ($args->{'species.only'} == 1) {
				$seq->[0] = ">$species";
			} else {
				$seq->[0] = ">$species $protein";
			}
			$seq->[0] =~ s/\s+/_/g;
			$seq->[0] =~ s/\'/p/g;
			foreach my $line ( @{ $seq } ) {
				say $out $line;
			}
			$total_sequences++;
			shift @{ $seq }; # discard title line
			$seq = join ('', @{ $seq });
			if (defined $args->{'tex.outfile'}) {
				#my $fasta;# = fasta2hash( $args->{'query.fasta'} );
				my $query_length = length $args->{'query.fasta'};
				#die "\"$args->{query}\" has 0 length from $args->{'query.fasta'}" if $query_length == 0;
				my $protein_name = $title;
				$protein_name =~ s/^\S+\s+//;
				$protein_name =~ s/[\[\]]+//g;
				$protein_name =~ s/_/ /g;
				my $space_species = $species;
				$space_species =~ s/\./ /g;
				$protein_name =~ s/$space_species\s?.*//ig;
				$all->{$species}{$protein}{evalue} = sprintf('%.2g', $all->{$species}{$protein}{evalue});
				push @data_table, [
					$species,
					$protein_name,
					length $seq,
					@{ $all->{$species}{$protein} }{'evalue', 'align_len', 'bit_score', 'identity', 'positive', 'score', 'gaps'}, # hash slice for simplicity, get the alignment data
					sprintf( '%.1f', 100 * $all->{$species}{$protein}{identity} / $query_length), # percent identity
					sprintf( '%.1f', 100 * $all->{$species}{$protein}{positive} / $query_length) # percent 
				];
			}
		}
	}
	close $out;
	say '@ ' . __FILE__ . '& line' . __LINE__;
	say 'wrote ' . colored(['black on_red'], $args->{'out.fasta'});
	if ($total_sequences < 2) {
		say colored(['yellow on_magenta'], "$args->{query} cannot run clustal, there's only 1 sequence");
		return {}
	}
	my %r;
	if (defined $args->{'tex.outfile'}) {
		@data_table = reverse @data_table; # the MSA image looks like this, the table should match
		write_2d_array_to_tex_tabular({
			data           => \@data_table,
			tex_filename   => $args->{'tex.outfile'},
			header         => ['Species', 'Hit Protein', 'Hit Length (a.a.)', @stat],
			'max.width'    => 30
		});
		$r{'tex.outfile'} = $args->{'tex.outfile'};
	}
	$args->{debug} = $args->{debug} // 0; # default for "debug" is off = 0
	if ($args->{debug} == 0) { # i.e. debug is off
#		say -s $args->{'out.fasta'};
		my $aln  = $args->{'out.fasta'};
		$aln =~ s/input\.fa$/output.aln/;
		my $tree = $args->{'out.fasta'};
		$tree =~ s/input\.fa$/output.newick/;
		my $log = $args->{'out.fasta'};
		$log =~ s/input\.fa$/clustal.log/;
#		say -s $args->{'out.fasta'};
		if ($aln eq $args->{'out.fasta'}) {
			die "\"$args->{'out.fasta'}\" cannot be identical to \"$aln\"; the files will overwrite";
		}
# include proteins that are absent from the clustal alignment as "-"
		$args->{debug} = $args->{debug} // 1;
		$args->{'run.clustal'} = $args->{'run.clustal'} // 1;
		say __FILE__ . ' @ LINE ' . __LINE__;
		if ($args->{'run.clustal'} > 0) {
			my $cmd = "$clustalo --threads 6 --in '$args->{'out.fasta'}' --outfmt=a2m --out '$aln' --guidetree-out='$tree' --force --log='$log' --dealign";
			say colored(['white on_magenta'], $cmd);
			say 'running clustalo...';
			execute($cmd);
			say 'finished clustalo...';
		}
		my $svg;
		if (defined $args->{'phylo.image.filename'}) {
			$svg = $args->{'phylo.image.filename'};
			my $ggtitle = '';
			$ggtitle = " + ggtitle(\"$args->{title}\")" if defined $args->{title};
			my $local_tree = $tree; # these are file names
			if (defined $args->{'query.species'}) {
				my $tree_data = file2string($tree);
				$tree_data    =~ s/$args->{query}/$args->{'query.species'}/g;
				foreach my $protein (keys %reformat) {
					$tree_data =~ s/$protein/$reformat{$protein}/g;
				}
				my ($tmp_tree, $tmpfilename) = tempfile(DIR => '/tmp', SUFFIX => '.newick', UNLINK => 0);
				say $tmpfilename;
				print $tmp_tree $tree_data;
				close $tmp_tree;
				$local_tree = $tmpfilename;
			}
			my $unlink = 0;
			my ($r, $r_filename) = tempfile(UNLINK => $unlink, DIR => '/tmp', SUFFIX => '.R');
			say $r 'library(ggtree)';
			say $r 'library(ggplot2)';
			say $r "tree <- read.tree('$local_tree')";
			$args->{'phylo.image.width'} = $args->{'phylo.image.width'} // 10;
			say $r "svg('$svg', width = $args->{'phylo.image.width'})";
			say $r 'options(ignore.negative.edge=TRUE)';
			say $r "ggtree(tree) + theme_tree2() + geom_tiplab(align=TRUE, linesize=0.5) +xlim(0,0.6) $ggtitle";
			say $r 'dev.off()';
			close $r;
			say colored(['magenta on_bright_yellow'], "R temp file: $r_filename") if $unlink == 0;
			say "running R to generate $svg";
			execute("Rscript $r_filename");
			say colored(['bright_green on_black'], "made $svg");
			$r{'phylo.image.filename'} = $svg;
		} else {
			$svg = $args->{'out.fasta'};
			$svg =~ s/fa$/phylo.svg/;
		}
		if (defined $args->{'msa.image.filename'}) {
			my $segments = '';
			if (defined $args->{'active.site.aa'}) {
				my $msa_query = fasta2hash( $aln, $args->{query} );
		# translate original sequence to its coordinate in multiple sequence alignment
				my %seq2msa; 
				my ($aa_count, $i) = (0, 0);
				foreach my $aa (split '', $msa_query) {
					if ($aa ne '-') {
						$seq2msa{$aa_count} = $i;
						$aa_count++;
					}
					$i++;
				}
				my %sites;
				while ( my ($aa,$nuc) = each %{ $args->{'active.site.aa'} }) {
					$sites{$aa} = $seq2msa{$nuc};
				}
				$segments = encode_json(\%sites);
			}
			my $title = '';
			$title = "--t \"$args->{title}\"" if defined $args->{title};
			$args->{'msa.xlabel'} = $args->{'msa.xlabel'} // 'Amino Acid Residue';
			$args->{'msa.ylabel'} = $args->{'msa.ylabel'} // 'Species';
			if ($segments ne '') {
				$segments = " --segments '$segments'";
			}
			my (@tmp_order, $msa_length);
			my $aln_hash = fasta2hash( $aln );
			foreach my $protein (keys %{ $aln_hash }, $args->{query}) {
				$msa_length = length $aln_hash->{$protein};
				last;
			}
#			p $aln_hash;
#			p @order, array_max => scalar @order;
#			p %reformat;
			foreach my $protein (grep {not defined $aln_hash->{$_} && $_ ne $args->{query}} @order) {
#				say "$protein isn't defined @ " . __FILE__ . ' line ' . __LINE__;
				$aln_hash->{$protein} = '-' x $msa_length;
			}
			my ($tmp, $fa_tempfile) = tempfile(DIR => '/tmp', UNLINK => 0, SUFFIX => '.fa');
			close $tmp;
=c			p $args->{'species.order'};
 [0] "DEG20010294",
 [1] "A.fumigatus",
 [2] "C.albicans",
 [3] "C.auris",
 [4] "C.tropicalis",
 [5] "C.parapsilosis",
 [6] "C.neoformans.B.3501A",
 [7] "C.neoformans.grubii.H99",
 [8] "C.neoformans.JEC21",
 [9] "H.capsulatum",
     (...skipping 2 items...)
=cut
			foreach my $key (@{ $args->{'species.order'} }) { # reformat
				my $new_key = $reformat{$key};
				$new_key =~ s/^(\w+)\.\s(\w+)/\$\\it\{$1. $2}\$/; # italicize genus & species
				say "$key\t$new_key";
				$aln_hash->{$new_key} = delete $aln_hash->{$key}; # rename, formatted
				push @tmp_order, $new_key;
			}
			$aln_hash->{$tmp_order[0]} = $aln_hash->{$args->{'query.species'}};
#			p @tmp_order, array_max => scalar @tmp_order;
#			p $aln_hash;
#			die;
			hash2fasta_file($aln_hash, $fa_tempfile, \@tmp_order);
			my $cmd = "python3 $py3 --f $fa_tempfile $title --o '$args->{'msa.image.filename'}' --x \"$args->{'msa.xlabel'}\" --y \"$args->{'msa.ylabel'}\" $segments";
			say colored(['white on_magenta'], $cmd);
			execute($cmd);
			say 'made ' . colored(['bright_blue on_black'], $args->{'msa.image.filename'});
			$r{'msa.image.filename'} = $args->{'msa.image.filename'};
		}
	}
	return \%r# https://4va.github.io/biodatasci/r-ggtree.html
}
=example execution
	msa_phylo_plot({
		'blast.json'           => $combined_blast,
		db                     => \%db,
		'msa.image.filename'   => $msa_image,
		'out.fasta'            => "$cwd/fa/clustal.$query.input.fa",
		'phylo.image.filename' => $phyl_image,
		query                  => $query,
		'query.fasta'          => "$home/identify.target/DEG/$row->[0].fa",
		'species.order'        => [@fungi, @hosts],
		'species.only'         => 1,
		title                  => $title
	});
=cut

sub comments {
	my $comment = 'This file was written by ' . getcwd() . '/' . __FILE__ . '. ';
	if (scalar @ARGV > 0) {
		$comment .= 'Command line arguments: ' . join (' ', @ARGV)
	} elsif (scalar @ARGV == 0) {
		$comment .= 'There were no command line arguments.'
	}
	$comment
}

sub list_regex_files {
	my $regex = shift;
	my $directory = '.';
	if (defined $_[0]) {
		$directory = shift;
	}
	my $case_sensitive = 'yes';
	my @files;
	opendir (my $dh, $directory);
	if ($case_sensitive eq 'yes') {
		$regex = qr/$regex/;
	} else {
		$regex = qr/$regex/i;
	}
	while (my $file = readdir $dh) {
		next if $file !~ $regex;
		next if $file =~ m/^\.{1,2}$/;
		my $f = "$directory/$file";
		if (-f $f) {
			if ($directory eq '.') {
				push @files, $file
			} else {
				push @files, $f
			}
		}
	}
	@files
}

sub plots_in_rows {
=example
each one of the sets "all", "hosts", and "pathogens" will become a row in the final figure
x axis is indicated by [0]
y axis is indicated by [1]
{
 data       {
     all         [
         [0] [
                 [0] 0,
                 [1] 1,
                 [2] 2,
                     (...skipping 1256 items...)
             ],
         [1] [
                 [0] 0,
                 [1] 0,
                 [2] 0,
                     (...skipping 1256 items...)
             ]
     ],
     hosts       [
         [0] [
                 [0] 0,
                 [1] 1,
                 [2] 2,
                     (...skipping 1256 items...)
             ],
         [1] [
                 [0] 0,
                 [1] 0,
                 [2] 0,
                     (...skipping 1256 items...)
             ]
     ],
 },
 filename   "svg/sneath.0729.colobar.heatmap.svg",
 title      "YEF3",
 xlabel     "MSA amino acid index",
}

or another example:
	plots_in_rows({
		filename	=> $similarity_svg,
		figwidth	=> 12,
		data		=> {
			hosts			=> [
				[0..$max_n],
				$group_similarity->{hosts}
			],
			pathogens	=> [
				[0..$max_n],
				$group_similarity->{pathogens}
			],
			domains		=> $domains
		},
		xlabel		=> 'MSA amino acid index',
		yaxis_off	=> ['domains'],
		ylabel		=> { # no ylabel for domains
			hosts			=> 'Φ',
			pathogens	=> 'Φ',
		},
		order			=> ['hosts', 'pathogens', 'domains'],
		xlim			=> [0,$max_n],
		ylim			=> [-0.05,1.05],
		title			=> $gene
	});
=cut
	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
	my ($args) = @_;
	my @reqd_args = (
		'data',		# hash of arrays of arrays
		'filename'	# include file type extension, e.g. ".svg"
	);
	my @undef_args = grep { !defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined.';
	}
	my @undef_grps = grep {not defined $args->{data}{$_}} (sort keys %{ $args->{data} });
	if (scalar @undef_grps > 0) {
		p $args->{data};
		p @undef_grps;
		die "the above groups are not defined in \"$current_sub\"";
	}
	my @defined_args = (
		'figheight', 'figwidth',
		'line.segments',	# call the helper function
		'order',
		'title',	# this is sent to suptitle; "title" won't work here
		'xlabel',# overall xlabel for all plots
		'xlim', 	# e.g. [0,1]
		'yaxis_off', # deactivate the y-axis for these groups, an array
		'ylabel',# a hash for each data group
		'ylim', 	# e.g. [0,1]
	@reqd_args);
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		print "the above arguments are not recognized.\n";
		p @defined_args;
		die 'The above args are accepted.'
	}
# https://stackoverflow.com/questions/45841786/creating-a-1d-heat-map-from-a-line-graph
	if (defined $args->{ylabel}) {
		my $ref = ref $args->{ylabel};
		unless ($ref eq 'HASH') {
			die "ylabel arg is a \"$ref\", while it must be a hash for each group in \"data\" for $current_sub";
		}
	}
	my @order;
	if (defined $args->{order}) {
		@order = @{ $args->{order} };
	} else {
		@order = sort {lc $a cmp lc $b} keys %{ $args->{data} };
	}
	$args->{figheight} = $args->{figheight} // 4.8;
	$args->{figwidth}  = $args->{figwidth}  // 6.4;
	my $unlink = 0;
	my ($py, $temp_py) = tempfile(DIR => '/tmp', SUFFIX => '.py', UNLINK => $unlink);
	say "temp file is $temp_py from " . __FILE__ . " in $current_sub" if $unlink == 0;
	say $py 'import matplotlib.pyplot as plt';
	my $nrows = scalar @order;
	my @ax = map {"ax$_"} 0..$#order;
	say $py 'fig, (' . join (',', @ax) . ") = plt.subplots(nrows = $nrows, sharex = True, figsize = ($args->{figwidth}, $args->{figheight}))";
	while (my ($n, $group) = each @order) {
#		say $py "ax$n.set_ylabel('$args->{ylabel}{$group}') " if defined $args->{ylabel}{$group};
		say $py "ax$n.set_title('$group')";
		if ((defined $args->{yaxis_off}) && (grep {$_ eq $group} @{ $args->{yaxis_off} })) {
			say $py "ax$n.get_yaxis().set_visible(False)";
		}
		if (ref $args->{data}{$group} eq 'ARRAY') { # simple 1 group plot
 # plotting conventionally
			say $py "ax$n.plot([" . join (',', @{ $args->{data}{$group}[0] }) .' ], [' . join (',', @{ $args->{data}{$group}[1] }) . '])';
			foreach my $arg (grep { defined $args->{$_} } ('xlim', 'ylim')) {
				say $py "ax$n.set_$arg(" . join (',', @{ $args->{$arg} }) . ')';
			}
		} elsif (
				(ref $args->{data}{$group} eq 'HASH') # line_segments or text
			) {
# plotting with helper functions
			if (defined $args->{data}{$group}{line_segments}) {
				my ($ymin, $ymax) = ('inf', '-inf');
				foreach my $segm (@{ $args->{data}{$group}{line_segments} }) {
					foreach my $pt (0,1) {
						$ymin = min($ymin, $segm->[$pt][1]);
						$ymax = max($ymax, $segm->[$pt][1]);
					}
				}
				line_segments({
					fh						=> $py,
					object				=> "axes$n",
					'line.segments'	=> $args->{data}{$group}{line_segments},
					ylim					=> [$ymin-1, $ymax+1]
				});
			}
			if (defined $args->{data}{$group}{text}) {
				text({
					fh			=> $py,
					text		=> $args->{data}{$group}{text},
					object	=> "axes$n",
				});
			}
		}
	}
	if (defined $args->{xlim}) {
		say $py "plt.xlim($args->{xlim}[0], $args->{xlim}[1])";
	}
	$args->{xlabel} = $args->{xlabel} // 'Amino Acid Residue (MSA Coordinate)';
	if (defined $args->{title}) {
		$args->{suptitle} = $args->{title};
	}
	foreach my $arg (grep {defined $args->{$_}} ('xlabel', 'suptitle')) {
		say $py "plt.$arg('$args->{$arg}')";
	}
	say $py 'fig.tight_layout(pad = 1.0)';
	say $py "plt.savefig('$args->{filename}', metadata={'Creator': 'made/written by " . getcwd() . "/$RealScript called using \"$current_sub\" in " . __FILE__ . "'})";
	close $py;
	execute("python3 $temp_py");
	say 'wrote ' . colored(['cyan on_bright_yellow'], "$args->{filename}");
	return $args->{filename};
}

sub seq_type ($array, $start_index = 1) {
=purpose
from CIAlign:
    Detects if an alignment is of nucleotides or amino acids using pre-built
    dictionarys of amino acid and nucleotide codes.
    Checks if arr contains characters that are not in the dictionary (not
    IUPAC)

    Parameters
    ----------
    $array
        array reference containing the alignment

    Returns
    -------
    string
    'aa' for amino acid and 'nt for nucleotide

this subroutine takes an array reference like this and an index to start
	[
    [0] ">NP_001150334.1 uncharacterized protein LOC100283964 [Zea mays]",
    [1] "MASMLAILRPSAPAPLAGRRARAAAPATARVALSSRSRYSSVRVSLGSEVAVGADALFADYKPTTAFLFPGQGAQTVGMG",
    [2] "AEAQSVPAATKLFNQANEILGYDLLDLCTNGPKEKLDSTMISQPAIYVTSLAAVEVLRARDGGQDVINSVDVTCGLSLGE",
    [3] "YTALAFAGAFSFEDGLKLVKLRGEAMQDASDAANSAMVSVIGLDSEKVQELCDAANDEVDENDRVQIANFLCPGNYAVSG",
    [4] "GVKGIEVVEAKAKSFKARMTVRLAVAGAFHTSFMQPAVSRLESALAATEIRTPRIPVISNVDAQPHSDPNTIKQILAQQV",
    [5] "TSPVQWETTVKNLMGKGLEISYELGPGKVIAGILKRINKGTSIENIGA"
] (8K)

=cut
=my ($args) = @_;
	unless (ref $args eq 'HASH') {
		my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
		die "args must be given as a hash ref, e.g. \"$current_sub({ filename => 'blah.xlsx' })\"
";
	}
	my @reqd_args = ('sequence');
	my @undef_args = grep { !defined $args->{$_}} @reqd_args;
	if (scalar @undef_args > 0) {
		p @undef_args;
		die 'the above args are necessary, but were not defined.';
	}
	my @defined_args = ('die', 'start_index', @reqd_args);
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		say 'the above arguments are not recognized.';
		p @defined_args;
		die 'The above args are accepted.'
	}
=cut
#	my $current_sub = (split(/::/,(caller(0))[3]))[-1];
	if (ref $array eq 'ARRAY') {
		$array = join ('', @{ $array });
	}
#	if (scalar @{ $array } == 0) {
#		die "$current_sub has an empty array";
#	}
	if (
			($start_index != 0) &&
			($start_index != 1)
		) {
		die "\$start_index = $start_index, which can only be 0 or 1";
	}
	my %nt = map {$_ => 1} qw(A C G N T); # "N" is unknown "N"ucleotide
	my %aa = map {$_ => 1} qw(A R N D C E Q G H I L K M F P S T W Y V);
	my ($nt_count, $aa_count, $total) = (0,0,0);
#	foreach my $r (1..$#$array) { # skip the title line, obviously
		foreach my $letter (split '', $array) {
			next if $letter eq '*';
			$letter = uc $letter;
			if (
					(not defined $aa{$letter})
					&&
					(not defined $nt{$letter})
				) {
				die "$letter is neither a defined nucleic acid, nor a defined amino acid";
			}
			$nt_count++ if defined $nt{$letter};
			$aa_count++ if defined $aa{$letter};
			$total++;
		}
#	}
	return 'nt' if $total == $nt_count;
	return 'aa' if $total == $aa_count;
	return 'other';# if $args->{'die'} > 0;
}
sub write_2d_array_to_tex_tabular {
	my ($args) = @_; # needs tex_filename, header, & 2D array
	my $current_sub = (split(/::/,(caller(0))[3]))[-1]; # https://stackoverflow.com/questions/2559792/how-can-i-get-the-name-of-the-current-subroutine-in-perl
	my @reqd_keys = ('data', 'tex_filename');
	my @undef_keys = grep {!defined $args->{$_}} @reqd_keys;
	if (scalar @undef_keys > 0) {
		p @undef_keys;
		die "the above keys are necessary for $current_sub, and are not defined.";
	}
	my @defined_args = (@reqd_keys,
		'alignment',		# e.g. "\centering"
		'bold.1st.col',	# 
		'check_header',
		'col_alignment',	# "c", "l", or "r", for going into declaring tables; all the same (edit later if desired)
		'comments',
		'format',
		'header', 'json.output', 'max.rows.tex',
		'max.width', 'size', 'worksheet.name', 'xlsx_filename');
	my @bad_args = grep { my $key = $_; not grep {$_ eq $key} @defined_args} keys %{ $args };
	if (scalar @bad_args > 0) {
		p @bad_args;
		say "the above arguments are not recognized by $current_sub.";
		p @defined_args;
		die "The above args are accepted by $current_sub"
	}
	if (grep {not defined $args->{$_}} keys %{ $args }) {
		p $args;
		die "undefined arguments were given to $current_sub";
	}
	if (scalar @{ $args->{data} } == 0) {
		die "An empty array of data was given to $current_sub";
	}
	if (not defined $args->{header}) {
		$args->{header} = shift @{ $args->{data} };
	}
	if (defined $args->{'json.output'}) {
		ref_to_json_file( $args->{data}, $args->{'json.output'} );
	}
	my $no_header_cols = scalar @{ $args->{header} };
	$args->{format} = $args->{format} // 0;
	$args->{tex_filename} =~ s/[\s:]+/_/g; # get rid of annoying characters in filenames
	my ($workbook, $bold, $worksheet);
	open my $tex, '>', $args->{tex_filename};
	$workbook = Excel::Writer::XLSX -> new($args->{xlsx_filename}) if defined $args->{xlsx_filename};
	$bold = $workbook->add_format(bold => 1) if defined $args->{xlsx_filename};
	$workbook -> set_properties(comments => comments()) if defined $args->{xlsx_filename};
	$args->{'worksheet.name'} = $args->{'worksheet.name'} // 'worksheet A';
	$worksheet = $workbook -> add_worksheet($args->{'worksheet.name'}) if defined $args->{xlsx_filename};
	$worksheet -> freeze_panes(1,0) if defined $args->{xlsx_filename}; # freeze header
	my @header = @{ $args->{header} };
	$worksheet -> write_row(0, 0, \@header, $bold) if defined $args->{xlsx_filename};
	if (defined $args->{header}) {
		foreach my $col (@{ $args->{header} }) {
			$col =~ s/([#_%&])/\\$1/g;
		}
	}
#----
	say $tex '%written by ' . getcwd() . '/' . $RealScript;
	if (defined $args->{comments}) {
		if ((ref $args->{comments}) eq '') {
			say $tex '%' . $args->{comments};
		} elsif (ref $args->{comments} eq 'ARRAY') {
			foreach my $c (@{ $args->{comments} }) {
				say $tex '%' . $c;
			}
		}
	}
	$args->{col_alignment} = $args->{col_alignment} // 'c';
	say $tex '\begin{tabular}{|' . ("$args->{col_alignment}|" x $no_header_cols) . '} \hline';
	say $tex $args->{size} if defined $args->{size};
	my $latex_header = join ('} & \textbf{', @{ $args->{header} });
	say $tex  '\textbf{' . $latex_header . '} \\\\ \hline';
	my $row = 1;
	$args->{check_header} = $args->{check_header} // 1;
	$args->{'max.rows.tex'} = $args->{'max.rows.tex'} // 'inf';
	$args->{'max.width'} = $args->{'max.width'} // 'inf';
	$args->{'bold.1st.col'} = $args->{'bold.1st.col'} // 'yes';
	foreach my $line (@{ $args->{data} }) {
		if (
				(scalar @{ $line } == 1)
				&&
				($line->[0] eq '\hline')
			) {
			say $tex '\hline';
			next;
		}
		my @line;
		if (($args->{check_header} > 0) && (scalar @{ $line } != $no_header_cols)) {
			say 'this line has ' . scalar @{ $line } . ' columns';
			p $line, array_max => scalar @{ $line };
			say 'but the header has ' . scalar @{ $args->{header} } . ' columns';
			p $args->{header}, array_max => scalar @{ $args->{header } };
			die "# of columns does not match " . join (',', @{ $args->{header} });
		}
		if (defined $args->{xlsx_filename}) {
			my @row = @{ $line };
			$worksheet->write($row, 0, $row[0], $bold);
			shift @row;
			$worksheet->write_row($row, 1, \@row);
		}
		$row++;
		next if $row > $args->{'max.rows.tex'}; # but keep writing xlsx
		foreach my $e (@{ $line }) {
			if (not defined $e) {
				p $args;
				die "undefined value found for $line";
			}
			my $element = $e;
			$element =~ s/([#_%&])/\\$1/g;
			if ($args->{format} != 0 && looks_like_number($element)) {
				$element = sprintf("%.4g", $element);
			}
			push @line, widecell2tabular($element, $args->{'max.width'}, $args->{col_alignment});
		}
		$line[0] = '\textbf{' . "$line[0]}" unless $args->{'bold.1st.col'} eq 'no';
		say $tex join (' & ', @line) . '\\\\';
	}
	if (
			(scalar @{ $args->{data}[-1] } == 1)
			&& # does the data end in an "\hline"?
			($args->{data}[-1][0] eq '\hline')
		) {
		say $tex '\end{tabular}';
	} else {
		say $tex '\hline \end{tabular}';
	}
	if (defined $args->{comment}) {
		say $tex '%' . $args->{comment};
	}
	close $tex;
	say 'wrote ' . colored(['black on_cyan'], $args->{tex_filename});
	say 'wrote ' . colored(['black on_bright_red'], $args->{xlsx_filename}) if defined $args->{xlsx_filename};
	return $args->{tex_filename};
}
=example
write_2d_array_to_tex_tabular({
	xlsx_filename  => 'xlsx/blast.non.yeast.hits.xlsx',
	data           => \@table,
	tex_filename   => 'tex/blast.non.yeast.hits.tex',
	header         => ['Protein', 'Hits', 'Name'],
	'max.rows.tex' => 20,
	'max.width'    => 90
});
=cut
1;
