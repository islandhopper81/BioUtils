package BioUtils::SeqSet::Diagnostics;

# This module can be used to find out more information on DNA sequences
# (especially) those with IUPAC characters.

use warnings;
use strict;
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use File::Temp qw(tempfile);
use Chart::Graph::Gnuplot qw(gnuplot);
use BioUtils::FastqIO 1.0.1;
use MyX::Generic 1.0.1;
use version; our $VERSION = qv('1.0.1');


{
	Readonly::Hash my %AMBIG_BASES => (
		R => 1,
		Y => 1,
		S => 1,
		W => 1,
		K => 1,
		M => 1,
		B => 1,
		D => 1,
		H => 1,
		V => 1,
		N => 1,
	);
	
	
    # Attributes #
    my %gap_chars_of;
	my %gap_count_of;
	my %gap_pos_dist_of;
	my %gaps_per_seq_of;
	my %ambig_base_count_of;
	my %ambig_base_pos_dist_of;
	my %ambig_base_comp_dist_of;
	my %ambig_bases_per_seq_of;
    
    # Setters #
    
    # Getters #
	sub get_gap_count;
	sub get_gap_pos_dist;
	sub get_gap_pos_dist_str;
	sub get_gaps_per_seq;
	sub get_gaps_per_seq_str;
	sub get_ambig_base_count;
	sub get_ambig_base_pos_dist;
	sub get_ambig_base_pos_dist_str;
	sub get_ambig_base_comp_dist;
	sub get_ambig_base_comp_dist_str;
	sub get_ambig_bases_per_seq;
	sub get_ambig_bases_per_seq_str;
    
    # Others #
	sub _init;
    sub run_diagnosis;
    sub _diagnose_seqeunce;
	sub _add_gap_at_pos;
	sub _add_ambig_base_at_pos;
    sub add_gap_character;
    sub remove_gap_character;
	sub has_gap_character;
	sub _increment_hash_value;
	
	# Output #
	sub print_gap_pos_dist_str;
	sub graph_gap_pos_dist;
	sub print_gaps_per_seq_str;
	sub graph_gaps_per_seq;
	sub print_ambig_base_pos_dist_str;
	sub graph_ambig_base_pos_dist;
	sub print_ambig_base_comp_dist_str;
	sub graph_ambig_base_comp_dist;
	sub print_ambig_bases_per_seq_str;
	sub graph_ambig_bases_per_seq;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # This is a hash that contains all the valid gap characters
        # Set the DEFUALT gap characters to '-'
        $gap_chars_of{ident $new_obj} = {'-' => 1};
		
		# Other attributes
		$new_obj->_init();
        
        return $new_obj;
    }
	
	sub _init {
		my ($self) = @_;
		
		$gap_count_of{ident $self} = 0;
		$gap_pos_dist_of{ident $self} = {};
		$gaps_per_seq_of{ident $self} = {};
		$ambig_base_count_of{ident $self} = 0;
		$ambig_base_pos_dist_of{ident $self} = {};
		$ambig_base_comp_dist_of{ident $self} = {};
		$ambig_bases_per_seq_of{ident $self} = {};
		
		return 1;
	}
    
    ###########
    # Setters #
    ###########
    
    
    ###########
    # Getters #
    ###########
    sub get_gap_count {
		my ($self) = @_;
		return $gap_count_of{ident $self};
	}
	
	sub get_gap_pos_dist {
		my ($self) = @_;
		return $gap_pos_dist_of{ident $self};
	}
	
	sub get_gap_pos_dist_str {
		my ($self) = @_;
		
		my $str = "";  # An empty string
		my $href = $gap_pos_dist_of{ident $self};
		foreach my $key ( sort {$a<=>$b} keys %{$href} ) {  # sort numerically 
			$str .= $key . "\t" . $href->{$key} . "\n";
		}
		
		return $str;
	}
	
	sub get_gaps_per_seq {
		my ($self) = @_;
		return $gaps_per_seq_of{ident $self};
	}
	
	sub get_gaps_per_seq_str {
		my ($self) = @_;
		
		my $str = "";  # An empty string
		my $href = $gaps_per_seq_of{ident $self};
		foreach my $key ( sort {$a<=>$b} keys %{$href} ) {  # sort numerically 
			$str .= $key . "\t" . $href->{$key} . "\n";
		}
		
		return $str;
	}
	
	sub get_ambig_base_count {
		my ($self) = @_;
		return $ambig_base_count_of{ident $self};
	}
	
	sub get_ambig_base_pos_dist {
		my ($self) = @_;
		return $ambig_base_pos_dist_of{ident $self};
	}
	
	sub get_ambig_base_pos_dist_str {
		my ($self) = @_;
		
		my $str = "";  # An empty string
		my $href = $ambig_base_pos_dist_of{ident $self};
		foreach my $key ( sort {$a<=>$b} keys %{$href} ) {  # sort numerically 
			$str .= $key . "\t" . $href->{$key} . "\n";
		}
		
		return $str;
	}
	
	sub get_ambig_base_comp_dist {
		my ($self) = @_;
		return $ambig_base_comp_dist_of{ident $self};
	}
	
	sub get_ambig_base_comp_dist_str {
		my ($self) = @_;
		
		my $str = "";  # An empty string
		my $href = $ambig_base_comp_dist_of{ident $self};
		foreach my $key ( sort keys %{$href} ) {  # sort alphabetically 
			$str .= $key . "\t" . $href->{$key} . "\n";
		}
		
		return $str;
	}
	
	sub get_ambig_bases_per_seq {
		my ($self) = @_;
		return $ambig_bases_per_seq_of{ident $self};
	}
	
	sub get_ambig_bases_per_seq_str {
		my ($self) = @_;
		
		my $str = "";  # An empty string
		my $href = $ambig_bases_per_seq_of{ident $self};
		foreach my $key ( sort {$a<=>$b} keys %{$href} ) {  # sort numerically 
			$str .= $key . "\t" . $href->{$key} . "\n";
		}
		
		return $str;
	}
	
    
    ##########
    # Others #
    ##########
    sub run_diagnosis {
        my ($self, $input_file) = @_;
        
        # Open the sequence file
        my $seqs_in;
        eval {
            $seqs_in = BioUtils::FastqIO->new( {
							stream_type => '<',
							file => $input_file,
						});
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            croak("$e->error()\n$e->trace()\n");
        }
        elsif ( $@ ) {
            croak($@);
        }
        
        # Loop through each sequence
        while ( my $seq_obj = $seqs_in->get_next_seq() ) {
            $self->diagnose_sequence($seq_obj->get_seq());
        }
		
		return 1;
    }
    
    sub diagnose_sequence {
        my ($self, $seq_str) = @_;
        
        # Go through each character in the sequence string
		my $pos = 0;
		my $gap_count = 0;
		my $ambig_base_count = 0;
        foreach my $char ( split //, $seq_str ) {
			$char = uc $char;
			
			# if the character is a gap
            if ( defined $gap_chars_of{ident $self}->{$char} ) {
                $self->_add_gap_at_pos($pos);
				$gap_count++;
            }
			
			# if the character is an IUPAC ambiguous base
			if ( defined $AMBIG_BASES{$char} ) {
				$self->_add_ambig_base_at_pos($char, $pos);
				$ambig_base_count++;
			}
			
            $pos++;
        }
		
		# Increment the gaps_per_seq and ambig_bases_per_seq distributions
		if ($gap_count > 0 ) {
			_increment_hash_value($gaps_per_seq_of{ident $self},
								  $gap_count);
		}
		if ( $ambig_base_count > 0 ) {
			_increment_hash_value($ambig_bases_per_seq_of{ident $self},
								  $ambig_base_count);
		}
    }
    
	sub _add_gap_at_pos {
		my ($self, $pos) = @_;
		
		# Increment the gap count
		$gap_count_of{ident $self}++;
		
		# Add the gap position to the gap_pos_dist
		_increment_hash_value($gap_pos_dist_of{ident $self}, $pos);
		
		return 1;
	}
	
	sub _add_ambig_base_at_pos {
		my ($self, $char, $pos) = @_;
		
		# Increment the ambiguous base count
		$ambig_base_count_of{ident $self}++;
		
		# Add the ambiguous base to the ambig_base_pos_dist
		_increment_hash_value($ambig_base_pos_dist_of{ident $self}, $pos);
		
		# Add the ambiguous base to the ambig_base_comp_dist
		_increment_hash_value($ambig_base_comp_dist_of{ident $self}, $char);
		
		return 1;
	}
	
    sub add_gap_character {
        my ($self, $char) = @_;
        
        # Do some error checks
        if ( length $char != 1 ) {
            croak("Gap character must be ONE character");
        }
        if ( $char !~ m/\S/ ) {
            croak("Gap character cannot be a white spacer");
        }
        
		# If the gap character was empty tell the user that it is now populated
		if (! keys %{$gap_chars_of{ident $self}} ) {
			carp("Gap character no longer empty after adding: $char");
		}
		
        # add the gap character
        $gap_chars_of{ident $self}->{$char} = 1;
        
        return 1;
    }
    
    sub remove_gap_character {
        my ($self, $char) = @_;
        
        # Do some error checks
        if ( length $char != 1 ) {
            croak("Gap character must be ONE character");
        }
        if ( $char !~ m/\S/ ) {
            croak("Gap character cannot be a white spacer");
        }
        
        # remove the gap character if it exists
        if ( $self->has_gap_character($char) ) {
            delete $gap_chars_of{ident $self}->{$char};
			
			# make sure a gap character is set in the hash
			if (! keys %{$gap_chars_of{ident $self}} ) {
				carp("All gap characters have been removed");
			}
        }
        
        return 1;
    }
	
	sub has_gap_character {
		my ($self, $char) = @_;
		
        # Do some error checks
        if ( length $char != 1 ) {
            croak("Gap character must be ONE character");
        }
        if ( $char !~ m/\S/ ) {
            croak("Gap character cannot be a white spacer");
        }
		
		if ( $gap_chars_of{ident $self}->{$char} ) {
			return 1; # TRUE
		}
		
		return 0;  # FALSE
	}
	
	sub _increment_hash_value {
		my ($href, $key) = @_;
		
		if ( defined $href->{$key} ) {
			$href->{$key}++;
		}
		else {
			$href->{$key} = 1;
		}
		
		return 1;
	}
	
	
	##########
	# Output #
	##########
	sub print_gap_pos_dist_str {
		my ($self, $file) = @_;
		
		if ( ! -f $file ) {
			croak("File: $file is not a file");
		}
		
		open my $OUT, ">$file" or croak("$file cannot be open\nERROR: $!\n");
		
		print $OUT $self->get_gap_pos_dist_str();
		
		close($OUT);
		
		return 1;
	}
	
	sub graph_gap_pos_dist {
		my ($self, $file_prefix) = @_;
		
		# print the gap_pos_dist to a temp file
		my ($fh, $temp_file) = tempfile();
		$self->print_gap_pos_dist_str($temp_file);
		
		# graph with gnuplot
		gnuplot(
				{'title' => 'Gaps per Position',
				'output type' => 'png',
				'output file' => "$file_prefix.png",
				'x-axis label' => "Position",
				'y-axis label' => "Count",
				'xrange' => "[0:]",
				'yrange' => "[0:]",
				'extra_opts' => "set key left top
								 set style fill solid 0.5 
								 set tics out nomirror
								 binwidth = 1
								 set boxwidth binwidth*0.9
								",
				},
				[{
				  'style' => 'boxes',
				  'title' => '',
				  'type' => 'file',
				  'using' => '1:2'
				  },
				 "$temp_file"
				],
			);
		
		return 1;
	}
	
	sub print_gaps_per_seq_str {
		my ($self, $file) = @_;
		
		if ( ! -f $file ) {
			croak("File: $file is not a file");
		}
		
		open my $OUT, ">$file" or croak("$file cannot be open\nERROR: $!\n");
		
		print $OUT $self->get_gaps_per_seq_str();
		
		close($OUT);
		
		return 1;
	}
	
	sub graph_gaps_per_seq {
		my ($self, $file_prefix) = @_;
		
		# print the gap_pos_dist to a temp file
		my ($fh, $temp_file) = tempfile();
		$self->print_gaps_per_seq_str($temp_file);
		
		# graph with gnuplot
		gnuplot(
				{'title' => 'Gaps per Sequence',
				'output type' => 'png',
				'output file' => "$file_prefix.png",
				'x-axis label' => "Number of Gaps",
				'y-axis label' => "Number of Seqs",
				'xrange' => "[0:]",
				'yrange' => "[0:]",
				'extra_opts' => "set key left top
								 set style fill solid 0.5 
								 set tics out nomirror
								 binwidth = 1
								 set boxwidth binwidth*0.9
								",
				},
				[{
				  'style' => 'boxes',
				  'title' => '',
				  'type' => 'file',
				  'using' => '1:2'
				  },
				 "$temp_file"
				],
			);
		
		return 1;
	}
	
	sub print_ambig_base_pos_dist_str {
		my ($self, $file) = @_;
		
		if ( ! -f $file ) {
			croak("File: $file is not a file");
		}
		
		open my $OUT, ">$file" or croak("$file cannot be open\nERROR: $!\n");
		
		print $OUT $self->get_ambig_base_pos_dist_str();
		
		close($OUT);
		
		return 1;
	}
	
	sub graph_ambig_base_pos_dist {
		my ($self, $file_prefix) = @_;
		
		# print the gap_pos_dist to a temp file
		my ($fh, $temp_file) = tempfile();
		$self->print_ambig_base_pos_dist_str($temp_file);
		
		# graph with gnuplot
		gnuplot(
				{'title' => 'Ambiguous Bases per Position',
				'output type' => 'png',
				'output file' => "$file_prefix.png",
				'x-axis label' => "Position",
				'y-axis label' => "Count",
				'xrange' => "[0:]",
				'yrange' => "[0:]",
				'extra_opts' => "set key left top
								 set style fill solid 0.5 
								 set tics out nomirror
								 binwidth = 1
								 set boxwidth binwidth*0.9
								",
				},
				[{
				  'style' => 'boxes',
				  'title' => '',
				  'type' => 'file',
				  'using' => '1:2'
				  },
				 "$temp_file"
				],
			);
		
		return 1;
	}
	
	sub print_ambig_base_comp_dist_str {
		my ($self, $file) = @_;
		
		if ( ! -f $file ) {
			croak("File: $file is not a file");
		}
		
		open my $OUT, ">$file" or croak("$file cannot be open\nERROR: $!\n");
		
		print $OUT $self->get_ambig_base_comp_dist_str();
		
		close($OUT);
		
		return 1;
	}
	
	sub graph_ambig_base_comp_dist {
		my ($self, $file_prefix) = @_;
		
		# print the gap_pos_dist to a temp file
		my ($fh, $temp_file) = tempfile();
		$self->print_ambig_base_comp_dist_str($temp_file);
		
		# graph with gnuplot
		gnuplot(
				{'title' => 'Ambiguous Bases Composition',
				'output type' => 'png',
				'output file' => "$file_prefix.png",
				'x-axis label' => "Position",
				'y-axis label' => "Count",
				'xrange' => "[0:]",
				'yrange' => "[0:]",
				'extra_opts' => "set key left top
								 set style fill solid 0.5 
								 set tics out nomirror
								 binwidth = 1
								 set boxwidth binwidth*0.9
								",
				},
				[{
				  'style' => 'boxes',
				  'title' => '',
				  'type' => 'file',
				  'using' => '2:xtic(1)'
				  },
				 "$temp_file"
				],
			);
		
		return 1;
	}
	
	sub print_ambig_bases_per_seq_str {
		my ($self, $file) = @_;
		
		if ( ! -f $file ) {
			croak("File: $file is not a file");
		}
		
		open my $OUT, ">$file" or croak("$file cannot be open\nERROR: $!\n");
		
		print $OUT $self->get_ambig_bases_per_seq_str();
		
		close($OUT);
		
		return 1;
	}
	
	sub graph_ambig_bases_per_seq {
		my ($self, $file_prefix) = @_;
		
		# print the gap_pos_dist to a temp file
		my ($fh, $temp_file) = tempfile();
		$self->print_ambig_bases_per_seq_str($temp_file);
		
		# graph with gnuplot
		gnuplot(
				{'title' => 'Ambiguous Bases per Sequence',
				'output type' => 'png',
				'output file' => "$file_prefix.png",
				'x-axis label' => "Number of Ambiguous Bases",
				'y-axis label' => "Number of Seqs",
				'xrange' => "[0:]",
				'yrange' => "[0:]",
				'extra_opts' => "set key left top
								 set style fill solid 0.5 
								 set tics out nomirror
								 binwidth = 1
								 set boxwidth binwidth*0.9
								",
				},
				[{
				  'style' => 'boxes',
				  'title' => '',
				  'type' => 'file',
				  'using' => '1:2'
				  },
				 "$temp_file"
				],
			);
		
		return 1;
	}
}



1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::SeqSet::Diagnostics - A module for accessing information about
DNA sequences


=head1 VERSION

This document describes BioUtils::SeqSet::Diagnostics version 1.0.1


=head1 SYNOPSIS

    use BioUtils::SeqSet::Diagnostics 1.0.1;
	
	my $my_diag_obj = BioUtils::SeqSet::Diagnostics->new();
	
	# Add/Remove a gap character [DEFAULT: -]
	$my_diag_obj->add_gap_char(".");
	$my_diag_obj->revmove_gap_char('-');
	
	# This runs the diagnosis.  The results are saved in the objects attributes
	$my_diag_obj->run_diagnosis($seqs_file);
	
	# Run the diagnosis on one sequence
	$my_diag_obj->diagnose_sequence("ATTG");
	
	# getters
	my $gap_count = $my_diag_obj->get_gap_count();
	my $gap_pos_dist_href = $my_diag_obj->get_gap_pos_dist();
	my $gap_pos_dist_str = $my_diag_obj->get_gap_pos_dist_str();
	my $gaps_per_seq_href = $my_diag_obj->get_gaps_per_seq();
	my $gaps_per_seq_str = $my_diag_obj->get_gaps_per_seq_str();
	my $ambig_base_count = $my_diag_obj->get_ambig_base_count();
	my $ambig_base_pos_dist_href = $my_diag_obj->get_ambig_base_pos_dist();
	my $ambig_base_copm_dist_href = $my_diag_obj->get_ambig_base_comp_dist();
	my $ambig_base_pos_dist_str = $my_diag_obj->get_ambig_base_pos_dist_str();
	my $ambig_base_comp_dist_str = $my_diag_obj->get_ambig_base_comp_dist_str();
	my $ambig_bases_per_seq_href = $my_diag_obj->get_ambig_bases_per_seq();
	my $ambig_bases_per_seq_str = $my_diag_obj->get_ambig_bases_per_seq_str();
	
	# output methods
	$my_diag_obj->print_gap_pos_dist_str($gap_pos_file);
	$my_diag_obj->print_gaps_per_seq_str($gas_per_seq_file);
	$my_diag_obj->print_ambig_base_pos_dist_str($ambig_base_pos_file);
	$my_diag_obj->print_ambig_base_comp_dist_str($ambig_base_comp_file);
	$my_diag_obj->print_ambig_bases_per_seq_str($ambig_bases_per_seq_file);
	
	# print graphs using gnuplot -- in png format
	$my_diag_obj->graph_gap_pos_dist($file_prefix);
	$my_diag_obj->graph_gaps_per_seq($file_prefix);
	$my_diag_obj->graph_ambig_base_pos_dist($file_prefix);
	$my_diag_obj->graph_ambig_base_comp_dist($file_prefix);
	$my_diag_obj->graph_ambig_bases_per_seq($file_prefix);
  
=head1 DESCRIPTION

This object is used to analysis the diagnostics of a set of sequences especially
sequences that contain gaps and ambiguous bases.  


=head1 CONFIGURATION AND ENVIRONMENT

gnuplot must be installed and locatable by Chart::Graph::Gnuplot qw(gnuplot)


=head1 DEPENDENCIES

	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	File::Temp qw(tempfile)
	Chart::Graph::Gnuplot qw(gnuplot)
	BioUtils::FastqIO 1.0.1
	MyX::Generic 1.0.1
	version


=head1 INCOMPATIBILITIES

	None reported.


=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: my $my_diag_obj = BioUtils::SeqSet::Diagnostics->new();
	Function: Create a new BioUtils::SeqSet::Diagnostics object
	Returns: Reference to a SeqeunceSetDiagnostics object
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _init

	Title: _init
	Usage: $my_diag_obj->_init()
	Function: Initializes/resets the attribute values
	Returns: NA
	Args: NA
	Throws: NA
	Comments: It is recommended to make a new object instead of using _init to
			  reset
	See Also: NA

=head2 get_gap_count

	Title: get_gap_count
	Usage: my $gap_count = $my_diag_obj->get_gap_count();
	Function: Returns the number of gaps in the sequence set
	Returns: int
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_gap_pos_dist

	Title: get_gap_pos_dist
	Usage: my $gap_pos_dist_href = $my_diag_obj->get_gap_pos_dist()
	Function: Returns the gaps per position in the sequence
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_gap_pos_dist_str

	Title: get_gap_pos_dist_str
	Usage: my $gap_pos_dist_str = $my_diag_obj->get_gap_pos_dist_str()
	Function: Returns the distribution as a string with KEY[tab]VALUE
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_gaps_per_seq

	Title: get_gaps_per_seq
	Usage: my $gaps_per_seq_href = $my_diag_obj->get_gaps_per_seq()
	Function: Returns the gaps per sequence for the set
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_gaps_per_seq_str

	Title: get_gaps_per_seq_str
	Usage: my $gap_per_seq_str = $my_diag_obj->get_gaps_per_seq_str()
	Function: Returns the distribution as a string with KEY[tab]VALUE
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_ambig_base_count

	Title: get_ambig_base_count
	Usage: my $ambig_base_count = $my_diag_obj->get_ambig_base_count()
	Function: Returns the number of ambiguous bases in the sequence set
	Returns: int
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_ambig_base_pos_dist

	Title: get_ambig_base_pos_dist
	Usage: my $ambig_base_pos_dist_href = $my_diag_obj->get_ambig_base_pos_dist()
	Function: Returns the ambiguous bases per position in the sequence
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_ambig_base_pos_dist_str

	Title: get_ambig_base_pos_dist_str
	Usage: my $ambig_base_pos_dist_str = $my_diag_obj->get_ambig_base_pos_dist_str()
	Function: Returns the distribution as a string with KEY[tab]VALUE
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_ambig_base_comp_dist

	Title: get_ambig_base_comp_dist
	Usage: my $ambig_base_comp_dist_href = $my_diag_obj->get_ambig_base_comp_dist()
	Function: Returns the ambiguous bases by composition that are in the sequence
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: For example, this will tell you there are 15 N's in the set
	See Also: NA

=head2 get_ambig_base_comp_dist_str

	Title: get_ambig_base_comp_dist_str
	Usage: my $ambig_base_comp_dist_str = $my_diag_obj->get_ambig_base_comp_dist_str()
	Function: Returns the distribution as a string with KEY[tab]VALUE
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_ambig_bases_per_seq

	Title: get_ambig_bases_per_seq
	Usage: my $ambig_bases_per_seq_href = $my_diag_obj->get_ambig_bases_per_seq()
	Function: Returns the ambiguous bases per sequence for the set
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_ambig_bases_per_seq_str

	Title: get_ambig_bases_per_seq_str
	Usage: my $ambig_bases_per_seq_str = $my_diag_obj->get_ambig_bases_per_seq_str()
	Function: Returns the distribution as a string with KEY[tab]VALUE
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 run_diagnosis

	Title: run_diagnosis
	Usage:  $my_diag_obj->run_diagnosis($input_file)
	Function: Goes through each sequence and stores the necessary info
	Returns: 1 on success
	Args: -input_file => a file with FASTQ formated sequences
	Throws: NA
	Comments: NA
	See Also: NA

=head2 diagnose_sequence

	Title: diagnose_sequence
	Usage:  $my_diag_obj->diagnose_sequence($seq)
	Function: Goes through each base of a sequence and checks for bases of interest
	Returns: 1 on success
	Args: -seq => a sequence as a string
	Throws: NA
	Comments: This method can be used to add sequence information to a
			  BioUtils::SeqSet::Diagnostics object one sequence at a time.
	See Also: NA

=head2 _add_gap_at_pos

	Title: _add_gap_at_pos
	Usage:  $my_diag_obj->_add_gap_at_pos($pos)
	Function: increments the gap count and adds the position to the gap pos dist
	Returns: 1 on success
	Args: -pos => the position of the gap in the sequence
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _add_ambig_base_at_pos

	Title: _add_ambig_base_at_pos
	Usage:  $my_diag_obj->_add_ambig_base_at_pos($base, $pos)
	Function: increments the ambig base count and adds the position to the ambig
			  base to the ambig base pos dist
	Returns: 1 on success
	Args: -base => the ambiguous base
		  -pos => the position of the ambiguous base in the sequence
	Throws: NA
	Comments: NA
	See Also: NA

=head2 add_gap_character

	Title: add_gap_character
	Usage:  $my_diag_obj->add_gap_character($char)
	Function: add a gap character to the set that are already checked for in
			  run_diagnosis()
	Returns: 1 on success
	Args: -char => a character to be added to the possible gap characters
	Throws: Warning if the gap characters goes from being empty
	Comments: I throw the above warning so that if you remove all the gap
			  characters it will be easy to tell that you have begun adding
			  back to their set.
	See Also: NA

=head2 remove_gap_character

	Title: remove_gap_character
	Usage:  $my_diag_obj->remove_gap_character($char)
	Function: removes a gap character from the set that checked for in
			  run_diagnosis()
	Returns: 1 on success
	Args: -char => a character to be removed from the possible gap characters
	Throws: Warning if all character are removed
	Comments: NA
	See Also: NA

=head2 has_gap_character

	Title: has_gap_character
	Usage:  $my_diag_obj->has_gap_character($char)
	Function: checks for the given character of possible gap characters checkd in
			  run_diagnosis()
	Returns: bool
	Args: -char => a character to check for
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _increment_hash_value

	Title: _increment_hash_value
	Usage:  _increment_hash_value($href, $key)
	Function: increments the value in $href at $key
	Returns: 1 on success
	Args: -href => an href
		  -key => the key to increment or add
	Throws: NA
	Comments: NA
	See Also: NA

=head2 print_gap_pos_dist_str

	Title: print_gap_pos_dist_str
	Usage:  $my_diag_obj->print_gap_pos_dist_str($file)
	Function: Prints the gap position distribution to the given file
	Returns: 1 on success
	Args: -file => a file to print the output
	Throws: Dies if given file does not exist
	Comments: NA
	See Also: NA

=head2 graph_gap_pos_dist

	Title: graph_gap_pos_dist
	Usage:  $my_diag_obj->graph_gap_pos_dist($file_prefix)
	Function: Creates a png file with the graphed distribution
	Returns: 1 on success
	Args: -file_prefix => file name prefix (suffix .png will be added
						  automatically)
	Throws: NA
	Comments: Requires gnuplot
	See Also: NA

=head2 print_gaps_per_seq_str

	Title: print_gaps_per_seq_str
	Usage:  $my_diag_obj->print_gaps_per_seq_str($file)
	Function: Prints the gaps per sequence distribution to the given file
	Returns: 1 on success
	Args: -file => a file to print the output
	Throws: Dies if given file does not exist
	Comments: NA
	See Also: NA

=head2 graph_gaps_per_seq

	Title: graph_gaps_per_seq
	Usage:  $my_diag_obj->graph_gaps_per_seq($file_prefix)
	Function: Creates a png file with the graphed distribution
	Returns: 1 on success
	Args: -file_prefix => file name prefix (suffix .png will be added
						  automatically)
	Throws: NA
	Comments: Requires gnuplot
	See Also: NA

=head2 print_ambig_base_pos_dist_str

	Title: print_ambig_base_pos_dist_str
	Usage:  $my_diag_obj->print_ambig_base_pos_dist_str($file)
	Function: Prints the ambiguous base position distribution to the given file
	Returns: 1 on success
	Args: -file => a file to print the output
	Throws: Dies if given file does not exist
	Comments: NA
	See Also: NA

=head2 graph_ambig_base_pos_dist

	Title: graph_ambig_base_pos_dist
	Usage:  $my_diag_obj->graph_ambig_base_pos_dist($file_prefix)
	Function: Creates a png file with the graphed distribution
	Returns: 1 on success
	Args: -file_prefix => file name prefix (suffix .png will be added
						  automatically)
	Throws: NA
	Comments: Requires gnuplot
	See Also: NA

=head2 print_ambig_base_comp_dist_str

	Title: print_ambig_base_comp_dist_str
	Usage:  $my_diag_obj->print_ambig_base_comp_dist_str($file)
	Function: Prints the ambiguous base composition distribution to the given file
	Returns: 1 on success
	Args: -file => a file to print the output
	Throws: Dies if given file does not exist
	Comments: NA
	See Also: NA

=head2 graph_ambig_base_comp_dist

	Title: graph_ambig_base_comp_dist
	Usage:  $my_diag_obj->graph_ambig_base_comp_dist($file_prefix)
	Function: Creates a png file with the graphed distribution
	Returns: 1 on success
	Args: -file_prefix => file name prefix (suffix .png will be added
						  automatically)
	Throws: NA
	Comments: Requires gnuplot
	See Also: NA

=head2 print_ambig_bases_per_seq_str

	Title: print_ambig_bases_per_seq_str
	Usage:  $my_diag_obj->print_ambig_bases_per_seq_str($file)
	Function: Prints the ambiguous base per sequence for the set
	Returns: 1 on success
	Args: -file => a file to print the output
	Throws: Dies if given file does not exist
	Comments: NA
	See Also: NA

=head2 graph_ambig_bases_per_seq

	Title: graph_ambig_bases_per_seq
	Usage:  $my_diag_obj->graph_ambig_bases_per_seq($file_prefix)
	Function: Creates a png file with the graphed distribution
	Returns: 1 on success
	Args: -file_prefix => file name prefix (suffix .png will be added
						  automatically)
	Throws: NA
	Comments: Requires gnuplot
	See Also: NA

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to Scott Yourstone
C<< <scott.yourstone81@gmail.com> >>


=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2012, Scott Yourstone C<< <scott.yourstone81@gmail.com> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut
