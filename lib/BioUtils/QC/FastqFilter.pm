package BioUtils::QC::FastqFilter;

use warnings;
use strict;

use version; our $VERSION = qv('1.0.1');

use Class::Std::Utils;
use Scalar::Util qw(looks_like_number);
use Getopt::Long;
use Carp qw(cluck carp croak);
use Pod::Usage;
use Statistics::Descriptive;
use Data::Dumper qw(Dumper);
use Readonly;
use Cwd;
use File::Basename;
use BioUtils::FastqIO 1.0.1;
use BioUtils::Codec::QualityScores qw(illumina_1_8_to_int);
use MyX::Generic;


{
    Readonly my $NEW_USAGE => q{ new( {fastq_file => ,
            output_dir => ,
            [min_len => ],
            [min_avg_qual => ],
            [min_base_qual => ],
            [allow_gaps => ],
            [allow_ambig_bases => ],
            [verbose => ],
            } ) };
    
    
    # Attributes #
    my %fastq_file_of;
    my %output_dir_of;
    my %min_len_of;
    my %min_avg_qual_of;
    my %min_base_qual_of;
    my %allow_gaps_of;
    my %allow_ambig_bases_of;
    my %verbose_of;

    # Setters #
    sub set_fastq_file;
    sub set_output_dir;
    sub set_min_len;
    sub set_min_avg_qual;
    sub set_min_base_qual;
    sub set_allow_gaps;
    sub set_allow_ambig_bases;
    sub set_verbose;
    
    # Getters #
    sub get_fastq_file;
    sub get_output_dir;
    sub get_min_len;
    sub get_min_avg_qual;
    sub get_min_base_qual;
    sub get_allow_gaps;
    sub get_allow_ambig_bases;
    sub get_verbose;
    
    # Others #
    sub filter;
    sub _init;
    sub _too_short;
    sub _below_avg_qual;
    sub _below_base_qual;
    sub _has_gaps;
    sub _has_ambig_bases;
    sub _get_file_prefix;
    sub _print_diagnostics;
    sub _Y_N_translate;

    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak('Constructor called on existing object instead of class')
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Initialize the objects attributes
        $new_obj->_init($arg_href);
        
        return $new_obj;
    }
    
    sub filter {
        my ($self) = @_;
        
        my $out_dir = $self->get_output_dir();
        my $input_file = $self->get_fastq_file();
        my $allow_gaps = $self->get_allow_gaps();
        my $allow_ambig_bases = $self->get_allow_ambig_bases();
        my $verbose = $self->get_verbose();
        
        # Get the input file prefix so I can build the output file(s).
        my $file_prefix = _get_file_prefix($input_file);
    
        
        # Build the FASTQ file in/out objects
        my ($in, $out, $filtered_out);
        
        eval {
            $in = BioUtils::FastqIO->new( {
                        stream_type => '<',
                        file => $input_file
                    });
            $out = BioUtils::FastqIO->new( {
                        stream_type => '>',
                        file => $out_dir . "/" .
                                $file_prefix .
                                '_unfiltered.fastq'
                    });
            if ( $verbose ) {
                $filtered_out = BioUtils::FastqIO->new( {
                                    stream_type => '>',
                                    file => $out_dir . "/" .
                                            $file_prefix .
                                            '_filtered.fastq'
                                });
            }
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            croak($e->message() . "\n" .
                  $e->trace() . "\n" .
                  "USAGE: " . $e->usage() . "\n"
                  );
        }
        elsif ( $@ ) {
            print $@;
        }
        
        # a table of why each sequence is filterted out (i.e. diagnostics)
        my @diag_table = ();
        
        # Read in the sequences from $input_file
        # AND run each filter test
        SEQ: while ( my $fastq_seq = $in->get_next_seq() ) {
            
            # A boolean to keep track if the sequence needs to be filtered out
            my $filter_flag = 0;
            
            # A data structure for storing this (a single) sequence's diagnostics
            my @diagnostics = ();
            push @diagnostics, $fastq_seq->get_id();
            
            # save the sequence and quals for testing
            my $seq = $fastq_seq->get_seq();
            my $quals_str = $fastq_seq->get_quals_str();
            
            # Decode the quals into ints and store in an aref
            my $int_quals_aref = illumina_1_8_to_int($quals_str);
            #print Dumper($int_quals_aref) . "\n";
            
            # Create a statistics descriptive object for qual score testing
            my $stat = Statistics::Descriptive::Sparse->new();
            $stat->add_data(@{$int_quals_aref});
            
            # Test length
            if ( _too_short($self->get_min_len(), $seq) ) {
                next SEQ if ( ! $verbose );
                
                # verbose operations
                push @diagnostics, 1;
                $filter_flag = 1;
            }
            else {
                push @diagnostics, 0;
            }
            
            # Test minimum average quality over whole read
            if ( _below_avg_qual($self->get_min_avg_qual(), $stat->mean()) ) {
                next SEQ if ( ! $verbose );
                
                # verbose operations
                push @diagnostics, 1;
                $filter_flag = 1;
            }
            else {
                push @diagnostics, 0;
            }
            
            # Test minimum base quality
            if ( _below_min_base_qual(
                                      $self->get_min_base_qual(),
                                      $stat->min())
                ) {
                next SEQ if ( ! $verbose );
                
                # verbose operations
                push @diagnostics, 1;
                $filter_flag = 1;
            }
            else {
                push @diagnostics, 0;
            }
            
            # Test that there are no gaps
            if ( _has_gaps($seq) and ! $allow_gaps ) {
                next SEQ if ( ! $verbose );
                
                # verbose operations
                push @diagnostics, 1;
                $filter_flag = 1;
            }
            else {
                push @diagnostics, 0;
            }
            
            # Test that there are no ambiguous bases
            if ( _has_ambiguous_bases($seq) and ! $allow_ambig_bases ) {
                next SEQ if ( ! $verbose );
                
                # verbose operations
                push @diagnostics, 1;
                $filter_flag = 1;
            }
            else {
                push @diagnostics, 0;
            }
            
            
            ### All the tests are done.  Now do some output.
            
            # if the sequence was flagged as filtered record the diagnostic info
            if ( $filter_flag ) {  
                push @diag_table, \@diagnostics;
            }
            
            # if the sequence is filtered out and we are in verbose mode print
            # the sequence.  Else just print the unfiltered seqs.  Unfiltered
            # means seqs that passed the tests.
            if ( $filter_flag and $verbose ) {
                $filtered_out->write_seq($fastq_seq);
            }
            else {
                $out->write_seq($fastq_seq);
            }
        }
        
        # Some clean up
        if ( $verbose ) {
            _print_diagnostics(\@diag_table, $out_dir, $file_prefix);
        }
        
        return \@diag_table;
    }
    
    sub _init {
        my ($self, $arg_href) = @_;
        
        $self->set_fastq_file($arg_href->{fastq_file});
        $self->set_output_dir($arg_href->{output_dir});
        $self->set_min_len($arg_href->{min_len});
        $self->set_min_avg_qual($arg_href->{min_avg_qual});
        $self->set_min_base_qual($arg_href->{min_base_qual});
        $self->set_allow_gaps($arg_href->{allow_gaps});
        $self->set_allow_ambig_bases($arg_href->{allow_ambig_bases});
        $self->set_verbose($arg_href->{verbose});
        
        return 1;
    }
    
    sub _too_short {
        my ($min_length, $seq) = @_;
        
        if ( ! defined $min_length ) {
            return 0;  # not too short
        }
    
        if ( length $seq >= $min_length ) {
            return 0;  # not too short
        }
        
        return 1;  # is too short
    }
    
    sub _below_avg_qual {
        my ($min_avg_qual, $mean) = @_;
    
        if ( $mean < $min_avg_qual ) {
            return 1;  # is below avg qual
        }
        
        return 0;  # is NOT below avg qual (i.e. okay to keep)
    }
    
    sub _below_min_base_qual {
        my ($min_base_qual, $min) = @_;
    
        if ( $min < $min_base_qual ) {
            return 1;  # is below the min base qual
        }
        
        return 0;
    }
    
    sub _has_gaps {
        my ($seq) = @_;
    
        # Right now the only gap characters I test for are '-' and '.'
        if ( $seq =~ m{-|\.}xms ) {
            return 1;  # has gaps
        }
        
        return 0;  # no gaps
    }
    
    sub _has_ambiguous_bases {
        my ($seq) = @_;
        
        # R    A or G
        # Y    C or T
        # S    G or C
        # W    A or T
        # K    G or T
        # M    A or C
        # B    C or G or T
        # D    A or G or T
        # H    A or C or T
        # V    A or C or G
        # N    any base
        
        if ( $seq =~ m{[RYSWKMBDHVN]}ixms
            ) {
            return 1;
        }
        
        return 0;
    }
    
    sub _get_file_prefix {
        my ($file) = @_;
        
        my $filename = fileparse($file);
        
        my @vals = split /\./, $filename;
        
        if ( scalar @vals == 0 ) {
            cluck( "Unrecognized file format: $file");
            return "seqs";
        }
        elsif ( scalar @vals <= 2 ) {
            return $vals[0];
        }
        else {
            pop @vals;
            my $prefix = join '.', @vals;
            return $prefix;
        }
    }
    
    sub _print_diagnostics {
        my ($aref, $output_dir, $file_prefix) = @_;

        my $file = $output_dir . "/" . $file_prefix . "_diagnostics.txt";
        open my $fh, '>', $file
            or croak("Cannot open file: $file");
        
        print $fh "#seq_id\tlength\taverage_qual\tbase_qual\tgaps\tambiguous_base\n";
        
        foreach my $vals ( @{$aref} ) {
            #my $line = (join "\t", @{$vals}) . "\n";
            print $fh (join "\t", @{$vals}) . "\n";
        }
        
        return 1;
    }
    
    sub set_fastq_file {
        my ($self, $file) = @_;
        
        if ( ! defined $file ) {
            croak("Undefined fastq_file\n$NEW_USAGE");
        }
        if ( ! -s $file ) {
            carp("Empty fastq_file: $file");
        }
        $fastq_file_of{ident $self} = $file;
        
        return 1;
    }
    
    sub set_output_dir {
        my ($self, $output_dir) = @_;

        if ( ! defined $output_dir ) {
            croak("Undefined output dir\n$NEW_USAGE");
        }
        elsif ( ! -d $output_dir ) {
            croak("Output dir doesn't exist: $output_dir");
        }
        $output_dir_of{ident $self} = $output_dir;
        
        return 1;
    }
    
    sub set_min_len {
        my ($self, $min_len) = @_;
        # if min_len is not defined then we want to keep entire seq

        if ( defined $min_len and $min_len < 0 ) {
            croak("min_len must be > 0");
        }
        $min_len_of{ident $self} = $min_len;
        
        return 1;
    }
    
    sub set_min_avg_qual {
        my ($self, $min_avg_qual) = @_;
        
        if ( ! defined $min_avg_qual ) {
            $min_avg_qual = 0;
        }
        if ( $min_avg_qual < 0 ) {
            croak("min_avg_qual must be > 0");
        }
        $min_avg_qual_of{ident $self} = $min_avg_qual;
        
        return 1;
    }
    
    sub set_min_base_qual {
        my ($self, $min_base_qual) = @_;
        
        # DEFAULT = 0
        if ( ! defined $min_base_qual ) {
            $min_base_qual = 0;
        }
        
        if ( $min_base_qual < 0 ) {
            croak("min_base_qual must be > 0");
        }
        
        $min_base_qual_of{ident $self} = $min_base_qual;
        
        return 1;
    }
    
    sub set_allow_gaps {
        my ($self, $allow_gaps) = @_;
        
        # DEFAULT = 0
        if ( ! defined $allow_gaps ) {
            $allow_gaps = 0;  
        }
        
        # translate Yes and No to 0 or 1
        if ( ! looks_like_number($allow_gaps) ) {
            $allow_gaps = _Y_N_translate($allow_gaps);
        }
        
        if ( ! ($allow_gaps == 0 or $allow_gaps == 1) ) {
            croak("allow_gaps must be 0 or 1");
        }
        
        $allow_gaps_of{ident $self} = $allow_gaps;
        
        return 1;
    }
    
    sub set_allow_ambig_bases {
        my ($self, $allow_ambig_bases) = @_;
        
        # DEFAULT = 0 (i.e. no)
        if ( ! defined $allow_ambig_bases ) {
            $allow_ambig_bases = 0;
        }
        
        # translate Yes and No to 0 or 1
        if ( ! looks_like_number($allow_ambig_bases) ) {
            $allow_ambig_bases = _Y_N_translate($allow_ambig_bases);
        }
        
        if ( ! ( $allow_ambig_bases == 0 or $allow_ambig_bases == 1) ) {
            croak("allow_ambig_bases must be 0 or 1");
        }
        
        $allow_ambig_bases_of{ident $self} = $allow_ambig_bases;
        
        return 1;
    }
    
    sub set_verbose {
        my ($self, $verbose) = @_;
        
        # DEFAULT = 0 (i.e. no)
        if ( ! defined $verbose ) {
            $verbose = 0;  # e.g. FALSE
        }
        
        # translate Yes and No to 0 or 1
        if ( ! looks_like_number($verbose) ) {
            $verbose = _Y_N_translate($verbose);
        }
        
        if ( ! ($verbose == 0 or $verbose == 1) ) {
            croak("verbose must be 0 or 1");
        }
        
        $verbose_of{ident $self} = $verbose;
        
        return 1;
    }
    
    sub get_fastq_file {
        my ($self) = @_;
        return $fastq_file_of{ident $self};
    }
    
    sub get_output_dir {
        my ($self) = @_;
        return $output_dir_of{ident $self};
    }
    
    sub get_min_len {
        my ($self) = @_;
        return $min_len_of{ident $self};
    }
    
    sub get_min_avg_qual {
        my ($self) = @_;
        return $min_avg_qual_of{ident $self};
    }
    
    sub get_min_base_qual {
        my ($self) = @_;
        return $min_base_qual_of{ident $self};
    }
    
    sub get_allow_gaps {
        my ($self) = @_;
        return $allow_gaps_of{ident $self};
    }
    
    sub get_allow_ambig_bases {
        my ($self) = @_;
        return $allow_ambig_bases_of{ident $self};
    }
    
    sub get_verbose {
        my ($self) = @_;
        return $verbose_of{ident $self};
    }
    
    sub _Y_N_translate {
        my ($val) = @_;
        
        my %legal_yes = map { $_ => 1 } qw (Y y YES Yes yes );
        
        my %legal_no = map { $_ => 1 } qw (N n NO No no);
        
        if ( defined $legal_yes{$val} ) {
            return 1;
        }
        elsif ( defined $legal_no{$val} ) {
            return 0;
        }
        else {
            croak("Bad Yes/No value: $val");
        }
    }
}


1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::QC::FastqFilter - Filters seqs in a Fastq file based on quality


=head1 VERSION

This document describes BioUtils::QC::FastqFilter version 1.0.1


=head1 SYNOPSIS

    use BioUtils::QC::FastqFilter;

    my $filter = BioUtils::QC::FastqFilter->new({
                    fastq_file => $fastq_file,
                    output_dir => $output_dir,
                    [min_len => $min_len],
                    [min_avg_qual => $min_avg_qual],
                    [min_base_qual => $min_base_qual],
                    [allow_gaps => $allow_gaps],
                    [allow_ambig_bases => $allow_ambig],
                    [verbose => 0],
                });

    $filter->filter();
  
  
=head1 DESCRIPTION

BioUtils::QC::FastqFilter is a module that can be used to filter sequences in a
Fastq file.  These sequences can be filtered based on their length, average
quality score, lowest base quality score, gaps, and ambiguous bases.


=head1 INTERFACE 

    filter

    # Setters #
    set_fastq_file
    set_output_dir
    set_min_len
    set_min_avg_qual
    set_min_base_qual
    set_allow_gaps
    set_allow_ambig_bases
    set_verbose

    # Getters #
    get_fastq_file
    get_output_dir
    get_min_len
    get_min_avg_qual
    get_min_base_qual
    get_allow_gaps
    get_allow_ambig_bases
    get_verbose

    # Private #
    _init
    _too_short
    _below_avg_qual
    _below_base_qual
    _has_gaps
    _has_ambig_bases
    _get_file_prefix
    _print_diagnostics
    _Y_N_translate


=head1 DIAGNOSTICS

=for author to fill in:
    List every single error and warning message that the module can
    generate (even the ones that will "never happen"), with a full
    explanation of each problem, one or more likely causes, and any
    suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back


=head1 CONFIGURATION AND ENVIRONMENT
  
BioUtils::QC::FastqFilter requires the system command 'grep' for testing.  Tests
will fail if 'grep' cannot be found as a system command.


=head1 DEPENDENCIES

    version
    Class::Std::Utils
    Scalar::Util qw(looks_like_number)
    Getopt::Long
    Carp qw(cluck carp croak)
    Pod::Usage
    Statistics::Descriptive
    Data::Dumper qw(Dumper)
    Readonly
    Cwd
    File::Basename
    BioUtils::FastqIO 1.0.1
    BioUtils::Codec::QualityScores qw(illumina_1_8_to_int)
    MyX::Generic


=head1 INCOMPATIBILITIES

    None reported.


=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: my $filter = BioUtils::QC::FastqFilter->new({
                            fastq_file => $fastq_file,
                            output_dir => $output_dir,
                            [min_len => $min_len],
                            [min_avg_qual => $min_avg_qual],
                            [min_base_qual => $min_base_qual],
                            [allow_gaps => $allow_gaps],
                            [allow_ambig_bases => $allow_ambig],
                            [verbose => 0],
                        });
	Function: Creates a new BioUtils::QC::FastqFilter object
	Returns: BioUtils::QC::FastqFilter
	Args: -fastq_file => full file path to fastq file with seqs to filter
          -output_dir => full path to output directory
          -min_len => min length of seqs to keep
          -min_avg_qual => min average quality value of seqs to keep
          -min_base_qual => min quality value for each base of seqs to keep
          -allow_gaps => keep seqs with gaps (0 == NO, 1 == YES)
          -allow_ambig_bases => keep seqs with ambiguous bases
                                (0 == NO, 1 == YES)
	Throws: NA
	Comments: Calls _init()
	See Also: _init()
    
=head2 filter

	Title: filter
	Usage: my $filter->filter();
	Function: Filters seqs in the objects fastq file.  Prints resutls to
              the objects output directory.
	Returns: Array Ref
	Args: NA
	Throws: MyX::Generic::Undef::Param
	Comments: The returned Array ref is mostly used for testing purposes
	See Also: NA
    
=head2 _init

	Title: _init
	Usage: _init($arg_href);
	Function: Sets the objects attribute variables based on the argument href
	Returns: 1 on successful completion
	Args: -arg_href => arguments passed to new
	Throws: NA
	Comments: DEFAULT VALUES:
                -min_len => NA (i.e undef)
                -min_avg_qual => 0
                -min_base_qual => 0
                -allow_gaps => 0 (i.e NO)
                -allow_ambig_bases => 0 (i.e. NO)
                -verbose => 0 (i.e. NO)
	See Also: NA
    
=head2 _too_short

	Title: _too_short
	Usage: _too_short($min_len, $seq);
	Function: Determines if seq is too short
	Returns: Bool
	Args: -min_len => the minimum allowed length
          -seq => the sequence to check
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 _below_avg_qual

	Title: _below_avg_qual
	Usage: _below_avg_qual($min_avg_qual, $mean);
	Function: Determines if seq is above min average quality
	Returns: Bool
	Args: -min_avg_qual => the minimum average quality for all bases
          -mean => the mean quality score for all bases
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _below_min_base_qual

	Title: _below_min_base_qual
	Usage: _below_min_base_qual($min_base_qual, $min);
	Function: Determines if all sequence bases are above min base qual
	Returns: Bool
	Args: -min_base_qual => the minimum allow base quality score
          -min => the lowest quality score from a sequence
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 _has_gaps

	Title: _has_gaps
	Usage: _has_gaps($seq);
	Function: Determines if the sequence has gaps
	Returns: Bool
	Args: -seq => the sequence to check
	Throws: NA
	Comments: Gap characters include '-' and '.'
	See Also: NA
    
=head2 _has_ambiguous_bases

	Title: _has_ambiguous_bases
	Usage: _has_ambiguous_bases($seq);
	Function: Determines if the sequence has ambiguous bases
	Returns: Bool
	Args: -seq => the sequence to check
	Throws: NA
	Comments: Ambiguous bases include R, Y, S, W, K, M, B, D, H, V, and N
	See Also: NA
    
=head2 _get_file_prefix

	Title: _get_file_prefix
	Usage: _get_file_prefix($file);
	Function: Gets the prefix of the file name
	Returns: Str
	Args: -file => the file path
	Throws: NA
	Comments: E.g. /Users/me/my_prefix.txt == my_prefix
	See Also: NA
    
=head2 _print_diagnostics

	Title: _print_diagnostics
	Usage: _print_diagnostics($aref, $output_dir, $file_prefix);
	Function: Print the diagnostics file
	Returns: 1 on successful completion
	Args: -aref => an array ref with the diagnostic information
          -output_dir => the output dir to print the file in
          -file_prefix => the prefix of the output file
	Throws: NA
	Comments: the output file is $output_dir/$prefix_diagnostics.txt
	See Also: NA
    
=head2 set_fastq_file

	Title: set_fastq_file
	Usage: $filter->set_fastq_file($file);
	Function: Sets the input fastq file
	Returns: 1 on successful completion
	Args: -file => a fastq file
	Throws: croak("Undefined fastq_file")
            carp("Empty fastq_file")
	Comments: NA
	See Also: NA
    
=head2 set_output_dir

	Title: set_output_dir
	Usage: $filter->set_output_dir($output_dir);
	Function: Sets the output directory
	Returns: 1 on successful completion
	Args: -dir => a path to where to print the ouput files
	Throws: croak("Undefined output dir")
            carp("Output dir doesn't exist")
	Comments: NA
	See Also: NA
    
=head2 set_min_len

	Title: set_min_len
	Usage: $filter->set_min_len($min_len);
	Function: Sets the min length of sequences allowed
	Returns: 1 on successful completion
	Args: -min_len => int of minimum sequence length
	Throws: croak("min_len must be > 0")
	Comments: NA
	See Also: NA
    
=head2 set_min_avg_qual

	Title: set_min_avg_qual
	Usage: $filter->set_min_avg_qual($min_avg_qual);
	Function: Sets the min average quality across all bases in the seq
	Returns: 1 on successful completion
	Args: -min_avg_qual => int of minimum average quality of all bases in seq
	Throws: croak("min_avg_qual must be > 0")
	Comments: NA
	See Also: NA
    
=head2 set_min_base_qual

	Title: set_min_base_qual
	Usage: $filter->set_min_base_qual($min_base_qual);
	Function: Sets the min quality for the lowest quality base in the seq
	Returns: 1 on successful completion
	Args: -min_base_qual => int of minimum base qual of lowest quality base
	Throws: croak("min_base_qual must be > 0")
	Comments: NA
	See Also: NA
    
=head2 set_allow_gaps

	Title: set_allow_gaps
	Usage: $filter->set_allow_gaps($allow_gaps);
	Function: Sets the boolean for allowing gaps in the sequence
	Returns: 1 on successful completion
	Args: -allow_gaps => 0 == don't allow; 1 == allow
	Throws: croak("allow_gaps must be 0 or 1")
	Comments: You can now also pass Y or N and it will be translated to 0 or 1
	See Also: NA
    
=head2 set_allow_ambig_bases

	Title: set_allow_ambig_bases
	Usage: $filter->set_allow_ambig_bases($allow_ambig_bases);
	Function: Sets the boolean for allowing ambiguous bases in the sequence
	Returns: 1 on successful completion
	Args: -allow_ambig_bases => 0 == don't allow; 1 == allow
	Throws: croak("allow_ambig_bases must be 0 or 1")
	Comments: You can now also pass Y or N and it will be translated to 0 or 1
	See Also: NA
    
=head2 set_verbose

	Title: set_verbose
	Usage: $filter->set_verbose($verbose);
	Function: Sets the boolean for verbose output
	Returns: 1 on successful completion
	Args: -verbose => 0 == normal mode; 1 == verbose mode
	Throws: croak("verbose must be 0 or 1")
	Comments: You can now also pass Y or N and it will be translated to 0 or 1
	See Also: NA
    
=head2 get_fastq_file

	Title: get_fastq_file
	Usage: my $fastq_file = $filter->get_fastq_file();
	Function: Gets the fastq file path
	Returns: String
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_output_dir

	Title: get_output_dir
	Usage: my $output_dir = $filter->get_output_dir();
	Function: Gets the fastq file path
	Returns: String
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_min_len

	Title: get_min_len
	Usage: my $min_len = $filter->get_min_len();
	Function: Gets the min_len attribute value
	Returns: Int
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_min_avg_qual

	Title: get_min_len
	Usage: my $min_avg_qual = $filter->get_min_avg_qual();
	Function: Gets the min_avg_qual attribute value
	Returns: Int
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_min_base_qual

	Title: get_min_base_qual
	Usage: my $min_base_qual = $filter->get_min_base_qual();
	Function: Gets the min_base_qual attribute value
	Returns: Int
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_allow_gaps

	Title: get_allow_gaps
	Usage: my $allow_gaps = $filter->get_allow_gaps();
	Function: Gets the allow_gaps attribute value
	Returns: Bool (0 or 1)
	Args: NA
	Throws: NA
	Comments: 0 == don't allow gaps; 1 == allow gaps
	See Also: NA
    
=head2 get_allow_ambig_bases

	Title: get_allow_ambig_bases
	Usage: my $allow_ambig_bases = $filter->get_allow_ambig_bases();
	Function: Gets the allow_ambig_bases attribute value
	Returns: Bool (0 or 1)
	Args: NA
	Throws: NA
	Comments: 0 == don't allow ambig bases; 1 == allow ambig bases
	See Also: NA
    
=head2 get_verbose

	Title: get_verbose
	Usage: my $verbose = $filter->get_verbose();
	Function: Gets the verbose attribute value
	Returns: Bool (0 or 1)
	Args: NA
	Throws: NA
	Comments: 0 == normal mode; 1 == verbose mode
	See Also: NA
    
=head2 _Y_N_translate

	Title: _Y_N_translate
	Usage: _Y_N_translate($val);
	Function: Translates Yes or No values into 0 or 1
	Returns: Bool (0 or 1)
	Args: -val => a Yes or No value
	Throws: croak("Bad Yes/No value:")
	Comments: 0 == No; 1 == Yes
              Valid Yes values include: Y, y, YES, Yes, yes
              Valid No values include: N, n, NO, No, no
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to Scott Yourstone
C<< <scott.yourstone81@gmail.com> >>

=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone C<< <scott.yourstone81@gmail.com> >>. All rights reserved.

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
