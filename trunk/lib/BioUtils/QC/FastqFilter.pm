package BioUtils::QC::FastqFilter;

use warnings;
use strict;

use version; our $VERSION = qv('1.2.1');

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
use BioUtils::FastqIO 1.2.1;
use BioUtils::Codec::QualityScores qw(illumina_1_8_to_int);
use MyX::Generic 0.0.1;


{
    Readonly my $NEW_USAGE => q{ new( {fastq_file => ,
            output_dir => ,
            [fastq_file2 => ,]
            [trim_to_base => ,]
            [min_len => ],
            [min_avg_qual => ],
            [min_base_qual => ],
            [min_c_score => ],
            [allow_gaps => ],
            [allow_ambig_bases => ],
            [verbose => ],
            } ) };
    
    
    # Attributes #
    my %fastq_file_of;
    my %fastq_file2_of;
    my %output_dir_of;
    my %trim_to_base_of;
    my %min_len_of;
    my %min_avg_qual_of;
    my %min_base_qual_of;
    my %min_c_score_of;
    my %allow_gaps_of;
    my %allow_ambig_bases_of;
    my %verbose_of;

    # Setters #
    sub set_fastq_file;
    sub set_fastq_file2;
    sub set_output_dir;
    sub set_trim_to_base;
    sub set_min_len;
    sub set_min_avg_qual;
    sub set_min_base_qual;
    sub set_min_c_score;
    sub set_allow_gaps;
    sub set_allow_ambig_bases;
    sub set_verbose;
    
    # Getters #
    sub get_fastq_file;
    sub get_fastq_file2;
    sub get_output_dir;
    sub get_trim_to_base;
    sub get_min_len;
    sub get_min_avg_qual;
    sub get_min_base_qual;
    sub get_min_c_score;
    sub get_allow_gaps;
    sub get_allow_ambig_bases;
    sub get_verbose;
    
    # Others #
    sub filter;
    sub filter_pairs;
    sub _make_io_objs;
    sub _make_pairs_io_objs;
    sub _trim;
    sub _test_seq;
    sub _init;
    sub _too_short;
    sub _below_avg_qual;
    sub _below_base_qual;
    sub _below_min_c_score;
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
        my $verbose = $self->get_verbose();
        my @diag_table = (); # a table for reasons for each sequences filtering
        my $file_prefix = _get_file_prefix($input_file); # for making output files
    
        # Build the FASTQIO objects
        my ($in, $out, $LQ) = _make_io_objs($input_file,
                                            $out_dir,
                                            $file_prefix,
                                            $verbose); 
        
        # Read in the sequences from $input_file
        # AND run each filter test
        SEQ: while ( my $fastq_seq = $in->get_next_seq() ) {
            
            $fastq_seq = $self->_trim($fastq_seq);
            
            my @diagnostics = ();
            my $seq_id = $fastq_seq->get_id();
            my $seq_header = $fastq_seq->get_header();
            my $seq_str = $fastq_seq->get_seq();
            my $quals_str = $fastq_seq->get_quals_str();
            
            # run tests
            my $LQ_flag = $self->_test_seq(\@diagnostics,
                                            $seq_id,
                                            $seq_header,
                                            $seq_str,
                                            $quals_str,
                                            $verbose);
            
            
            ### All the tests are done.  Now do some output.
            
            # if the sequence was flagged as filtered (ie LQ--low qual) record
            # the diagnostic info
            if ( $LQ_flag ) {  
                push @diag_table, \@diagnostics;
            }
            
            # if the sequence is filtered out and we are in verbose mode print
            # the sequence.  Else just print the HQ seqs.  HQ
            # means seqs that passed the tests.
            if ( $LQ_flag and $verbose ) {
                $LQ->write_seq($fastq_seq);
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
    
    sub filter_pairs {
        my ($self) = @_;
        
        my $out_dir = $self->get_output_dir();
        my $fwd_file = $self->get_fastq_file();
        my $rev_file = $self->get_fastq_file2();
        my $verbose = $self->get_verbose();
        my @diag_table = (); # a table for reasons for each sequences filtering
        my $fwd_prefix = _get_file_prefix($fwd_file); # for making output files
        my $rev_prefix = _get_file_prefix($rev_file); # for making output files
        
        # Build the FASTQIO objects
        my ($fwd_in, $fwd_out, $fwd_LQ,
            $rev_in, $rev_out, $rev_LQ) = _make_pairs_io_objs(
                                                        $fwd_file,
                                                        $rev_file,
                                                        $out_dir,
                                                        $fwd_prefix,
                                                        $rev_prefix,
                                                        $verbose);
        
        # Read in the sequences and run each filter test
        SEQ: while ( my $fwd_seq = $fwd_in->get_next_seq() ) {
            # check that the rev seq is defined
            my $rev_seq = $rev_in->get_next_seq();
            if ( $rev_seq == 0) {
                # when 0 is returned by FastqIO get_next_seq() was called when
                # no more sequences remain.
                croak ("Mismatched Pairs: fewer rev seqs");
            }
            
            # check that the fwd and rev ids match
            # I comment this out because there are times when the names do not
            #   match and I dno't want to overload the output with this warning
            #if ( $fwd_seq->get_id() ne $rev_seq->get_id() ) {
            #    carp("Fwd and Rev IDs do NOT match:\n" .
            #         "\tFWD: " . $fwd_seq->get_id() . "\n" .
            #         "\tREV: " . $rev_seq->get_id() . "\n");
            #}
            
            # trim
            $fwd_seq = $self->_trim($fwd_seq);
            $rev_seq = $self->_trim($rev_seq);
            
            my @diagnostics = ();
            my $seq_id = $fwd_seq->get_id() . "-" . $rev_seq->get_id();
            my $seq_header = $fwd_seq->get_header() . " | " . $rev_seq->get_header();
            my $seq_str = $fwd_seq->get_seq() . $rev_seq->get_seq();
            my $quals_str = $fwd_seq->get_quals_str() .
                            $rev_seq->get_quals_str();
            
            ### run tests ###
            my $LQ_flag = $self->_test_seq(\@diagnostics,
                                               $seq_id,
                                               $seq_header,
                                               $seq_str,
                                               $quals_str,
                                               $verbose);
            
            ### Now output ###
            # if the sequence was flagged as filtered record the diagnostic info
            if ( $LQ_flag) {
                push @diag_table, \@diagnostics;
            }
            
            # if the sequence is filtered out and we are in verbose mode print
            # the sequence.  Else just print the HQ seqs.  HQ
            # means seqs that passed the tests.
            if ( $LQ_flag and $verbose ) {
                $fwd_LQ->write_seq($fwd_seq);
                $rev_LQ->write_seq($rev_seq);
            }
            else {
                $fwd_out->write_seq($fwd_seq);
                $rev_out->write_seq($rev_seq);
            }
            
        }
        
        # check to see that all the rev seqs were used
        if ( $rev_in->get_next_seq() != 0 ) {
            # when 0 is returned by FastqIO get_next_seq() was called when
            # no more sequences remain.
            croak("Mismatched Pairs: fewer fwd seqs");
        }
        
        ### Some cleanup ###
        if ( $verbose ) {
            # NOTE: the diagnostics output file is made based on the fwd
            # file prefix
            _print_diagnostics(\@diag_table, $out_dir, $fwd_prefix);
        }
        
        return \@diag_table;
    }
    
    sub _make_io_objs {
        my ($input_file, $out_dir, $file_prefix, $verbose) = @_;
        
        my ($in, $out, $LQ);
        
        eval {
            $in = BioUtils::FastqIO->new( {
                        stream_type => '<',
                        file => $input_file
                    });
            $out = BioUtils::FastqIO->new( {
                        stream_type => '>',
                        file => $out_dir . "/" .
                                $file_prefix .
                                '_HQ.fastq'
                    });
            if ( $verbose ) {
                $LQ = BioUtils::FastqIO->new( {
                                    stream_type => '>',
                                    file => $out_dir . "/" .
                                            $file_prefix .
                                            '_LQ.fastq'
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
        
        return ($in, $out, $LQ);
    }
    
    sub _make_pairs_io_objs {
        my ($fwd_file, $rev_file, $out_dir,
            $fwd_prefix, $rev_prefix, $verbose) = @_;
        
        my ($fwd_in, $fwd_out, $fwd_LQ,
            $rev_in, $rev_out, $rev_LQ);
        
        eval {
            $fwd_in = BioUtils::FastqIO->new( {
                        stream_type => '<',
                        file => $fwd_file
                    });
            $fwd_out = BioUtils::FastqIO->new( {
                        stream_type => '>',
                        file => $out_dir . "/" .
                                $fwd_prefix .
                                '_HQ.fastq'
                    });
            if ( $verbose ) {
                $fwd_LQ = BioUtils::FastqIO->new( {
                                    stream_type => '>',
                                    file => $out_dir . "/" .
                                            $fwd_prefix .
                                            '_LQ.fastq'
                                });
            }
            $rev_in = BioUtils::FastqIO->new( {
                        stream_type => '<',
                        file => $rev_file
                    });
            $rev_out = BioUtils::FastqIO->new( {
                        stream_type => '>',
                        file => $out_dir . "/" .
                                $rev_prefix .
                                '_HQ.fastq'
                    });
            if ( $verbose ) {
                $rev_LQ = BioUtils::FastqIO->new( {
                                    stream_type => '>',
                                    file => $out_dir . "/" .
                                            $rev_prefix .
                                            '_LQ.fastq'
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
        
        return ($fwd_in, $fwd_out, $fwd_LQ,
                $rev_in, $rev_out, $rev_LQ);
    }
    
    sub _trim {
        my ($self, $fastq_seq) = @_;
        
        if ( ! defined $self->get_trim_to_base() ) { return $fastq_seq; }
        if ( $self->get_trim_to_base() =~ m/na/i ) { return $fastq_seq; }
        
        # The trim_front subroutine returns a FastqSeq object of the portion
        #   that was trimmed off the front.  So I resent $fastq_seq to that
        #   portion.  This is trimming off the first $self->get_trim_to_base()
        #   bases and returning that as the new fastq_seq.  See doc for
        #   BioUtils::FastqSeq for more info.

        return $fastq_seq->trim_front($self->get_trim_to_base());
    }
    
    sub _test_seq {
        my ($self, $diag_aref, $id, $header, $seq_str, $quals_str, $verbose) = @_;
        
        # A boolean to keep track if the sequence needs to be filtered out
        # LQ stands for low quality
        my $LQ_flag = 0;

        # A data structure for storing this (a single) sequence's diagnostics
        push @$diag_aref, $id;
        
        # Decode the quals into ints and store in an aref
        my $int_quals_aref = illumina_1_8_to_int($quals_str);
        #print Dumper($int_quals_aref) . "\n";
        
        # Create a statistics descriptive object for qual score testing
        my $stat = Statistics::Descriptive::Sparse->new();
        $stat->add_data(@{$int_quals_aref});
        
        # Test length
        if ( _too_short($self->get_min_len(), $seq_str) ) {
            return $LQ_flag if ( ! $verbose );
            
            # verbose operations
            push @$diag_aref, 1;
            $LQ_flag = 1;
        }
        else {
            push @$diag_aref, 0;
        }
        
        # Test minimum average quality over whole read
        if ( _below_avg_qual($self->get_min_avg_qual(), $stat->mean()) ) {
            return 1 if ( ! $verbose );
            
            # verbose operations
            push @$diag_aref, 1;
            $LQ_flag = 1;
        }
        else {
            push @$diag_aref, 0;
        }
        
        # Test minimum base quality
        if ( _below_min_base_qual(
                                  $self->get_min_base_qual(),
                                  $stat->min())
            ) {
            return 1 if ( ! $verbose );
            
            # verbose operations
            push @$diag_aref, 1;
            $LQ_flag = 1;
        }
        else {
            push @$diag_aref, 0;
        }
        
        # Test min c-score
        if ( _below_min_c_score($self->get_min_c_score(), $header) ) {
            return 1 if ( ! $verbose);
            
            # verbose operations
            push @$diag_aref, 1;
            $LQ_flag = 1;
        }
        else {
            push @$diag_aref, 0;
        }
        
        # Test that there are no gaps
        if ( _has_gaps($seq_str) and ! $self->get_allow_gaps() ) {
            return 1 if ( ! $verbose );
            
            # verbose operations
            push @$diag_aref, 1;
            $LQ_flag = 1;
        }
        else {
            push @$diag_aref, 0;
        }
        
        # Test that there are no ambiguous bases
        if ( _has_ambiguous_bases($seq_str) and
            ! $self->get_allow_ambig_bases ) {
            return 1 if ( ! $verbose );
            
            # verbose operations
            push @$diag_aref, 1;
            $LQ_flag = 1;
        }
        else {
            push @$diag_aref, 0;
        }
        
        return $LQ_flag;
    }
    
    sub _init {
        my ($self, $arg_href) = @_;
        
        $self->set_fastq_file($arg_href->{fastq_file});
        $self->set_fastq_file2($arg_href->{fastq_file2});
        $self->set_output_dir($arg_href->{output_dir});
        $self->set_trim_to_base($arg_href->{trim_to_base});
        $self->set_min_len($arg_href->{min_len});
        $self->set_min_avg_qual($arg_href->{min_avg_qual});
        $self->set_min_base_qual($arg_href->{min_base_qual});
        $self->set_min_c_score($arg_href->{min_c_score});
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
    
    sub _below_min_c_score {
        my ($min_c_score, $header) = @_;
        
        if ( $header =~ m/\|/ ) {
            # this means it is a pair
            my @parts = split /\|/, $header;
            foreach (@parts) {
                if ( $_ =~ m/c_score\s*=\s*(\d+\.*\d*)/ ) {
                    if ( $1 < $min_c_score ) {
                        return 1;  # below the min c_score
                    }
                }
            }
        }
        else {
            if ( $header =~ m/c_score\s*=\s*(\d+\.*\d*)/ ) {
                if ( $1 < $min_c_score ) {
                    return 1;  # below the min c_score
                }
            }
        }
        
        return 0;  # c_score is either okay or not present
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
        
        print $fh "#seq_id\tlength\taverage_qual\tbase_qual\tc_score\tgaps\tambiguous_base\n";
        
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
    
    sub set_fastq_file2 {
        my ($self, $file) = @_;
        
        if ( defined $file and ! -s $file ) {
            carp("Empty fastq_file2: $file");
        }
        $fastq_file2_of{ident $self} = $file;
        
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
    
    sub set_trim_to_base {
        my ($self, $last_base) = @_;

        if ( defined $last_base and $last_base =~ m/na/i ) {
            $last_base = undef;
        }
        
        if ( defined $last_base and $last_base < 0 ) {
            croak("trim_to_base must be > 0");
        }
        
        $trim_to_base_of{ident $self} = $last_base;
        
        return 1;
    }
    
    sub set_min_len {
        my ($self, $min_len) = @_;
        # if min_len is not defined then we want to keep entire seq
        
        if ( defined $min_len and $min_len =~ m/na/i ) { $min_len = undef; }

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
    
    sub set_min_c_score {
        my ($self, $min_c_score) = @_;
        
        # DEFAULT = 0
        if ( ! defined $min_c_score ) {
            $min_c_score = 0;
        }
        
        if ( $min_c_score < 0 ) {
            croak("min_c_score must be > 0");
        }
        
        $min_c_score_of{ident $self} = $min_c_score;
        
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
    
    sub get_fastq_file2 {
        my ($self) = @_;
        return $fastq_file2_of{ident $self};
    }
    
    sub get_output_dir {
        my ($self) = @_;
        return $output_dir_of{ident $self};
    }
    
    sub get_trim_to_base {
        my ($self) = @_;
        return $trim_to_base_of{ident $self};
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
    
    sub get_min_c_score {
        my ($self) = @_;
        return $min_c_score_of{ident $self};
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

This document describes BioUtils::QC::FastqFilter version 1.2.1


=head1 SYNOPSIS

    use BioUtils::QC::FastqFilter;

    my $filter = BioUtils::QC::FastqFilter->new({
                    fastq_file => $fastq_file,
                    output_dir => $output_dir,
                    [trim_to_base => $last_base],
                    [fastq_file2 => $rev_file],
                    [min_len => $min_len],
                    [min_avg_qual => $min_avg_qual],
                    [min_base_qual => $min_base_qual],
                    [min_c_score => $min_c_score],
                    [allow_gaps => $allow_gaps],
                    [allow_ambig_bases => $allow_ambig],
                    [verbose => 0],
                });

    $filter->filter();
  
  
=head1 DESCRIPTION

BioUtils::QC::FastqFilter is a module that can be used to filter sequences in a
Fastq file.  These sequences can be filtered based on their length, average
quality score, lowest base quality score, gaps, and ambiguous bases.

Sequences can also be filtered by pairs.  For example, if two files are provided
via the attributes fastq_file and fastq_file2, FastqFilter will assume these
files are matching pairs.  For implementation simplicity sequences and quality
values from these pairs are concatenated together.  Test are done on these
concatenated sequences and are NOT done in the individual FWD and REV reads.
For example, if the user provides a pair of fastq files and wants to filter out
any pair that has more than 150bp, the pairs are concatenated and subsiquently
checked to ensure the concatenated sequence is less than 150bp.  Currently it is
not possible to apply filters baesed on individual FWD or REV reads.


=head1 INTERFACE 

    filter
    filter_pairs

    # Setters #
    set_fastq_file
    set_fastq_file2
    set_output_dir
    set_trim_to_base
    set_min_len
    set_min_avg_qual
    set_min_base_qual
    set_min_c_score
    set_allow_gaps
    set_allow_ambig_bases
    set_verbose

    # Getters #
    get_fastq_file
    get_fastq_file2
    get_output_dir
    get_trim_to_base
    get_min_len
    get_min_avg_qual
    get_min_base_qual
    get_min_c_score
    get_allow_gaps
    get_allow_ambig_bases
    get_verbose

    # Private #
    _make_io_objs
    _make_pairs_io_objs
    _trim
    _test_seq
    _init
    _too_short
    _below_avg_qual
    _below_base_qual
    _below_min_c_score
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

=item C<< MyX::Generic::Undef::Param >>

Thrown when trying to access undefined paramters

=item C<< Mismatch Pairs: fewer rev seqs >>

Croaks when filter_pairs is called and there are more FWD seqs than REV seqs

=item C<< Mismatch Pairs: fewer fwd seqs >>

Croaks when filter_pairs is called and there are more REV seqs than FWD seqs

=item C<< Fwd and Rev IDs do NOT match >>

Carps when two FWD and REV read IDs do not match.  You can ignore these warnings
if you know the reads match

=item C<< Undefined fastq_file >>

Croaks when setting an undefined fastq file

=item C<< Empty fastq_file >>

Carps when setting an empty fastq file

=item C<< Undefined fastq_file2 >>

Croaks when setting an undefined fastq file for the matching pairs file

=item C<< Empty fastq_file2 >>

Carps when setting an empty fastq file for the matching pairs file

=item C<< Undefined output dir >>

Croaks when setting an undefined output directory

=item C<< Outut dir doesn't exist >>

Croaks when the output dir doesn't exist

=item C<< min_len must be > 0 >>

Croaks when setting min_len to less than 0

=item C<< min_avg_qual must be > 0 >>

Croaks when setting min_avg_qual to less than 0

=item C<< min_base_qual must be > 0 >>

Croaks when setting min_base_qual to less than 0

=item C<< min_c_score must be > 0 >>

Croaks when setting min_c_score to less than 0

=item C<< allow_gaps must be 0 or 1 >>

Croaks when setting allow_gaps to something other than 0 or 1

=item C<< allow_ambig_bases must be 0 or 1 >>

Croaks when setting allow_ambig_bases to something other than 0 or 1

=item C<< verbose must be 0 or 1 >>

Croaks when setting verbose to something other than 0 or 1

=item C<< Bad Yes/No value >>

Croaks when transforming a Yes or No parameter into a boolean value

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
    BioUtils::FastqIO 1.2.1
    BioUtils::Codec::QualityScores qw(illumina_1_8_to_int)
    MyX::Generic


=head1 INCOMPATIBILITIES

    None reported.


=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: my $filter = BioUtils::QC::FastqFilter->new({
                            fastq_file => $fastq_file,
                            fastq_file2 => $fastq_file2,
                            output_dir => $output_dir,
                            [trim_to_base => $last_base],
                            [min_len => $min_len],
                            [min_avg_qual => $min_avg_qual],
                            [min_base_qual => $min_base_qual],
                            [min_c_score => $min_c_score],
                            [allow_gaps => $allow_gaps],
                            [allow_ambig_bases => $allow_ambig],
                            [verbose => 0],
                        });
	Function: Creates a new BioUtils::QC::FastqFilter object
	Returns: BioUtils::QC::FastqFilter
	Args: -fastq_file => full file path to fastq file with seqs to filter
          -fastq_file2 => full file path to fastq file with matching pairs
          -output_dir => full path to output directory
          -min_len => min length of seqs to keep
          -min_avg_qual => min average quality value of seqs to keep
          -min_base_qual => min quality value for each base of seqs to keep
          -min_c_score => min consensus score.  See below for details
          -allow_gaps => keep seqs with gaps (0 == NO, 1 == YES)
          -allow_ambig_bases => keep seqs with ambiguous bases
                                (0 == NO, 1 == YES)
	Throws: See _init() and setter methods
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
    
=head2 filter_pairs

	Title: filter_pairs
	Usage: my $filter->filter_pairs();
	Function: Filters seqs in the objects fastq file.  Prints resutls to
              the objects output directory.  Filtering is done based on
              concatenating the pairs sequences and quality scores.  However,
              seperate files are output for FWD and REV reads.
	Returns: Array Ref
	Args: NA
	Throws: MyX::Generic::Undef::Param
            craok(Mismatched Pairs: fewer rev seqs)
            croak(Mismatched Pairs: fewer fwd seqs)
            carp(Fwd and Rev IDs do NOT match)
	Comments: The returned Array ref is mostly used for testing purposes
	See Also: NA
    
=head2 _make_io_objs

	Title: _make_io_objs
	Usage: my ($in, $out, $LQ) = _make_io_objs(
                                                $input_file,
                                                $out_dir,
                                                $file_prfix,
                                                $verbose
                                              );
	Function: Instantiates the FastqIO objects used in filter
	Returns: Array of three FastqIO objects
	Args: input_file => path to the input fastq file
          out_dir => path to the output directory
          file_prefix => prefix of the input fastq file
          verbose => verbose attribute
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 _make_pairs_io_objs

	Title: _make_pairs_io_objs
	Usage: my ($fwd_in, $fwd_out, $fwd_LQ,
               $rev_in, $rev_out, $rev_LQ) =
                    $self->_make_pairs_io_objs(
                                                    $fwd_file,
                                                    $rev_file,
                                                    $out_dir,
                                                    $fwd_prefix,
                                                    $rev_prefix,
                                                    $verbose
                                                );
	Function: Instantiates the FastqIO objects used in filter
	Returns: Array of six FastqIO objects
	Args: fwd_file => path to the input fwd fastq file
          rev_file => path to the input rev fastq file
          out_dir => path to the output directory
          fwd_prefix => prefix of the input fwd fastq file
          rev_prefix => prefix of the input rev fastq file
          verbose => verbose attribute
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 _trim

	Title: _trim
	Usage: my $fastq_seq_obj = $self->_trim($fastq_seq_obj);
	Function: Trims the fastq_seq_obj to the base stored in trim_to_base
              attribute.
	Returns: BioUtils::FastqSeq object where the sequence and qual values are of
             length trim_to_base, the attribute stored in this object
	Args: fastq_seq_obj => BioUtils::Fastq object
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 _test_seq

	Title: _test_seq
	Usage: my $LQ_flag = $self->_test-seq(
                                        $diag_aref,
                                        $id,
                                        $seq_str,
                                        $quals_str,
                                        $verbose);
	Function: Tests sequences to see if they should be filtered out and why
	Returns: Boolean informing wither to filter out (1 == filter out)
	Args: diag_aref => an array ref with diagnostic information about why a seq
                       is filtered
          id => sequence id
          seq_str => sequence string
          quals_str => quality values string
          verbose => verbose attribute
	Throws: NA
	Comments: Currently the quality values must be in the illumina 1.8 format.
              When _test_seq is run on pairs the sequence and quality string
              should be concatenated before running _test_seq.
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
    
=head2 _below_min_c_score

	Title: _below_min_c_score
	Usage: _below_min_c_score($min_c_score, $header);
	Function: Determines if the sequence is above the minimum c_score in header
	Returns: Bool
	Args: -min_c_score => the minimum allowed c_score
          -header => the header which may contain the c_score value
	Throws: NA
	Comments: A c_score is a metric for grading a consensus sequnece.  It says
              how similar the sequences that were used to build the consensus
              are to each other.  It can have values between 0 and 40 because it
              is based on quality scores.  If a fastq sequence has a c_score it
              should be declared in the header using this format: c_score=40.0.
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

=head2 set_fastq_file2

	Title: set_fastq_file2
	Usage: $filter->set_fastq_file2($rev_file);
	Function: Sets the input fastq file of reverse matched pairs
	Returns: 1 on successful completion
	Args: -rev_file => a fastq file with matching pairs to the first fastq_file
	Throws: croak("Undefined fastq_file2")
            carp("Empty fastq_file2")
	Comments: NA
	See Also: NA
    
=head2 set_output_dir

	Title: set_output_dir
	Usage: $filter->set_output_dir($output_dir);
	Function: Sets the output directory
	Returns: 1 on successful completion
	Args: -dir => a path to where to print the ouput files
	Throws: croak("Undefined output dir")
            croak("Output dir doesn't exist")
	Comments: NA
	See Also: NA
    
=head2 set_trim_to_base

	Title: set_trim_to_base
	Usage: $filter->set_trim_to_base($last_base);
	Function: Before running filtering protcol trim the sequence to the given
              index.
	Returns: 1 on successful completion
	Args: -last_base => int representing the index of the last base to keep
	Throws: croak("min_len must be > 0")
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
    
=head2 set_min_c_score

	Title: set_min_c_score
	Usage: $filter->set_min_c_score($min_c_score);
	Function: Sets the min c_score
	Returns: 1 on successful completion
	Args: -min_c_score => int of minimum c_score
	Throws: croak("min_c_score must be > 0")
	Comments: NA
	See Also: _below_min_c_score
    
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
    
=head2 get_fastq_file2

	Title: get_fastq_file2
	Usage: my $fastq_file2 = $filter->get_fastq_file2();
	Function: Gets the fastq file path for the matched pairs files
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

=head2 get_trim_to_base

	Title: get_trim_to_base
	Usage: my $last_base = $filter->get_trim_to_base();
	Function: Gets the trim_to_base attribute value
	Returns: Int
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
    
=head2 get_min_c_score

	Title: get_min_c_score
	Usage: my $min_c_score = $filter->get_min_c_score();
	Function: Gets the min_c_score attribute value
	Returns: Digit
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

Copyright (c) 2013, Scott Yourstone
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


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
