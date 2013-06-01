package BioUtils::QC::ContaminantFilter;

use warnings;
use strict;
use Carp qw(carp croak);
use Class::Std::Utils;
use Readonly;
use List::MoreUtils qw(any);
use MyX::Generic 1.0.2;
use BioUtils::FastaSeq 1.0.2;
use BioUtils::FastaIO 1.0.2;
use File::Temp qw(tempfile);
use version; our $VERSION = qv('1.0.2');

{
    Readonly my $NEW_USAGE => q{ new( {params_file => } ) };
    
    # Attributes #
    my %blast_db_of;
    my %query_file_of;
    my %output_dir_of;
    my %output_prefix_of;
    my %eval_of;
    my %perc_iden_of;
    my %output_fmt_of;
    my %max_targets_of;
    my %otu_table_of;
    
    # Setters #
    sub set_blast_db;
    sub set_query_file;
    sub set_output_dir;
    sub set_output_prefix;
    sub set_eval;
    sub set_perc_iden;
    sub set_output_fmt;
    sub set_max_targets;
    sub set_otu_table;
    
    # Getters #
    sub get_blast_db;
    sub get_query_file;
    sub get_output_dir;
    sub get_output_prefix;
    sub get_eval;
    sub get_perc_iden;
    sub get_output_fmt;
    sub get_max_targets;
    sub get_otu_table;
    
    # Others #
    sub run_filter;
    sub _init;
    sub _run_blast;
    sub _parse_blast_file;
    sub _print_results;
    sub _sequence_printing;
    sub _otu_table_printing;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Initialize the optional objects attributes        
        $new_obj->_init($arg_href);
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_blast_db {
        my ($self, $db) = @_;
        
        if ( ! -s "$db.nhr" ) {
            system("makeblastdb -in $db -dbtype nucl");
        }
        
        $blast_db_of{ident $self} = $db;
        return 1;
    }
    
    sub set_query_file {
        my ($self, $file) = @_;
        
        $query_file_of{ident $self} = $file;
        return 1;
    }
    
    sub set_output_dir {
        my ($self, $dir) = @_;
        
        $output_dir_of{ident $self} = $dir;
        return 1;
    }
    
    sub set_output_prefix {
        my ($self, $prefix) = @_;
        
        $output_prefix_of{ident $self} = $prefix;
        return 1;
    }
    
    sub set_eval {
        my ($self, $e) = @_;
        
        $eval_of{ident $self} = $e;
        return 1;
    }
    
    sub set_perc_iden {
        my ($self, $p) = @_;
        
        $perc_iden_of{ident $self} = $p;
        return 1;
    }
    
    sub set_output_fmt {
        my ($self, $fmt) = @_;
        
        $output_fmt_of{ident $self} = $fmt;
        return 1;
    }
    
    sub set_max_targets {
        my ($self, $t) = @_;
        
        $max_targets_of{ident $self} = $t;
        return 1;
    }
    
    sub set_otu_table {
        my ($self, $t) = @_;
        
        $otu_table_of{ident $self} = $t;
        return 1;
    }
    

    ###########
    # Getters #
    ###########
    sub get_blast_db {
        my ($self) = @_;
        return $blast_db_of{ident $self};
    }
    
    sub get_query_file {
        my ($self) = @_;
        return $query_file_of{ident $self};
    }
    
    sub get_output_dir {
        my ($self) = @_;
        return $output_dir_of{ident $self};
    }
    
    sub get_output_prefix {
        my ($self) = @_;
        return $output_prefix_of{ident $self};
    }
    
    sub get_eval {
        my ($self) = @_;
        return $eval_of{ident $self};
    }
    
    sub get_perc_iden {
        my ($self) = @_;
        return $perc_iden_of{ident $self};
    }
    
    sub get_output_fmt {
        my ($self) = @_;
        return $output_fmt_of{ident $self};
    }
    
    sub get_max_targets {
        my ($self) = @_;
        return $max_targets_of{ident $self};
    }
    
    sub get_otu_table {
        my ($self) = @_;
        return $otu_table_of{ident $self};
    }
    
    
    ##########
    # Others #
    ##########
    sub run_filter {
        my ($self) = @_;
        
        my $blast_output_file = $self->_run_blast();
        my $contaminant_names_href = $self->_parse_blast_file($blast_output_file);
        $self->_print_results($contaminant_names_href);
        
        return 1;
    }
    
    sub _init {
        my ($self, $arg_href) = @_;
        
        if ( defined $arg_href->{blast_db} ) {
            $self->set_blast_db($arg_href->{blast_db});
        }
        if ( defined $arg_href->{query_file} ) {
            $self->set_query_file($arg_href->{query_file});
        }
        if ( defined $arg_href->{output_dir} ) {
            $self->set_output_dir($arg_href->{output_dir});
        }
        if ( defined $arg_href->{output_prefix} ) {
            $self->set_output_prefix($arg_href->{output_prefix});
        }
        if ( defined $arg_href->{eval} ) {
            $self->set_eval($arg_href->{eval});
        }
        if ( defined $arg_href->{perc_iden} ) {
            $self->set_perc_iden($arg_href->{perc_iden});
        }
        if ( defined $arg_href->{output_fmt} ) {
            $self->set_output_fmt($arg_href->{output_fmt});
        }
        else {
            $self->set_output_fmt(6);
        }
        if ( defined $arg_href->{max_targets} ) {
            $self->set_max_targets($arg_href->{max_targets});
        }
        else {
            $self->set_max_targets(1);
        }
        if ( defined $arg_href->{otu_table} ) {
            $self->set_otu_table($arg_href->{otu_table});
        }
        
        return 1;
    }
    
    sub _run_blast {
        my ($self) = @_;
        
        my ($fh, $output_file) = tempfile();
        close($fh);
        
        my $database = $self->get_blast_db();
        my $seqs_file = $self->get_query_file();
        my $evalue = $self->get_eval();
        my $perc_identity = $self->get_perc_iden();
        my $output_fmt = $self->get_output_fmt();
        my $max_targets = $self->get_max_targets();

        
        my $command = "blastn " .
                      "-db $database " .
                      "-query $seqs_file " .
                      "-out $output_file ";
        
        # add in optional paramters
        if ( defined $evalue ) {
            $command .= "-evalue $evalue ";
        }
        if ( defined $perc_identity ) {
            $command .= "-perc_identity $perc_identity ";
        }
        if ( defined $output_fmt ) {
            $command .= "-outfmt $output_fmt ";
        }
        if ( defined $max_targets ) {
            $command .= "-max_target_seqs $max_targets";
        }
        
        print "BLAST Command:\n $command\n";
        system("$command");
        
        return $output_file;
    }
    
    sub _parse_blast_file {
        my ($self, $file) = @_;
        
        # This method will only work when out_fmt is 6
        
        my %contaminant_names = ();  # NAME => contaminantX
        
        open my $BLS, "<", "$file" or
            croak("Cannot open temp BLAST output file: $file");
        
        foreach my $line ( <$BLS> ) {
            chomp $line;
            my @values = split /\t/, $line;
            $contaminant_names{$values[0]} = $values[1];
        }
        
        return \%contaminant_names;
    }
    
    sub _print_results {
        my ($self, $contaminant_names_href) = @_;
        
        $self->_sequence_printing($contaminant_names_href);
        
        if ( defined $self->get_otu_table() ) {
            $self->_otu_table_printing($contaminant_names_href);
        }
        
        return 1;
    }

    sub _sequence_printing {
        my ($self, $contam_names_href) = @_;
        
        my $seq_file = $self->get_query_file();
        my $fasta_in = BioUtils::FastaIO->new({
                            stream_type => '<',
                            file => $seq_file}
                        );
        
        # open the two output fasta files
        my $dir = $self->get_output_dir();
        my $prefix = $self->get_output_prefix();
        
        my $contam_output = $dir . "/" . $prefix . "_contaminants.fasta";
        my $non_contam_output = $dir . "/" . $prefix . "_non_contaminants.fasta";
    
        my $contam_out = BioUtils::FastaIO->new({
                                stream_type => '>',
                                file => $contam_output}
                            );
        my $non_contam_out = BioUtils::FastaIO->new({
                                    stream_type => '>',
                                    file => $non_contam_output}
                                );
        
        while ( my $seq = $fasta_in->get_next_seq() ) {
            my @values = split /\s/, $seq->get_header();
            my $contam_name = $contam_names_href->{$values[0]};
            if ( $contam_names_href->{$values[0]} ) {
                $seq->set_header($seq->get_header() .
                                 " " .
                                 $contam_names_href->{$values[0]}
                                );
                $contam_out->write_seq($seq);
            }
            else {
                $non_contam_out->write_seq($seq);
            }
        }
        
        return 1;
    }
    
    sub _otu_table_printing {
        my ($self, $contam_names_href) = @_;
        
        my $otu_table_file = $self->get_otu_table();
        my $dir = $self->get_output_dir();
        my $prefix = $self->get_output_prefix();
        
        # open the input OTU table
        open my $OTU_IN, "<", "$otu_table_file" or
            croak("Cannot open file: $otu_table_file");
        
        # open the two output OTU table files
        my $file1 = $dir . "/" . $prefix . "_contaminants_otu_table.txt";
        open my $OTU_CON, ">", "$file1" or
            croak("Cannot open file: $file1");
        
        my $file2 = $dir . "/" . $prefix . "_non_contaminants_otu_table.txt";
        open my $OTU_NON, ">", "$file2" or
            croak("Cannot open file: $file2");
        
        # foreach OTU
        LINE: foreach my $line ( <$OTU_IN> ) {
            chomp $line;
            my @values = split /\t/, $line;
            
            # print header lines to both OTU tables
            if ( $line =~ m/^#/ ) {
                print $OTU_CON $line, "\n";
                print $OTU_NON $line, "\n";
                next LINE;
            }
            
            if ( defined $contam_names_href->{"OTU_" . $values[0]} ) {
                print $OTU_CON $line, "\n";
            }
            else {
                print $OTU_NON $line, "\n";
            }
        }
        
        close($OTU_IN);
        close($OTU_CON);
        close($OTU_NON);
        
        return 1;
    }
}


1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::QC::ContaminantFilter - Identifies and removes sequence contaminants


=head1 VERSION

This document describes BioUtils::QC::ContaminantFilter version 1.0.2


=head1 Included Modules

    Carp qw(carp croak)
    Class::Std::Utils
    Readonly
    List::MoreUtils qw(any)
    MyX::Generic
    BioUtils::FastaSeq 1.0.2
    BioUtils::FastaIO 1.0.2
    File::Temp qw(tempfile)
    version


=head1 SYNOPSIS

    use BioUtils::QC::ContaminantFilter;
    
    # set the parameters
    my $contam_filter = BioUtils::QC::ContaminantFilter->new({
                            blast_db => $blast_db,
                            query_file => $query_file,
                            output_dir => $output_dir,
                            output_prefix => $output_prefix,
                            eval => $eval,
                            perc_iden => $perc_iden,
                            output_fmt => $output_fmt,
                            max_targets => $max_targets,
                            otu_table => $otu_table_file,
                        });
    
    # run the filtering process
    $contam_filter->filter();
  
  
=head1 DESCRIPTION

    After running OTUpipe we want to remove any contaminant sequences such as
    chloroplasts or mithochondria that may confound our abundance estimates.  We
    use this object to help with that.  You can optionally include a OTU table
    file.  Both the input fasta file and the OTU table are split into two output
    files.  One file has all the information about the contaminants and the
    other file has all the information about the non-contaminants.  The
    directory where these output files are located is specified in the params
    file provided by the user.


=head1 METHODS

    # Setters #
    set_blast_db
    set_query_file
    set_output_dir
    set_output_prefix
    set_eval
    set_perc_iden
    set_output_fmt
    set_max_targets
    set_otu_table
    
    # Getters #
    get_params_file
    get_blast_db
    get_query_file
    get_output_dir
    get_output_prefix
    get_eval
    get_perc_iden
    get_output_fmt
    get_max_targets
    get_otu_table
    
    # Others #
    run_filter
    _init
    _run_blast
    _parse_blast_file
    _print_results
    _sequence_printing
    _otu_table_printing

=head1 DIAGNOSTICS

=over

=item C<< Undefined parameter value >>

A required parameter value was not set correctly.

=item C<< Cannot open temp BLAST output file >>

The BLAST output file cannot be opened (probably because it was not created).

=item C<<Cannot open file: [file] >>

A required file cannot be opened.

=back


=head1 CONFIGURATION AND ENVIRONMENT

    BLAST 2.2.25+ or a compatable version must be installed.  The blastn and
    makeblastdb commands must be avaliable as system commands (i.e. they must
    be found in the system PATH variable).


=head1 DEPENDENCIES

 BLAST 2.2.25+


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS DESCRIPTION

=head2 new

    Title: new
    Usage: Sample->new({blast_db => $blast_db,
                        query_file => $query_file,
                        output_dir => $output_dir,
                        output_prefix => $output_prefix,
                        eval => $eval,
                        perc_iden => $perc_iden,
                        output_fmt => $output_fmt,
                        max_targets => $max_targets,
                        otu_table => $otu_table_file,
                        });
    Function: Creates a new BioUtils::QC::ContaminantFilter object
    Returns: BioUtils::QC::ContaminantFilter
    Args: -blast_db => full path and name of the blast database
	Throws: MyX::Generic::Undef::Param
	Comments: e.g. /Home/Me/db/my_blast_db
	See Also: NA

=head2 run_filter

    Title: run_filter
    Usage: $contam_filter->run_filter();
    Function: Filters out OTUs from the fasta and OTU table that are contaminants
    Returns: 1 on successful completion
    Args: NA
	Throws: NA
	Comments: There are four output files: 2 fasta files and 2 otu table files.
              Each pair of output files consists of both a contaminant file and 
              a non-contaminant file.
	See Also: NA
    
=head2 set_blast_db

    Title: set_blast_db
    Usage: $contam_filter->set_blast_db($blast_db);
    Function: Sets the blast database path and name
    Returns: 1 on successful completion
    Args: blast_db => the path and name of the blast database
	Throws: NA
	Comments: set_blast_db checks for the *.nhr blast database file by looking
              for "$blast_db.nhr".  If that file cannot be found set_blast_db
              assumes that the $blast_db value is a fasta file and creates a
              blast database using that file and the makeblastdb system command.
	See Also: NA
    
=head2 set_query_file

    Title: set_query_file
    Usage: $contam_filter->set_query_file($query_file);
    Function: Sets the file with the query sequences
    Returns: 1 on successful completion
    Args: query_file => the query sequence file path
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 set_output_dir

    Title: set_output_dir
    Usage: $contam_filter->set_output_dir($output_dir);
    Function: Sets the output directory path
    Returns: 1 on successful completion
    Args: output_dir => directory file path
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 set_output_prefix

    Title: set_output_prefix
    Usage: $contam_filter->set_output_prefix($output_prefix);
    Function: Sets the output prefix
    Returns: 1 on successful completion
    Args: output_prefix => prefix string
	Throws: NA
	Comments: Output files are named by prepending the prefix
	See Also: NA
    
=head2 set_eval

    Title: set_eval
    Usage: $contam_filter->set_eval($eval);
    Function: Sets the min eval for blast to display
    Returns: 1 on successful completion
    Args: output_dir => eval as the blastn program would accept
	Throws: NA
	Comments: If not specified then blastn default is used
	See Also: NA
    
=head2 set_perc_iden

    Title: set_perc_iden
    Usage: $contam_filter->set_perc_iden($perc_iden);
    Function: Sets the min precent identity for blast to display
    Returns: 1 on successful completion
    Args: perc_iden => percent identity as the blastn program would accept
	Throws: NA
	Comments: If not specified then blastn default is used
	See Also: NA
    
=head2 set_output_fmt

    Title: set_output_fmt
    Usage: $contam_filter->set_output_fmt($output_fmt);
    Function: Sets the blast output format
    Returns: 1 on successful completion
    Args: output_fmt => output format as the blastn program would accept
	Throws: NA
	Comments: If not specified then blastn default is used
	See Also: NA
    
=head2 set_max_targets

    Title: set_max_targets
    Usage: $contam_filter->set_max_targets($max_targets);
    Function: Sets the maximum number of hits for blast to display
    Returns: 1 on successful completion
    Args: max_targets => maximum targets as the blastn program would accept
	Throws: NA
	Comments: If not specified then blastn default is used
	See Also: NA
    
=head2 set_otu_table

    Title: set_otu_table
    Usage: $contam_filter->set_otu_table($otu_table);
    Function: Sets the otu table file path
    Returns: 1 on successful completion
    Args: otu_table => file path to otu table file
	Throws: NA
	Comments: The OTU table is an optional input
	See Also: NA
    
=head2 get_blast_db

    Title: get_blast_db
    Usage: my $blast_db = $contam_filter->get_blast_db();
    Function: Gets the blast database path and name
    Returns: Str
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_query_file

    Title: get_query_file
    Usage: my $query_file = $contam_filter->get_query_file();
    Function: Gets the query file path
    Returns: Str
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_output_dir

    Title: get_output_dir
    Usage: my $output_dir = $contam_filter->get_output_dir();
    Function: Gets the output dir path
    Returns: Str
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_output_prefix

    Title: get_output_prefix
    Usage: my $output_prefix = $contam_filter->get_output_prefix();
    Function: Gets the output prefix string
    Returns: Str
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_eval

    Title: get_eval
    Usage: my $eval = $contam_filter->get_eval();
    Function: Gets the min evalue blast parameter
    Returns: Str
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_perc_iden

    Title: get_perc_iden
    Usage: my $perc_iden = $contam_filter->get_perc_iden();
    Function: Gets the min precent identity blast parameter
    Returns: Int
    Args: NA
	Throws: NA
	Comments: This is a whole number (e.g. 97 is 97% identity)
	See Also: NA
    
=head2 get_output_fmt

    Title: get_output_fmt
    Usage: my $output_fmt = $contam_filter->get_output_fmt();
    Function: Gets the output format blast parameter
    Returns: Int
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_max_targets

    Title: get_max_targets
    Usage: my $max_targets = $contam_filter->get_max_targets();
    Function: Gets the max number of targets (i.e. hits) blast parameter
    Returns: Int
    Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_otu_table

    Title: get_otu_table
    Usage: my $otu_table = $contam_filter->get_otu_table();
    Function: Gets the otu table file path
    Returns: Str
    Args: NA
	Throws: NA
	Comments: May return undef if this parameter is not specified.  It is an
              optional parameter.
	See Also: NA
    
=head2 _init

    Title: _init
    Usage: $contam_filter->_init($arg_href);
    Function: Sets the parameter values based on the arg_href
    Returns: 1 on successful completion
    Args: arg_href => a hash ref with the parameter values
	Throws: NA
	Comments: Acceptable paramter values include:
                blast_db
                query_file
                output_dir
                output_prefix
                eval
                perc_iden
                output_fmt
                max_targets
                otu_table
              All other parameter values are ignored
	See Also: NA
    
=head2 _run_blast

    Title: _run_blast
    Usage: $contam_filter->_run_blast();
    Function: Runs command line blast program using the parameters stored in
              this object.
    Returns: Blast output file
    Args: NA
	Throws: NA
	Comments: The blast output file is a tempfile.  It is returned so it can be
              parsed by _parse_blast_file
	See Also: NA
    
=head2 _parse_blast_file

    Title: _parse_blast_file
    Usage: $contam_filter->_parse_blast_file();
    Function: Parses the blast output
    Returns: Hash ref of contaminant names
    Args: NA
	Throws: NA
	Comments: The contaminant names are used by _print_results
	See Also: NA
    
=head2 _print_results

    Title: _print_results
    Usage: $contam_filter->_print_results();
    Function: Prints the output files
    Returns: 1 on successful completion
    Args: NA
	Throws: NA
	Comments: an OTU table output file is only printed if an otu table parameter
              is defined in the object.  Outupt files are printed in the
              directory specified by output_dir with a prefix string specified
              by output_prefix.
	See Also: _sequence_printing, _otu_table_printing
    
=head2 _sequence_printing

    Title: _sequence_printing
    Usage: $contam_filter->_sequence_printing();
    Function: Prints the sequence output files
    Returns: 1 on successful completion
    Args: NA
	Throws: NA
	Comments: Prints two files: a contaminant file and a non_contaminant file.
              Both are in fasta format.
	See Also: NA
    
=head2 _otu_table_printing

    Title: _otu_table_printing
    Usage: $contam_filter->_otu_table_printing();
    Function: Prints the OTU table output files
    Returns: 1 on successful completion
    Args: NA
	Throws: NA
	Comments: Prints two files: a contaminant file and a non_contaminant file.
              Both are tab delimited otu table text files.
	See Also: NA
    

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to the author


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