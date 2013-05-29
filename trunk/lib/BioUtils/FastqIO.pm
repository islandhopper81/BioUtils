package BioUtils::FastqIO;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use MyX::Generic 1.0.2;
use BioUtils::FastqSeq 1.0.2;

use version; our $VERSION = qv('1.0.2');

{
    # Global variables
    Readonly my $NEW_USAGE => q{ new( {file => , stream_type => < | >, } ) };
    
    # Class Attributes
    my %stream_type_of;
    my %file_of;
    my %fh_of;
    
    # Public Class Methods
    sub get_next_seq;
    sub write_seq;
    
    # Private Class Methods
    sub _open_fh;
    sub DESTROY;
    
    # Constructor
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Check that the required parameters are defined (only seq and quals_str)
        if ( any { !defined $_ } $arg_href->{file}, $arg_href->{stream_type} ) {
            MyX::Generic::Undef::Param->throw(
                                                   error => "Undefined parameter",
                                                   usage => $NEW_USAGE,
                                                   );
        }
        
        # Initialize the objects attributes
        $stream_type_of{ident $new_obj} = $arg_href->{stream_type};
        $file_of{ident $new_obj} = $arg_href->{file};
        
        # Open and store the file handle
        $fh_of{ident $new_obj} = $new_obj->_open_fh($arg_href->{file});
    
        return $new_obj;
    }

    sub get_next_seq {
        my ($self) = @_;
        
        my $fh = $fh_of{ident $self};  # stored to avoid repeated calls to ident
        
        # Check for some error in the file handle... do I really want to do these each time I call get_next_seq??
        if ( eof $fh == 1 ) {
            return 0;  # Indicates the end of the file
        }
        
        # Store the lines from the file
        my $header = readline $fh;
        my $seq = readline $fh;
        readline $fh; # This is the quals header which I don't really need
        my $quals_str = readline $fh;
        chomp ($header, $seq, $quals_str);
        
        # remove @ symbol from header if needed
        if ( $header =~ m/@(.*)/ ) {
            $header = $1;
        }
        
        # Create a new BioUtils::FastqSeq
        my $fastq_seq = undef;
        eval {
            $fastq_seq = BioUtils::FastqSeq->new({
                            header => $header,
                            seq => $seq,
                            quals_str => $quals_str});
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            $e->rethrow();
        }        
        
        return $fastq_seq;
    }
    
    sub write_seq {
        my ($self, $fastq_seq) = @_;
        
        # check if the object is a fastq_seq
        
        # check if the values are all defined
        
        # Get the header
        my $header = undef;
        eval {
            $header = $fastq_seq->get_header();
        };
        if ( my $e = MyX::Generic::Undef::Attribute->caught() ) {
            MyX::Generic::Undef::Attribute->throw(
                error => 'Printing a BioUtils::FastqSeq without a header'
            );
        }
        
        # Get the base sequence string
        my $seq = $fastq_seq->get_seq();
        
        # Get the quals string
        my $quals_str = $fastq_seq->get_quals_str();
        
        # Print the gathered information
        my $fh = $fh_of{ident $self};
        print $fh   qw{@}, $header, "\n",
                                    $seq, "\n",
                                    qw{+}, $header, "\n",
                                    $quals_str, "\n"
                                    ;
        
        return 1;
    }
    
    sub _open_fh {
        my ($self, $file) = @_;
        
        open my $fh, $stream_type_of{ident $self}, $file or
            croak("Cannot open file: $file\nERROR: $!\n");

        return $fh;
    }
    
    sub DESTROY {
        my ($self) = @_;
        
        delete $stream_type_of{ident $self};
        delete $file_of{ident $self};
        delete $fh_of{ident $self};
    }
}

1;
__END__


=head1 NAME

BioUtils::FastqIO - An object for parsing FASTQ formated sequence files

=head1 VERSION

This documentation refers to BioUtils::FastqIO version 1.0.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(croak)
    MyX::Generic
    BioUtils::FastqSeq

=head1 Inherit

    NA

=head1 SYNOPSIS

    use BioUtils::FastqIO;
    
    # Build a parser for reading FASTQ files
    my $in = BioUtils::FastqIO->new( {stream_type => '<', file => $file} );
    
    # Read in all the sequences
    while ( my $seq = $in->get_next_seq() ) {
        # Do something with $seq
    }
    
    
    # Build a parser for writing FASTQ files
    my $out = BioUtils::FastqIO->new( {stream_type => '>', file => $file} );
    
    # Write a BioUtils::FastqSeq object to $out
    my $fastq_seq = BioUtils::FastqSeq->new( {
                        header => 'Seq1',
                        seq => 'ATCG',
                        quals_str => '9999'}
                    );
    $out->write_seq($fastq_seq);
    
    

=head1 DESCRIPTION

BioUtils::FastqIO is a classes used to parse FASTQ formated files and write fastq
sequences to FASTQ formated files.  A FASTQ formated file has entries that span
four lines as show in the following example:

@Seq1
ATCG
+Seq1
9999

The first line is the headers line and is started with the '@' symbol.  This
line may include any type of characters up until the newline character.  It
should be a unique identifier of its sequence inside its file.

The second line is of course the sequence.  Right now there are no restrictions
on the types of characters that can be used in a sequence.  If the user has
restrictions on chacters or is worried about illegal characters, they can
look for illegal characters before opperating on/with a character.  In the
future a robust parser may be provided that will do this while reading in a
FASTQ file.

The thrid line is the same as the sequence header (first line) execpt it begins
with the '+' symbol.  Some FASTQ formated files contain only a '+' at this line
to save on file size.  This format can be parsed by the parser, but on writing
a BioUtils::FastqSeq to a file it will contain the full header preceeded by the
'+' symbol.

The fourth and final line are the quality values.  Quality values are encoded
under the Illumina-1.8+ quality value system.  For a short description of
encodings see the Wikipedia page FASTQ Format
(L<en.wikipedia.org/wiki/FASTQ_format>).  The quality value line should be
the same length as the sequence string.  In other words each character in the
quality value string encodes a quality value score.  

=head1 METHODS

    new
    get_next_seq
    write_seq
    _open_fh
    DESTROY


=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: BioUtils::FastqIO->new({stream_type => '<', file => $file});
    Function: Creates a new BioUtils::FastqIO object
    Returns: BioUtils::FastqIO
    Args: -stream_type => '<' for input OR '>' for output
          -file => a string path where the input OR output file is located
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA

=head2 get_next_seq
    
    Title: get_next_seq
    Usage: while ( my $fastq_seq = $in_parser->get_next_seq() ) {...}
    Function: Gets the next sequence in a FASTQ formated file
    Returns: BioUtils::FastqSeq
    Args: NA
    Throws: MyX::Generic::Undef::Param
    Comments: This function will most commonly be encased in a while loop.  
    See Also: NA

=head2 write_seq
    
    Title: write_seq
    Usage: $parser_out->write_seq($fastq_seq);
    Function: Writes a BioUtils::FastqSeq object to a file in FASTQ format
    Returns: 1 on success
    Args: -fastq_seq => a BioUtils::FastqSeq object
    Throws: MyX::Generic::Undef::Attribute
    Comments: NA
    See Also: NA

=head2 DESTROY
    
    Title: DESTROY
    Usage: $my_fastq_seq->DESTROY();
    Function: Deletes attribute hashes
    Returns: NA
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA


=head1 CONFIGURATION AND ENVIRONMENT
  
BioUtils::FastqIO requires no configuration files or environment variables.


=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
   

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
