package BioUtils::FastqSeq;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use Scalar::Util qw(looks_like_number);
use MyX::Generic 1.0.11;
use BioUtils::FastaSeq 1.0.11;

use version; our $VERSION = qv('1.0.11');

{
    Readonly my $ASCII_OFFSET => 33;
    Readonly my $ASCII_MAX => 126;  # This is a limitation of perls chr command
    Readonly my $NEW_USAGE => q{ new( {header =>, seq =>, quals_str => } ) };
    
    # Class Attributes
    my %header_of;
    my %seq_of;
    my %quals_str_of;
    
    # Public Class Methods
    sub get_header;
    sub get_id;
    sub get_seq;
    sub get_quals_str;
    sub get_quals_aref;
    sub get_qual_at;
    sub set_header;
    sub set_seq;
    sub set_quals_str;
    sub to_FastaSeq;
    sub trim_back;
    sub trim_front;
    sub rev_comp;
    sub rev;
    sub comp;
    
    # Private Class Methods
    sub _dec_to_encoding;
    sub _encoding_to_dec;
    
    # Constructor
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{seq}, $arg_href->{quals_str} ) {
            MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => $NEW_USAGE,
            );
        }
        
        # Initialize the objects attributes
        $header_of{ident $new_obj} = $arg_href->{header};  # an optional parameter
        $seq_of{ident $new_obj} = $arg_href->{seq};
        $quals_str_of{ident $new_obj} = $arg_href->{quals_str};
    
        return $new_obj;
    }
    
    sub get_header {
        my ($self) = @_;
        
        # check if the header is defined
        if ( ! defined $header_of{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                error => 'Undefined header'
            );
        }
        
        return $header_of{ident $self};
    }
    
    sub get_id {
        my ($self) = @_;
        
        # check if the header is defined
        if ( ! defined $header_of{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                error => 'Undefined header'
            );
        }
        
        # Mach the header up until the first white space
        if ( $header_of{ident $self} =~ m/(\S+)/ ) {
            return $1;
        }
        
        return undef;
    }

    sub get_seq {
        my ($self) = @_;
        return $seq_of{ident $self};
    }
        
    sub get_quals_str {
        my ($self) = @_;
        return $quals_str_of{ident $self};
    }
    
    sub get_quals_aref {
        my ($self) = @_;
        
        # First split the quals string and then change each value to a decimal
        #my @quals = map { _encoding_to_dec($_) } split //, $quals_str_of{ident $self};
        my @quals = split //, $quals_str_of{ident $self};
        
        return \@quals;
    }
    
    sub get_qual_at {
        my ($self, $index) = @_;
        
        # Check if the index is out of bounds
        if ( $index < 0 or $index >= length $quals_str_of{ident $self} ) {
            MyX::Generic::OutOfBounds->throw(
                error => 'Index out of bounds',
                index => $index,
            );
        }
        
        my $qual = _encoding_to_dec(substr $quals_str_of{ident $self}, $index, 1);
        
        return $qual;
    }
    
    sub set_header($) {
        my ($self, $header) = @_;
        $header_of{ident $self} = $header;
        return 1;
    }
        
    sub set_seq($) {
        my ($self, $seq) = @_;
        $seq_of{ident $self} = $seq;
        return 1;
    }
    
    sub set_quals_str($){
        my ($self, $quals_str) = @_;
        $quals_str_of{ident $self} = $quals_str;
        return 1;
    }
    
    sub to_FastaSeq {
        my ($self) = @_;
        
        my $header = $self->get_header();
        $header =~ s/@//g;
        my $fasta_seq = BioUtils::FastaSeq->new({
                            header => $header,
                            seq => $self->get_seq()
                        });
        
        return $fasta_seq;
    }
    
    sub trim_front {
        my ($self, $len) = @_;
        
        if ( ! defined $len ) {
            my $seq_obj = BioUtils::FastqSeq->new({
                header => $self->get_header(),
                seq => '',
                quals_str => '',
            });
        
            return $seq_obj;
        }
        
        if ( ! looks_like_number($len) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => "trim_front requires digit > 0",
                value => $len,
            );
        }
        
        if ( $len < 0 ) {
            MyX::Generic::Digit::TooSmall->throw(
                error => "trim_front requires digit > 0",
                value => $len,
                MIN => 0,
            )
        }
        
        my $seq = $self->get_seq();
        my $qual = $self->get_quals_str();
        my $keep_seq = substr $seq, -((length $seq) - $len);
        my $keep_qual = substr $qual, -((length $seq) - $len);
        my $trim_seq = substr $seq, 0, $len;
        my $trim_qual = substr $qual, 0, $len;
        
        my $trimmed_seq_obj = BioUtils::FastqSeq->new({
            header => $self->get_header(),
            seq => $trim_seq,
            quals_str => $trim_qual,
        });
        
        $self->set_seq($keep_seq);
        $self->set_quals_str($keep_qual);
        
        return $trimmed_seq_obj;
    }
    
    sub trim_back {
        my ($self, $len) = @_;
        
        if ( ! defined $len ) {
            my $seq_obj = BioUtils::FastqSeq->new({
                header => $self->get_header(),
                seq => '',
                quals_str => '',
            });
            
            return $seq_obj;
        }
        
        if ( ! looks_like_number($len) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => "trim_back requires digit > 0",
                value => $len,
            );
        }
        
        if ( $len < 0 ) {
            MyX::Generic::Digit::TooSmall->throw(
                error => "trim_back requires digit > 0",
                value => $len,
                MIN => 0,
            )
        }
        
        my $seq = $self->get_seq();
        my $qual = $self->get_quals_str();
        my $trim_seq = substr $seq, -$len;
        my $trim_qual = substr $qual, -$len;
        my $keep_seq = substr $seq, 0, (length $seq) - $len;
        my $keep_qual = substr $qual, 0, (length $seq) - $len;
        
        my $trimmed_seq_obj = BioUtils::FastqSeq->new({
            header => $self->get_header(),
            seq => $trim_seq,
            quals_str => $trim_qual,
        });
        
        $self->set_seq($keep_seq);
        $self->set_quals_str($keep_qual);
        
        return $trimmed_seq_obj;
    }
    
    sub rev_comp {
        my ($self) = @_;
        
        $self->rev();
        $self->comp();
        
        return 1;
    }
    
    sub rev {
        my ($self) = @_;
        
        $self->set_seq(scalar reverse $self->get_seq());
        $self->set_quals_str(scalar reverse $self->get_quals_str());
        
        return 1;
    }
    
    sub comp {
        my ($self) = @_;
        
        my $str = $self->get_seq();
        $str =~ tr/ACGTacgt/TGCAtgca/;
        $self->set_seq($str);
        
        return 1;
    }
    
    sub _dec_to_encoding {
        my ($dec)  = @_;
        
        # encode the decimal in Illumina-1.8+ encoding format
        my $encoding = chr ($dec + $ASCII_OFFSET);
        
        return $encoding;
    }
    
    sub _encoding_to_dec {
        my ($encoding) = @_;
        
        # decode the Illumina-1.8+ encoding to a decimal
        my $dec = (ord $encoding) - $ASCII_OFFSET;
    }
    
    sub DESTROY {
        my ($self) = @_;
        
        delete $header_of{ident $self};
        delete $seq_of{ident $self};
        delete $quals_str_of{ident $self};
    }
}

1;
__END__


#######
# POD #
#######
=head1 BioUtils::FastqSeq

FastqSeq - A data structure to store a sequence string and quality string

=head1 VERSION

This documentation refers to FastqSeq version 1.0.11.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(croak)
    Scalar::Util qw(looks_like_number)
    MyX::Generic 1.0.11
    BioUtils::FastaSeq 1.0.11

=head1 Inherit

    NA

=head1 SYNOPSIS

    use BioUtils::FastqSeq;
    my $fastq_seq = BioUtils::FastqSeq->new( {
                        header => 'Seq1',
                        seq => 'ATCG',
                        quals_str => '9876'
                    } );
    
    my $header = $fastq_seq->get_header();
    my $id = $fastq_seq->get_id();
    my $seq_str = $fastq_seq->get_seq();
    my $quals_str = $fastq_seq->get_quals_str();
    my $quals_aref = $fastq_seq->get_quals_aref();
    my $qual_at_index = $fastq_seq->get_qual_at($index);
    my $encoding = $fastq_seq->get_qual_encoding();
    
    my $new_seq = "ATCG";
    my $new_quals = "efcf";
    $fastq_seq->set_quals_encoding('L')
    $fastq_seq->set_seq($new_seq);
    $fastq_seq->set_quals($new_quals);
    $fastq_seq->set_header("seq1");
    
    # convert to FastaSeq Object
    my $fasta_seq = $fastq_seq->to_FastaSeq()
    
    # trim
    my $trimmed_portion = $fastq_seq->trim_front(2);
    my $trimmed_portion = $fastq_seq->trim_back(2);
    
    # reverse the sequence
    $fasta_seq->rev();
    
    # complement the sequence
    $fasta_seq->comp();
    
    # reverse complement the sequence
    $fasta_seq->rev_comp();
    

=head1 DESCRIPTION

FastqSeq is a data structure designed to store a biological sequence and its
corresponding quality values.  The sequence and quality values are stored as
strings.

There are several differnt quality value encodings that have been used to
represent base quality scores.  This object assumes quality values are encoded
in the Illumina-1.8+ format.  For a short description of encodings see the
Wikipedia page FASTQ Format (http://en.wikipedia.org/wiki/FASTQ_format).

Various getter and setter subroutines are provided.  Note that the quality
values can be accessed as a string or array reference.

=head1 METHODS

=over

    new
    get_header
    get_id
    get_seq
    get_quals_str
    get_quals_aref
    get_qual_at
    set_header
    set_seq
    set_quals_str
    to_FastaSeq
    trim_front
    trim_back
    rev
    comp
    rev_comp
    _dec_to_encode
    _encode_to_dec
    DESTROY
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: BioUtils::FastqSeq->new({
                seq => $seq,
                quals_str => $quals_str,
                header => $header
            });
    Function: Creates a new FastqSeq object
    Returns: BioUtils::FastqSeq
    Args: -seq => a string representing the sequence
          -quals_str => a string of quality values in Illimina-1.8+ encoding
          -[header] => a string representing the header
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA

=head2 get_header

    Title: get_header
    Usage: my $header = $my_fastq_seq->get_header();
    Function: Gets the header attribute if defined
    Returns: String
    Args: None
    Throws: MyX::Generic::Undef::Attribute
    Comments: NA
    See Also: NA

=head2 get_id

    Title: get_id
    Usage: my $id = $my_fastq_seq->get_id();
    Function: Gets the id by parsing the header attribute (if defined)
    Returns: String
    Args: None
    Throws: MyX::Generic::Undef::Attribute
    Comments: NA
    See Also: NA

=head2 get_seq

    Title: get_seq
    Usage: my $seq = $my_fastq_seq->get_seq();
    Function: Gets the sequence string stored in the FastqSeq object
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_quals_str

    Title: get_quals_str
    Usage: my $quals_str = $my_fastq_seq->get_quals_str();
    Function: Gets the string of quality values in Illumina-1.8+ encoding
    Returns: String
    Args: None
    Throws: NA
    Comments: See http://en.wikipedia.org/wiki/FASTQ_format for encoding
              information.
    See Also: NA
    
=head2 get_quals_aref

    Title: get_quals_aref
    Usage: my $quals_aref = $my_fastq_seq->get_quals_aref();
    Function: Gets the array reference of the quality values stored in the
              FastqSeq object.  The values in the array are in decimal format
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: Values in the array are in decimal format.
    See Also: NA
    
=head2 get_qual_at

    Title: get_qual_at
    Usage: my $qual_at_index = $my_fastq_seq->get_qual_at($index);
    Function: Gets the quality value at the given index in decimal format
    Returns: int
    Args: -index => index of quality value to return
    Throws: NA
    Comments: Return value is in decimal format.
    See Also: NA

=head2 set_header

    Title: set_header
    Usage: $my_fastq_seq->set_header("seq1");
    Function: Sets the header value in the FastqSeq object
    Returns: 1 on successful completion
    Args: -seq => a string representing the header
    Throws: NA
    Comments: NA
    See Also: NA
   
=head2 set_seq

    Title: set_seq
    Usage: $my_fastq_seq->set_seq("ATCG");
    Function: Sets the sequence value in the FastqSeq object
    Returns: 1 on successful completion
    Args: -seq => a string representing the sequence
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_quals_str
    
    Title: set_quals_str
    Usage: $my_fastq_seq->set_quals_str($quals_str);
    Function: Sets the quality string for the sequence in the FastqSeq object
    Returns: 1 on successful completion
    Args: -quals_str => a string of Illumina-1.8+ encoded quality values
    Throws: NA
    Comments: See http://en.wikipedia.org/wiki/FASTQ_format for encoding
              information.
    See Also: NA
    
=head2 to_FastaSeq
    
    Title: to_FastaSeq
    Usage: my $fasta_seq = $my_fastq_seq->to_FastaSeq();
    Function: Converts the FastqSeq object as a FastaSeq object
    Returns: BioUtils::FastaSeq
    Args: NA
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastaSeq
    
=head2 trim_front
    
    Title: trim_front
    Usage: my $trimmed_portion = $my_fastq_seq->trim_front($num);
    Function: Trims X bases off the front of the FastqSeq object
    Returns: BioUtils::FastqSeq
    Args: -num => the number of bases to trim off the front
    Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
    Comments: The part that is trimmed off is returned.  To ignore/trash that
              simply call the method like $my_fastq_seq->trim_front(2) (i.e.
              don't store the return value).
    See Also: NA
    
=head2 trim_back
    
    Title: trim_back
    Usage: my $trimmed_portion = $my_fastq_seq->trim_back($num);
    Function: Trims X bases off the back of the FastqSeq object
    Returns: BioUtils::FastqSeq
    Args: -num => the number of bases to trim off the back
    Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
    Comments: The part that is trimmed off is returned.  To ignore/trash that
              simply call the method like $my_fastq_seq->trim_back(2) (i.e.
              don't store the return value).
    See Also: NA
    
=head2 rev
    
    Title: rev
    Usage: $my_fastq_seq->rev();
    Function: Reverse the sequence
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 comp
    
    Title: comp
    Usage: $my_fastq_seq->comp();
    Function: Compelements the sequence
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 rev_comp
    
    Title: rev_comp
    Usage: $my_fastq_seq->rev_comp();
    Function: Reverse compelements the sequence
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 _dec_to_encoding
    
    Title: _dec_to_encoding
    Usage: _dec_to_encoding($dec);
    Function: Encodes the given decimal in Illiina-1.8+ encoding format
    Returns: char
    Args: -dec => an interger
    Throws: NA
    Comments: NA
    See Also: NA
  
=head2 _encoding_to_dec
    
    Title: _encoding_to_dec
    Usage: _encoding_to_dec($encoding);
    Function: Decodes the Illumina-1.8+ format encoding into a decimal
    Returns: int
    Args: -encoding => a char in Illumina-1.8+ encoding format
    Throws: NA
    Comments: NA
    See Also: NA
  
=head2 DESTROY
    
    Title: DESTROY
    Usage: $my_fastq_seq->DESTROY();
    Function: Deteles attribute hashes
    Returns: NA
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head1 CONFIGURATION AND ENVIRONMENT

FastaSeq requires no configuration files or environment variables.


=head1 DEPENDENCIES

    None.


=head1 INCOMPATIBILITIES

    None reported.


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-fastaseq@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.


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




