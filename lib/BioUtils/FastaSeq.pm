package BioUtils::FastaSeq;

use warnings;
use strict;

use Carp qw(croak carp);
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Scalar::Util qw(looks_like_number);
use MyX::Generic 0.0.1;
use BioUtils::MyX::Fasta;
use UtilSY qw(:all);
use version; our $VERSION = qv('1.2.1');

{
    Readonly my $NEW_USAGE => q{ new( {header =>, seq => } ) };
    Readonly::Hash my %CODON_TBL => ('TCA'=>'S','TCC'=>'S','TCG'=>'S',
        'TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y',
        'TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_',
        'TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P',
        'CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q',
        'CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I',
        'ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T',
        'ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S',
        'AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V',
        'GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D',
        'GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G',
        'GGT'=>'G');
    
    # Class Attributes
    my %header_of;
    my %seq_of;
    
    # Public Class Methods
    sub get_header;
    sub get_id;
    sub get_seq;
    sub set_header;
    sub set_seq;
    sub trim_front;
    sub trim_back;
    sub substr;
    sub rev;
    sub comp;
    sub rev_comp;
    sub translate;
    
    # Private Class Methods

    
    # Constructor
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{seq} ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Initialize the objects attributes
        $header_of{ident $new_obj} = $arg_href->{header};  # an optional parameter
        $seq_of{ident $new_obj} = $arg_href->{seq};
    
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
    
    sub set_header {
        my ($self, $header) = @_;
        $header_of{ident $self} = $header;
        return 1;
    }
    
    sub set_seq {
        my ($self, $seq) = @_;
        $seq_of{ident $self} = $seq;
        return 1;
    }
    
    sub trim_front {
        my ($self, $len) = @_;
        
        if ( ! defined $len ) {
            my $seq_obj = BioUtils::FastaSeq->new({
                header => $self->get_header(),
                seq => '',
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
        my $keep_seq = CORE::substr $seq, -((length $seq) - $len);
        my $trim_seq = CORE::substr $seq, 0, $len;
        
        my $trimmed_seq_obj = BioUtils::FastaSeq->new({
            header => $self->get_header(),
            seq => $trim_seq,
        });
        
        $self->set_seq($keep_seq);
        
        return $trimmed_seq_obj;
    }
    
    sub trim_back {
        my ($self, $len) = @_;
        
        if ( ! defined $len ) {
            my $seq_obj = BioUtils::FastaSeq->new({
                header => $self->get_header(),
                seq => '',
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
        my $trim_seq = CORE::substr $seq, -$len;
        my $keep_seq = CORE::substr $seq, 0, (length $seq) - $len;
        
        my $trimmed_seq_obj = BioUtils::FastaSeq->new({
            header => $self->get_header(),
            seq => $trim_seq,
        });
        
        $self->set_seq($keep_seq);
        
        return $trimmed_seq_obj;
    }
    
    sub substr {
        my ($self, $start, $end) = @_;
        
        if ( ! defined $start ) {
            MyX::Generic::Undef::Param->throw(
                error => 'start parameter not defined',
                usage => '$fasta_seq->substr(10, 20)'
            );
        }
        if ( ! defined $end ) {
            MyX::Generic::Undef::Param->throw(
                error => 'end parameter not defined',
                usage => '$fasta_seq->substr(10, 20)'
            );
        }
        if ( ! looks_like_number($start) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => "substr start requires digit > 0",
                value => $start,
            );
        }
        if ( ! looks_like_number($end) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => "substr end requires digit > 0",
                value => $end,
            );
        }
        if ( $start < 0 ) {
             MyX::Generic::Digit::TooSmall->throw(
                error => "start requires digit > 0",
                value => $start,
                MIN => 0,
            )
        }
        if ( $end < 0 ) {
             MyX::Generic::Digit::TooSmall->throw(
                error => "end requires digit > 0",
                value => $end,
                MIN => 0,
            )
        }
        if ( $start > $end ) {
            MyX::Generic::Digit::TooBig->throw(
                error => "end is larger than start",
                value => $end,
                MAX => "must be larger than start parameter"
            )
        }
        
        my $seq = $self->get_seq();
        my $len = $end - $start + 1;
        my $seq_substr = CORE::substr $seq, $start, $len;
        
        my $seq_substr_obj = BioUtils::FastaSeq->new({
            header => $self->get_header(),
            seq => $seq_substr,
        });
        
        return $seq_substr_obj;
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
        return 1;
    }
    
    sub comp {
        my ($self) = @_;
        
        my $str = $self->get_seq();
        $str =~ tr/ACGTacgt/TGCAtgca/;
        $self->set_seq($str);
        
        return 1;
    }
    
    sub translate {
        my ($self, $is_gene) = @_;
        
        # set the default is_gene value to TRUE
        if ( ! defined $is_gene ) {
            $is_gene = 1; # TRUE
        }
        $is_gene = to_bool($is_gene);
        
        my $str = $self->get_seq();
        my $len = length $str;
        my $aa_str = "";
        my $aa;
        my $codon;
        
        # make sure this sequence is divisible by 3.  Only send a warning if it
        # is not
        if ( $len % 3 != 0 ) {
            warn "Seq (" . $self->get_id() . ") is not divisible by 3.  Trailing bases are ignored."
        }
        
        for ( my $i = 0; $i < $len-2; $i = $i + 3 ) {
            $codon = CORE::substr $str, $i, 3;
            $codon = uc $codon;
            
            if ( $aa = $CODON_TBL{$codon} ) {
                $aa_str .= $aa;
            }
            elsif ( $codon =~ m/N/g ) {
                # Any codon that has an N is an unknown codon
                $aa_str .= 'X';
            }
            else {
                BioUtils::MyX::Fasta::BadCodon->throw(
                    error => "Bad codon: $codon",
                    codon => $codon
                )
            }
        }
        
        # remove any trailing stop codon symbols (ie "_")
        $aa_str =~ s/_$//;
        
        # make sure the first codon is M
        # "Alternate start codons are still translated as Met when they are at
        # the start of a protein (even if the codon encodes a different amino
        # acid otherwise)"
        # -- wikipedia "Start Codon"
        # There are two conditions in which the first codon is M:
        # 1. if it starts with an unknown codon (ie the codon has an N in it)
        # 2. if the caller specifies the sequence is not a gene using the
        #    is_gene parameter
        if ( $aa_str !~ m/^M|X/ and $is_gene ) {
            $aa_str =~ s/^./M/;
        }
        
        $self->set_seq($aa_str);
        
        return 1;
    }
    
    sub DESTROY {
        my ($self) = @_;
        
        delete $header_of{ident $self};
        delete $seq_of{ident $self};
    }
}


1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::FastaSeq - A data structure for sequences


=head1 VERSION

This document describes BioUtils::FastaSeq version 1.2.1


=head1 SYNOPSIS

    use BioUtils::FastaSeq;

    my $fasta_seq = BioUtils::FastaSeq->new({
        header => 'Seq1',
        seq => 'ATCG'
    });
    
    my $header = $fasta_seq->get_header();
    my $id = $fasta_seq->get_id();
    my $seq_str = $fasta_seq->get_seq();
    
    $fasta_seq->set_seq("ATCT");
    $fasta_seq->set_header("new_header");
    
    # trim
    my $trimmed_portion = $fasta_seq->trim_front(2);
    my $trimmed_portion = $fasta_seq->trim_back(2);
    
    # get a substr of the sequence
    my $substr = $fasta_seq->get_substr(2,4);
    
    # reverse the sequence
    $fasta_seq->rev();
    
    # complement the sequence
    $fasta_seq->comp();
    
    # reverse complement the sequence
    $fasta_seq->rev_comp();
    
    # translate the sequence
    $fasta_seq->translate();
  
  
=head1 DESCRIPTION

    BioUtils::FastaSeq is a very basic data structure for storing and retrieving
    sequence information.  For speed enhancements it does very few checks for
    correctness.


=head1 METHODS

    get_header
    get_id
    get_seq
    set_header
    set_seq
    trim_front
    trim_back
    substr
    rev
    comp
    rev_comp
    translate


=head1 DIAGNOSTICS

=over

=item C<< Undefined parameter value >>

A parameter value was either undefined or was missing from the arguments

=item C<< Undefined header >>

Trying to access a header value that is undefined.  A sequence header is not
a required attribute for a BioUtils::FastaSeq.

=back


=head1 CONFIGURATION AND ENVIRONMENT

BioUtils::FastaSeq requires no configuration files or environment variables.


=head1 DEPENDENCIES

    Carp qw(croak carp)
    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Scalar::Util qw(looks_like_number)
    MyX::Generic 1.2.1
    BioUtils::MyX::Fasta
    version our $VERSION = qv('1.2.1')
    UstilSY qw(:all)


=head1 INCOMPATIBILITIES

    None reported.


=head1 METHODS DESCRIPTION

=head2 new

    Title: new
    Usage: BioUtils::FastaSeq->new({header => "my_seq", seq => "ATCG"});
    Function: Creates a new BioUtils::FastaSeq object
    Returns: BioUtils::FastaSeq
    Args: -[header] => a string representing the header information
          -seq => a string representing the sequence
    Throws: MyX::Generic::Undef::Param
    Comments: The header argument is not required
    See Also: NA


=head2 get_header

    Title: get_header
    Usage: my $header = $my_fasta_seq->get_header();
    Function: Gets the header attribute if defined
    Returns: String
    Args: None
    Throws: MyX::Generic::Undef::Attribute
    Comments: NA
    See Also: NA

=head2 get_id

    Title: get_id
    Usage: my $id = $my_fasta_seq->get_id();
    Function: Gets the id by parsing the header attribute (if defined)
    Returns: String
    Args: None
    Throws: MyX::Generic::Undef::Attribute
    Comments: NA
    See Also: NA

=head2 get_seq

    Title: get_seq
    Usage: my $seq = $my_fasta_seq->get_seq();
    Function: Gets the sequence string stored in the BioUtils::FastaSeq object
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_header

    Title: set_header
    Usage: $my_fasta_seq->set_header("my_seq");
    Function: Sets the header value in the BioUtils::FastaSeq object
    Returns: 1 on successful completion
    Args: -header => a string representing the header
    Throws: NA
    Comments: NA
    See Also: NA
   
=head2 set_seq

    Title: set_seq
    Usage: $my_fasta_seq->set_seq("ATCG");
    Function: Sets the sequence value in the BioUtils::FastaSeq object
    Returns: 1 on successful completion
    Args: -seq => a string representing the sequence
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 trim_front
    
    Title: trim_front
    Usage: my $trimmed_portion = $my_fasta_seq->trim_front($num);
    Function: Trims X bases off the front of the FastaSeq object
    Returns: BioUtils::FastaSeq
    Args: -num => the number of bases to trim off the front
    Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
    Comments: The part that is trimmed off is returned.  To ignore/trash that
              simply call the method like $my_fasta_seq->trim_front(2) (i.e.
              don't store the return value).
    See Also: NA
    
=head2 trim_back
    
    Title: trim_back
    Usage: my $trimmed_portion = $my_fasta_seq->trim_back($num);
    Function: Trims X bases off the back of the FastaSeq object
    Returns: BioUtils::FastaSeq
    Args: -num => the number of bases to trim off the back
    Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
    Comments: The part that is trimmed off is returned.  To ignore/trash that
              simply call the method like $my_fasta_seq->trim_back(2) (i.e.
              don't store the return value).
    See Also: NA
    
=head2 substr

    Title: substr
    Usage: my $substr_seq = $my_fasta_seq->substr($start, $end);
    Function: Get a substring as a new FastaSeq object
    Returns: BioUtils::FastaSeq
    Args: -start => start position of substring
          -end => end position of substring
    Throws: MyX::Generic::Undef:Param
            MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
            MyX::Generic::Digit::TooBig
    Comments: The start and end positions assume the first base is at position
              0.  So if the sequence is "ATCGATCG" and the start position is 2
              and the end position is 5 the returned BioUtils::FastaSeq object
              has a sequence of CGA
    See Also: NA
    
=head2 rev
    
    Title: rev
    Usage: $my_fasta_seq->rev();
    Function: Reverse the sequence
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 comp
    
    Title: comp
    Usage: $my_fasta_seq->comp();
    Function: Compelements the sequence
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 rev_comp
    
    Title: rev_comp
    Usage: $my_fasta_seq->rev_comp();
    Function: Reverse compelements the sequence
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 translate
    
    Title: translate
    Usage: $my_fasta_seq->translate(is_gene);
    Function: Translates a nucleotide sequence to amino acids
    Returns: 1 on successful completion
    Args: -is_gene => boolean specifing if sequence is a gene (DEFAULT: T)
    Throws: BioUtils::MyX::Fasta::BadCodon
    Comments: This assumes the nucleotide sequence stored in this object is
              alread in frame.  It also throws a warning if there are trailing
              bases.  Those trailing bases are ignored (ie dropped, ie not
              translated).  This function also assumes the sequence is a
              nucleotide sequence.  If it is not you will likely throw a
              BioUtils::MyX::Fasta::BadCodon error.  Any trailing stop codon
              values are removed.  So if your protein is "AIN_" the "_"
              character is automatically removed.
              
              "Alternate start codons are still translated as Met when they are
              at the start of a protein (even if the codon encodes a different
              amino acid otherwise)"
              -- wikipedia "Start Codon"
              
              To get the amino acid sequence without changing the first codon to
              "M" (ie without assuming th sequence is a gene) set the "is_gene"
              parameter to "F".
    See Also: NA

=head2 DESTROY
    
    Title: DESTROY
    Usage: $my_fasta_seq->DESTROY();
    Function: Deteles attribute hashes
    Returns: NA
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA


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
