package BioUtils::FastaSeq;

use warnings;
use strict;

use Carp qw(croak carp);
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use MyX::Generic 1.0.1;
use version; our $VERSION = qv('1.0.1');

{
    Readonly my $NEW_USAGE => q{ new( {header =>, seq => } ) };
    
    # Class Attributes
    my %header_of;
    my %seq_of;
    
    # Public Class Methods
    sub get_header;
    sub get_id;
    sub get_seq;
    sub set_header;
    sub set_seq;
    
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
    
    sub set_seq($) {
        my ($self, $seq) = @_;
        $seq_of{ident $self} = $seq;
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

This document describes BioUtils::FastaSeq version 1.0.1


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

    None.


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
