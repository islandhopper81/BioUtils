package BioUtils::ConsensusBuilder::FastqConsensus;

use strict;
use warnings;

use Class::Std::Utils;
use Carp;
use Readonly;
use List::MoreUtils qw(any);
use Scalar::Util qw(looks_like_number);
use version; our $VERSION = qv('1.2.0');
use BioUtils::FastqSeq 1.2.0;
use base qw(BioUtils::FastqSeq);  # inherits from FastqSeq

{
    Readonly my $NEW_USAGE => q{ new( {seq => ,
                                       quals_str => ,
                                       [c_score => ,]
                                       } ) };

    # Attributes #
    my %c_score_of;
    
    ###########
    # Setters #
    ###########
    sub set_c_score;
    
    ###########
    # Getters #
    ###########
    sub get_c_score;
    
    ##########
    # Others #
    ##########
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{seq},
                               $arg_href->{quals_str},
            ) {
            MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => $NEW_USAGE,
            );
        }
        
        # Bless a scalar to instantiate an object
        my $new_obj = $class->SUPER::new($arg_href);
        if ( defined $arg_href->{c_score} ) {
            $new_obj->set_c_score($arg_href->{c_score});
        }
        
        return $new_obj;
    }
    
    ###########
    # Setters #
    ###########
    sub set_c_score {
        my ($self, $c_score) = @_;
        
        if ( ! looks_like_number($c_score) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => 'Must be a digit',
                value => $c_score,
            );
        }
        
        $c_score_of{ident $self} = $c_score;
        return 1;
    }
    
    ###########
    # Getters #
    ###########
    sub get_c_score {
        my ($self) = @_;
        
        if ( ! defined $c_score_of{ident $self} ) {
            MyX::Generic::Undef::Param->throw(
                error => 'Getting undefined c_score attribute value',
            )
        }
        
        return $c_score_of{ident $self};
    }
}


1;
__END__

#######
# POD #
#######
=head1 BioUtils::ConsensusBuilder::FastqConsensus

FastqConsensus - A BioUtils::FastqSeq that is specifically named as a consensus
to avoid confusion

=head1 VERSION

This documentation refers to FastqConsensus version 1.2.0.

=head1 Included Modules

    BioUtils::FastqSeq

=head1 Inherit

    BioUtils::FastqSeq

=head1 SYNOPSIS
    
    # NOTE: This is simply a BioUtils::FastqSeq that is renamed explicitly as a
    #       consensus to avoid confusion
    
    use BioUtils::ConsensusBuilder::FastqConsensus;
    my $fastq_con = BioUtils::ConsensusBuilder::FastqConsensus->new(
        {seq => $seq, quals_str => $quals_str}
    );
    
    my $seq_str = $fastq_con->get_seq();
    my $quals_str = $fastq_con->get_quals_str();
    my $quals_aref = $fastq_con->get_quals_aref();
    
    my $new_seq = "ATCG";
    $fastq_con->set_seq($new_seq);

    $fastq_con->set_quals_str('9876');

=head1 DESCRIPTION

FastqConsensus is simply a BioUtils::FastqSeq that is explicitly named as a
consensus to avoid confusion.  See BioUtils::FastqSeq for details on methods.

=head1 METHODS

=over

    new
    set_c_score
    get_c_score
    
=back

    
=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: BioUtils::ConsensusBuilder::FastqConsensus->new(
                {seq => $seq, quals_str => $quals_str}
            );
    Function: Creates a new FastqConsensus object
    Returns: BioUtils::FastqSeq
    Args: -seq => a string representing the sequence
          -quals_str => a string of quality values in Illimina-1.8+ encoding
    Throws: MyX::Generic::Digit::MustBeDigit
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_c_score
    
    Title: set_c_score
    Usage: my $fastqCon = set_c_score(40)
    Function: Sets the consensus score value.  The consensus score is a grade of
              how similar the sequences making up this consensus are.
    Returns: 1 on success
    Args: -c_score => an int or float
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA
    
=head2 get_c_score
    
    Title: get_c_score
    Usage: my $c_score = $fastqCon->get_c_score()
    Function: Gets the consensus score value.  The consensus score is a grade of
              how similar the sequences making up this consensus are.
    Returns: Digit
    Args: NA
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA

    

=head1 BUGS AND LIMITATIONS

No bugs have been reported.


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
