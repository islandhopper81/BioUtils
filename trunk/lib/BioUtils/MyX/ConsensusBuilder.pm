package BioUtils::MyX::ConsensusBuilder;


use Exception::Class (
    'BioUtils::MyX::ConsensusBuilder' => {
    },
    
    'BioUtils::MyX::ConsensusBuilder::SeqsNotSqr' => {
        isa => 'BioUtils::MyX::ConsensusBuilder',
        fields => [ 'alignment_len', 'seq_len' ],
    },
    
    'BioUtils::MyX::ConsensusBuilder::QualsNotSqr' => {
        isa => 'BioUtils::MyX::ConsensusBuilder',
        fields => [ 'alignment_len', 'quals_len'],
    },
    
    'BioUtils::MyX::ConsensusBuilder::NoSeqs' => {
        isa => 'BioUtils::MyX::ConsensusBuilder',
    },
    
    'BioUtils::MyX::ConsensusBuilder::TooFewSeqs' => {
        isa => 'BioUtils::MyX::ConsensusBuilder',
    },
    
);

1;
__END__


#######
# POD #
#######
=head1 NAME

BioUtils::MyX::ConsensusBuilder - A hierarchy of exceptions that can be used
when working with ConsensusBuilder objects.

=head1 VERSION

This documentation refers to BioUtils::MyX::ConsensusBuilder version 1.2.0.

=head1 Included Modules

    Exception::Class

=head1 Inherit

    NA

=head1 SYNOPSIS

    # Throw a Fastq exception
    use BioUtils::MyX::ConsensusBuilder 1.2.0;
    if ( ... ) {   # Some code looking for an error
        BioUtils::MyX::ConsensusBuilder->throw(
            error => 'A ConsesusBuilder exception'
        );
    }
    
    # In caller catch the Generic exception
    eval { ... };
    if ( my $e = BioUtils::MyX::ConsensusBuilder->caught() ) {
        # Do something to handle the exception like print an error message
        print $e->error(), " via package ", $e->package(), " at ", $e->file,
            " line ", $e->line();
    }
    

=head1 DESCRIPTION

BioUtils::MyX::ConsensusBuilder holds a hierarchy of exception classes that are
associated with ConsensusBuilder objects.  

For more information what can be done when throwing and catching an exception
see Exception::Class and Exception::Class::Base.

=head1 CLASSES

=over

    BioUtils::MyX::ConsensusBuilder
    BioUtils::MyX::ConsensusBuilder::SeqsNotSqr
    BioUtils::MyX::ConsensusBuilder::QualsNotSqr
    BioUtils::MyX::ConsensusBuilder::NoSeqs
    BioUtils::MyX::ConsensusBuilder::TooFewSeqs
    
    
=back

=head1 CLASSES DESCRIPTION

=head2 BioUtils::MyX::ConsensusBuilder
    
    Title: BioUtils::MyX::ConsensusBuilder
    Throw Usage: BioUtils::MyX::ConsensusBuilder->throw(
                    error => 'Any ConsensusBuilder error message'
                );
    Catch Usage: if ( my $e = BioUtils::MyX::ConsensusBuilder->caught() ) { ... }
    Function: Throw/Catch a BioUtils::MyX::ConsensusBuilder exception
    Fields: error => an error message
    Inherits: NA
    Comments: NA
    See Also: NA

=head2 BioUtils::MyX::ConsensusBuilder::SeqsNotSqr

    Title: BioUtils::MyX::ConsensusBuilder::SeqsNotSqr
    Throw Usage: BioUtils::MyX::ConsensusBuilder::SeqsNotSqr->throw(
                    alignment_len => $alignment_len,
                    seq_len => $seq_len,
                 );
    Catch Usage: if ( my $e = BioUtils::MyX::ConsensusBuilder::SeqsNotSqr->caught() )
                    { ... }
    Function: Throw/Catch a BioUtils::MyX::ConsensusBuilder::SeqsNotSqr exception
              when the sequences in an alignment are not square
    Fields: alignment_len => The alignment length set by the length of the first
                             sequence.
            seq_len => The length of the sequence that doesn't match
    Inherits: BioUtils::MyX::ConsensusBuilder
    Comments: NA
    See Also: NA
    
=head2 BioUtils::MyX::ConsensusBuilder::QualsNotSqr

    Title: BioUtils::MyX::ConsensusBuilder::QualsNotSqr
    Throw Usage: BioUtils::MyX::ConsensusBuilder::QualsNotSqr->throw(
                    alignment_len => $alignment_len,
                    quals_len = $quals_len,
                 );
    Catch Usage: if ( my $e = BioUtils::MyX::ConsensusBuilder::QualsNotSqr->caught() )
                    { ... }
    Function: Throw/Catch a BioUtils::MyX::ConsensusBuilder::QualsNotSqr exception
              when the quality scores in an alignment are not square
    Fields: alignment_len => The alignment length set by the length of the first
                             sequence.
            quals_len => The length of the quality score string that doesn't match
    Inherits: BioUtils::MyX::ConsensusBuilder
    Comments: NA
    See Also: NA
    
=head2 BioUtils::MyX::ConsensusBuilder::NoSeqs

    Title: BioUtils::MyX::ConsensusBuilder::NoSeqs
    Throw Usage: BioUtils::MyX::ConsensusBuilder::NoSeqs->throw();
    Catch Usage: if ( my $e = BioUtils::MyX::ConsensusBuilder::NoSeqs->caught() )
                    { ... }
    Function: Throw/Catch a BioUtils::MyX::ConsensusBuilder::NoSeqs exception
              when there are no sequence object to use in building a consensus
    Fields: NA -- Only generic error fields defined by Exception::Class
    Inherits: BioUtils::MyX::ConsensusBuilder
    Comments: NA
    See Also: NA
    
=head2 BioUtils::MyX::ConsensusBuilder::TooFewSeqs

    Title: BioUtils::MyX::ConsensusBuilder::TooFewSeqs
    Throw Usage: BioUtils::MyX::ConsensusBuilder::TooFewSeqs->throw();
    Catch Usage: if ( my $e = BioUtils::MyX::ConsensusBuilder::TooFewSeqs->caught() )
                    { ... }
    Function: Throw/Catch a BioUtils::MyX::ConsensusBuilder::TooFewSeqs exception
              when all the sequences are different lengths.  This means that a
              consensus cannot be built without first building a MSA, which
              effectively makes all the seqs the same length.
    Fields: NA -- Only generic error fields defined by Exception::Class
    Inherits: BioUtils::MyX::ConsensusBuilder
    Comments: NA
    See Also: NA


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




