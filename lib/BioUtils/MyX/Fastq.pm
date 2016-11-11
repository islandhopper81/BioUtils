package BioUtils::MyX::Fastq;


use Exception::Class (
    'BioUtils::MyX::Fastq' => {
    },
    
    'BioUtils::MyX::Fastq::BadHeaderFormat' => {
        isa => 'BioUtils::MyX::Fastq',
        fields => [ 'header' ],
    },
    
);

1;
__END__


#######
# POD #
#######
=head1 NAME

BioUtils::MyX::Fastq - A hierarchy of exceptions that can be used when working
with fastq sequences and corresponding objects.

=head1 VERSION

This documentation refers to BioUtils::MyX::Fastq version 1.2.1.

=head1 Included Modules

    Exception::Class

=head1 Inherit

    NA

=head1 SYNOPSIS

    # Throw a Fastq exception
    use BioUtils::MyX::Fastq 1.2.1;
    if ( ... ) {   # Some code looking for an error
        BioUtils::MyX::Fastq->throw(
            error => 'A fastq exception'
        );
    }
    
    # In caller catch the Generic exception
    eval { ... };
    if ( my $e = BioUtils::MyX::Fastq->caught() ) {
        # Do something to handle the exception like print an error message
        print $e->error(), " via package ", $e->package(), " at ", $e->file,
            " line ", $e->line();
    }
    

=head1 DESCRIPTION

BioUtils::MyX::Fastq holds a hierarchy of exception classes that are associated
with fastq sequences and their corresponding classes.  

For more information what can be done when throwing and catching an exception
see Exception::Class and Exception::Class::Base.

=head1 CLASSES

=over

    BioUtils::MyX::Fastq
    BioUtils::MyX::Fastq::BadHeaderFormat
    
    
=back

=head1 CLASSES DESCRIPTION

=head2 BioUtils::MyX::Fastq
    
    Title: BioUtils::MyX::Fastq
    Throw Usage: BioUtils::MyX::Fastq->throw(
                    error => 'Any fastq error message'
                );
    Catch Usage: if ( my $e = BioUtils::MyX::Fastq->caught() ) { ... }
    Function: Throw/Catch a BioUtils::MyX::Fastq exception
    Fields: error => an error message
    Inherits: NA
    Comments: NA
    See Also: NA

=head2 BioUtils::MyX::Fastq::BadHeaderFormat

    Title: BioUtils::MyX::Fastq::BadHeaderFormat
    Throw Usage: BioUtils::MyX::Fastq::BadHeaderFormat->throw(
                    header => $header
                 );
    Catch Usage: if ( my $e = BioUtils::MyX::Fastq::BadHeaderFormat->caught() )
                    { ... }
    Function: Throw/Catch a BioUtils::MyX::Fastq::BadHeaderFormat exception when
              a header with an incorrect format is encounted.
    Fields: header => The header with incorrect format
    Inherits: BioUtils::MyX::Fastq
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




