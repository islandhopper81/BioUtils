package BioUtils::Align::FastaSeq;


use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use MyX::Generic 0.0.1;
use version; our $VERSION = qv('1.2.0');

use BioUtils::FastaSeq;
use base qw(BioUtils::FastaSeq);

{
    Readonly my $NEW_USAGE => q{ new( { header => ,
                                        seq => ,
                                        start => ,
                                        end => ,
                                       } ) };

    # Attributes #
    my %start_of;
    my %end_of;
    
    # Setters #
    sub set_start;
    sub set_end;
    
    # Getters #
    sub get_start;
    sub get_end;
    
    
    # Others #
    sub _init;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
            
        # Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{start},
                               $arg_href->{end},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Bless a scalar to instantiate an object
        my $new_obj = $class->SUPER::new($arg_href);
        
        $new_obj->set_start($arg_href->{start});
        $new_obj->set_end($arg_href->{end});
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_start {
        my ($self, $start) = @_;
        $start_of{ident $self} = $start;
        return 1;
    }
    
    sub set_end {
        my ($self, $end) = @_;
        $end_of{ident $self} = $end;
        return 1;
    }
    
    
    
    ###########
    # Getters #
    ###########
    sub get_start {
        my ($self) = @_;
        return $start_of{ident $self};
    }
    
    sub get_end {
        my ($self) = @_;
        return $end_of{ident $self};
    }
    
    
    
    ##########
    # Others #
    ##########
}

1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::Align::FastaSeq - A data structure to holding an aligned fasta seq


=head1 VERSION

This document describes BioUtils::Align::FastaSeq version 1.2.0


=head1 SYNOPSIS

    use BioUtils::Align::FastaSeq
	
	my $aln_seq = BioUtils::Align::FastaSeq->new({
				header => "seq1",
				seq => "A-TCG",
				start => 10,
                end => 14
			});
	
	# Set new parameters
	$aln_seq->set_start(11);
    $aln_seq->set_end(15);
    
    # Get the data stored in the object
    # use the same getters as in BioUtils::FastaSeq and
    my $start = $aln_seq->get_start();
    my $end = $aln_seq->get_end();

  
=head1 DESRCIPTION

BioUtils::Align::FastaSeq is a data structure for storing a sequence that has
been aligned.  It inherits from BioUtils::FastaSeq.  The extra information
stored in these objects are the start and end positions of the aligned portion
of the sequence with repsect to the original sequence.


=head1 DEPENDENCIES

	version
	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	Scalar::Util qw(looks_like_number)
	MyX::Generic 1.2.0
    
    BioUtils::FastaSeq


=head1 INCOMPATIBILITIES

	None reported.

=head1 METHODS
	
	get_start
    get_end
    set_start
    set_end


=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: 	my $nw = BioUtils::Align::FastaSeq->new({
						header => "seq1",
						seq => "ATC-T",
						start => 10,
                        end => 14,
					});
	Function: Create a new BioUtils::Align::FastaSeq object
	Returns: Reference to a BioUtils::Align::FastaSeq object
	Args: -header => the sequence header
		  -seq => the sequence string
		  -start => the position in the original sequence where the
                        alignment starts
          -end => the position in the original sequence where the alignment ends
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA

=head2 set_start

	Title: set_start
	Usage: $aln_seq->set_start($start)
	Function: Sets the start of the alignment with respect to the original seq
	Returns: 1 on successful completion
	Args: start => a number for the start position
	Throws: NA
	Comments: 0 offset.  There can be gaps in the aligned sequence so the length
              of the sequecne does not necessarily correspond to the start and
              end positions.
	See Also: NA
    
=head2 set_end

	Title: set_end
	Usage: $aln_seq->set_end($end)
	Function: Sets the end of the alignment with respect to the original seq
	Returns: 1 on successful completion
	Args: end => a number for the end position
	Throws: NA
	Comments: 0 offset.  There can be gaps in the aligned sequence so the length
              of the sequecne does not necessarily correspond to the start and
              end positions.
	See Also: NA
	
=head2 get_start
	
	Title: get_start
	Usage: my $start = $aln_seq->get_start()
	Function: gets the start position of the alignment with respect to the
              original sequence.
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_end
	
	Title: get_end
	Usage: my $end = $aln_seq->get_end()
	Function: gets the end position of the alignment with respect to the
              original sequence.
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to author


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
