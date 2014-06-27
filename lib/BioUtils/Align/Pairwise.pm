package BioUtils::Align::Pairwise;

use warnings;
use strict;

use version; our $VERSION = qv('1.0.11');
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use Scalar::Util qw(looks_like_number);

use MyX::Generic 1.0.11;


{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new( {
										match_score => ,
										mismatch_score => ,
										gap_score => ,
										} ) };
	Readonly my $MATCH_SCORE => 5;
	Readonly my $MISMATCH_SCORE => -1;
	Readonly my $GAP_SCORE => -2;
	
    # Attributes #
	my %match_score_of;
	my %mismatch_score_of;
	my %gap_score_of;
    
    # Subroutines #
	sub _init;
	sub get_match_score;
	sub get_mismatch_score;
	sub get_gap_score;
	
	sub set_match_score;
	sub set_mismatch_score;
	sub set_gap_score;
    
    
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
		
		# Initialize the parameters
		$new_obj->_init($arg_href);
        
        return $new_obj;
    }
    
    ###############
	# Subroutines #
	###############
	sub _init {
		my ($self, $arg_href) = @_;
		
		$self->set_match_score($arg_href->{match_score});
		$self->set_mismatch_score($arg_href->{mismatch_score});
		$self->set_gap_score($arg_href->{gap_score});
		
		return 1;
	}
	
	sub set_match_score {
		my ($self, $match_score) = @_;
		
		if ( ! defined $match_score ) {
			$match_score_of{ident $self} = $MATCH_SCORE;
			return 1;
		}
		
		if ( ! looks_like_number($match_score) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Match Score must be a digit",
				value => $match_score,
			);
		}
		
		$match_score_of{ident $self} = $match_score;
		
		return 1;
	}
	
	sub set_mismatch_score {
		my ($self, $mismatch_score) = @_;
		
		if ( ! defined $mismatch_score ) {
			$mismatch_score_of{ident $self} = $MISMATCH_SCORE;
			return 1;
		}
		
		if ( ! looks_like_number($mismatch_score) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Mismatch Score must be a digit",
				value => $mismatch_score,
			);
		}
		
		$mismatch_score_of{ident $self} = $mismatch_score;
		
		return 1;
	}
	
	sub set_gap_score {
		my ($self, $gap_score) = @_;
		
		if ( ! defined $gap_score ) {
			$gap_score_of{ident $self} = $GAP_SCORE;
			return 1;
		}
		
		if ( ! looks_like_number($gap_score) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Gap Score must be a digit",
				value => $gap_score,
			);
		}
		
		$gap_score_of{ident $self} = $gap_score;
		
		return 1;
	}
	
	sub get_match_score {
		my ($self) = @_;
		
		return $match_score_of{ident $self};
	}
	
	sub get_mismatch_score {
		my ($self) = @_;
		
		return $mismatch_score_of{ident $self};
	}
	
	sub get_gap_score {
		my ($self) = @_;
		
		return $gap_score_of{ident $self};
	}

}
	
	

1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::Align::Pairwise - An abstract parent module for SW and NW


=head1 VERSION

This document describes BioUtils::Align::Pairwise version 1.0.11


=head1 SYNOPSIS

    See documentation for SW and NW child classes.

  
=head1 DESRCIPTION

An abstract parent module for SW and NW 


=head1 DEPENDENCIES

	version
	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	Scalar::Util qw(looks_like_number)

	MyX::Generic 1.0.11


=head1 INCOMPATIBILITIES

	None reported.

=head1 METHODS
	
	_init
	get_match_score
	get_mismatch_score
	get_gap_score
	
	set_match_score
	set_mismatch_score
	set_gap_score


=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: 	See SW and NW child classes
	Function: See SW and NW child classes
	Returns: See SW and NW child classes
	Args: See SW and NW child classes
	Throws: See SW and NW child classes
	Comments: Because this module is considered abstract, new should not be
			  called on it explicitely.  Call new via the SW and NW child
			  classes.
	See Also: BioUtils::Align::Pairwise::SW
			  BioUtils::Align::Pairwise::NW

=head2 _init

	Title: _init
	Usage: $nw->_init($arg_href)
	Function: Initializes the parameters
	Returns: 1 on completion
	Args: arg_href => arguments hash reference
	Throws: MyX::Generic::Digit::MustBeDigit
	Comments: NA
	See Also: NA

=head2 set_match_score

	Title: set_match_score
	Usage: $obj->set_match_score($match_score)
	Function: Sets the match score parameter
	Returns: 1 on successful completion
	Args: match_score => a number for the match score parameter
	Throws: MyX::Generic::Digit::MustBeDigit
	Comments: NA
	See Also: NA
	
=head2 set_mismatch_score

	Title: set_mismatch_score
	Usage: $obj->set_mismatch_score($mismatch_score)
	Function: Sets the mismatch score parameter
	Returns: 1 on successful completion
	Args: mismatch_score => a number for the mismatch score parameter
	Throws: MyX::Generic::Digit::MustBeDigit
	Comments: NA
	See Also: NA
	
=head2 set_gap_score

	Title: set_gap_score
	Usage: $obj->set_gap_score($gap_score)
	Function: Sets the gap score parameter
	Returns: 1 on successful completion
	Args: gap_score => a number for the gap score parameter
	Throws: MyX::Generic::Digit::MustBeDigit
	Comments: NA
	See Also: NA
	
=head2 get_match_score
	
	Title: get_match_score
	Usage: my $match_score = $obj->get_match_score()
	Function: gets the current match score parameter
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 get_mismatch_score
	
	Title: get_mismatch_score
	Usage: my $mismatch_score = $obj->get_mismatch_score()
	Function: gets the current mismatch score parameter
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_gap_score
	
	Title: get_gap_score
	Usage: my $gap_score = $obj->get_gap_score()
	Function: gets the current gap score parameter
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
