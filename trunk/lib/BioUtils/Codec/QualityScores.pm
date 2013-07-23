package BioUtils::Codec::QualityScores;

use warnings;
use strict;
use Readonly;
use Class::Std::Utils;
use MyX::Generic 1.0.5;
use Exporter qw( import );
our @EXPORT_OK = qw( int_to_illumina_1_8 illumina_1_8_to_int );
use version; our $VERSION = qv('1.0.5');


{
    # Global variables
	Readonly my $NEW_USAGE => q{ new() };
    Readonly my $MAX_INT => 93;
    Readonly my $MIN_INT => 0;
    Readonly my $ILL_1_8_BASE => 0;
    Readonly my $ILL_1_8_OFFSET => 33;
	
	# Class Attributes
	my %int_to_illumina_1_8_lookup_of;
	my %illumina_1_8_to_int_lookup_of;
	
	# Public Class Methods
	sub int_to_illumina_1_8_lookup;
	sub illumina_1_8_to_int_lookup;
	
	# Private Class Methods
	sub _set_int_to_illumina_1_8_lookup;
	sub _set_illumina_1_8_to_int_lookup;

    # Module Subroutines -- these can be exported w/o calling new
    sub int_to_illumina_1_8;
    sub illumina_1_8_to_int;
	

	# Constructor
	sub new {
		my ($class) = @_;
		
		# Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Initialize the objects attributes
        $int_to_illumina_1_8_lookup_of{ident $new_obj} =
			$new_obj->_set_int_to_illumina_1_8_lookup();
		$illumina_1_8_to_int_lookup_of{ident $new_obj} =
			$new_obj->_set_illumina_1_8_to_int_lookup();
		
        return $new_obj;
	}
	
	sub int_to_illumina_1_8_lookup {
		my ($class, $ints) = @_;
		
		my $illumina_str;
        my @ints_arr;
        
        # add to @ints_arr based on what was passed in (either int or aref).
        if ( ref $ints ) {
            @ints_arr = @$ints;
        }
        else {
            push @ints_arr, $ints;
        }
        
        # Convert int(s) to string of illumina scores
        foreach my $int ( @ints_arr ) {
            # make sure the int is in bounds
            if ( $int > $MAX_INT or
                 $int < $MIN_INT
                ) {
                MyX::Generic::OutOfBounds->throw(
					error => 'Int out of ASCII bounds'
				);
            }
            
            $illumina_str .= $int_to_illumina_1_8_lookup_of{ident $class}->{$int};
        }
        
        return $illumina_str;
	}
	
	sub illumina_1_8_to_int_lookup {
		my ($class, $chars) = @_;
		
		my @chars_arr = split //, $chars;
        my @ints_arr;  # the return values
        
        foreach my $char ( @chars_arr ) {
            my $int = $illumina_1_8_to_int_lookup_of{ident $class}->{$char};
            
            # Make sure the int is in bounds
            if ( ! defined $int or
				 $int > $MAX_INT or
                 $int < $MIN_INT
                ) {
                MyX::Generic::OutOfBounds->throw(
					error => 'Char out of ASCII bounds'
				);
            }
            
            push @ints_arr, $int;
        }
        
        if ( scalar @ints_arr == 1 ) {
            return $ints_arr[0]; # return as an int and not an arref
        }
        
        # Return an array ref of ints
        return \@ints_arr;
	}
	
	sub _set_int_to_illumina_1_8_lookup {
		my ($self) = @_;
		
		my %lookup = ();
		foreach my $int ( $MIN_INT..$MAX_INT ) {
			$lookup{ $int } = int_to_illumina_1_8($int);
		}
		
		return \%lookup;
	}
	
	sub _set_illumina_1_8_to_int_lookup {
		my ($self) = @_;
		
		my %lookup = ();
		foreach my $int ( $MIN_INT..$MAX_INT ) {
			$lookup{ int_to_illumina_1_8($int) } = $int;
		}
		
		return \%lookup;
	}
    

    sub int_to_illumina_1_8 {
        my ($ints) = @_;
        
        my $illumina_str;
        my @ints_arr;
        
        # add to @ints_arr based on what was passed in (either int or aref).
        if ( ref $ints ) {
            @ints_arr = @$ints;
        }
        else {
            push @ints_arr, $ints;
        }
        
        # Convert int(s) to string of illumina scores
        foreach my $int ( @ints_arr ) {
            # make sure the int is in bounds
            if ( $int > $MAX_INT or
                 $int < $MIN_INT
                ) {
                MyX::Generic::OutOfBounds->throw(
					error => 'Int out of ASCII bounds'
				);
            }
            
            $illumina_str .= chr($int + $ILL_1_8_OFFSET);
        }
        
        return $illumina_str;
    }
    
    sub illumina_1_8_to_int {
        my ($chars) = @_;
        
        my @chars_arr = split //, $chars;
        my @ints_arr;  # the return values
        
        foreach my $char ( @chars_arr ) {
            my $int = (ord $char) - $ILL_1_8_OFFSET;
            
            # Make sure the int is in bounds
            if ( $int > $MAX_INT or
                 $int < $MIN_INT
                ) {
                MyX::Generic::OutOfBounds->throw(
					error => 'Char out of ASCII bounds'
				);
            }
            
            push @ints_arr, $int;
        }
        
        if ( scalar @ints_arr == 1 ) {
            return $ints_arr[0]; # return as an int and not an arref
        }
        
        # Return an array ref of ints
        return \@ints_arr;
    }
}


1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::Codec::QualityScores - Convert sequences quality score integers to
encoded formats


=head1 VERSION

This document describes BioUtils::Codec::QualityScores version 1.0.5


=head1 SYNOPSIS

    use BioUtils::Codec::QualityScores;
    use BioUtils::Codec::QualityScores
		qw(int_to_illumina_1_8 illumina_1_8_to_int);
    
    # Convert an int value to illumina-1.8 quality encoding
    eval {
        my $encoding = int_to_illumina_1_8(40);
    };
    if ( my $e = MyX::Generic::OutOfBounds->caught() ) {
        # handle exception
    }
    
    # Convert a single illumina-1.8 quality encoding to an int
    eval {
        my $int_score = illumina_1_8_to_int('!');
    };
    if ( my $e = MyX::Generic::UnrecignizedChar->caught() ) {
        # handle exception
    }
    
    # Convert an array reference of int values to a string of illumina-1.8
    # quality encodings
    my $int_aref = [40,40,40,40];
    eval {
        my $encoding_str = int_to_illumina_1_8($int_aref);
    }
    if ( my $e = MyX::Generic::OutOfBounds->caught() ) {
        # handle exception
    }
    
    # Convert a string of illumina-1.8 quality encodings to an array reference
    # of ints
    my $encoding_str = "IIII";
    eval {
        my $int_scores_aref = illumina_1_8_to_int($encoding_str);
    };
    if ( my $e = MyX::Generic::UnrecignizedChar->caught() ) {
        # handle exception
    }
    
  
=head1 DESCRIPTION

This perl library can be used either as a module where subroutines are exported
or as an object where an obect can be made and methods called on that object.
The exportable subroutines include int_to_illumina_1_8 and illumina_1_8_to_int.

This library has a collection of subroutines to convert between integers and
encoded quality scores.  There are several different encoding that have been
used to encode quality scores.  Right now this module only converts between
integers and the most recent illumina quality encoding (illumina-1.8).  This
encoding is synonymous with the original Sanger encoding called Phred+33.

I experimented with adding lookup methods where a table of values is stored in a
hash.  I thought this would be faster than converting each value on the fly.
After benchmarking both methods it turns out they are about equal.  The lookup
might be slightly faster, but not enough to make much of a difference.


=head1 METHODS 

	int_to_illumina_1_8_lookup
	illumina_1_8_to_int_lookup
	_set_int_to_illumina_1_8_lookup
	_set_illumina_1_8_to_int_lookup
	int_to_illumina_1_8
	illumina_1_8_to_int


=head1 CONFIGURATION AND ENVIRONMENT
  
BioUtils::Codec::QualityScores requires no configuration files or environment
variables.


=head1 DEPENDENCIES

	Readonly
	Class::Std::Utils
	MyX::Generic
	Exporter qw( import )
	version


=head1 INCOMPATIBILITIES

	None reported.


=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: my $codec = BioUtils::Codec::QualityScores->new();
    Function: Creates a new Codec::QualityScores object
    Returns: Codec::QualityScores
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 int_to_illumina_1_8_lookup
    
    Title: int_to_illumina_1_8_lookup
    Usage: $codec->int_to_illumina_1_8_lookup(40);
    Function: Looks up the illumina_1.8 quality encoding for the given int
    Returns: String
    Args: -int => the int to look up
    Throws: MyX::Generic::OutOfBounds
    Comments: Must be between 0 and 93
    See Also: NA
	
=head2 illumina_1_8_to_int_lookup
    
    Title: illumina_1_8_to_int_lookup
    Usage: $codec->illumina_1_8_to_int_lookup('K');
    Function: Looks up the int value for the given illumina_1.8 quality encoding
    Returns: Int
    Args: -qual => the encoded quality value to look up
    Throws: MyX::Generic::UnrecignizedChar
    Comments: Must be between ! and ~
    See Also: NA

=head2 _set_int_to_illumina_1_8_lookup
    
    Title: _set_int_to_illumina_1_8_lookup
    Usage: $codec->_set_int_to_illumina_1_8_lookup();
    Function: Create the lookup table for int to illumina_1_8
    Returns: Href
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 _set_illumina_1_8_to_int_lookup
    
    Title: _set_illumina_1_8_to_int_lookup
    Usage: $codec->_set_illumina_1_8_to_int_lookup();
    Function: Create the lookup table for illumina_1_8 to int
    Returns: Href
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 int_to_illumina_1_8

	Title: int_to_illumina_1_8
	Usage: int_to_illumina_1_8($int) OR int_to_illumina_1_8($int_aref);
	Function: Convert an int to an encoded illumina-1.8 quality score
	Returns: Str
	Args: -int => a single int value
          -int_aref => an array reference of int values
	Throws: MyX::Generic::OutOfBounds
	Comments: NA
	See Also: NA

=head2 illumina_1_8_to_int

	Title: illumina_1_8_to_int
	Usage: illumina_1_8_to_int($char);
	Function: Convert a character to an encoded illumina-1.8 quality score
	Returns: int or array reference of ints
	Args: -char => a string (any length) of encoded illumina-1.8 quality scores
	Throws: MyX::Generic::UnrecignizedChar
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
