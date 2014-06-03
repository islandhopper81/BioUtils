package BioUtils::Codec::IUPAC;

use warnings;
use strict;
use Readonly;
use Class::Std::Utils;
use MyX::Generic 1.0.9;
use Exporter qw( import );
our @EXPORT_OK = qw( nuc_str_to_iupac iupac_to_nuc_str );
use version; our $VERSION = qv('1.0.9');


{
    # Global variables

    # Module Subroutines -- these can be exported w/o calling new
	sub nuc_str_to_iupac;
	sub iupac_to_nuc_str;
	
	# Private Subroutines
	sub _get_sorted;
	
	# Subroutines #
	sub nuc_str_to_iupac {
		my ($str) = @_;
		
		if ( ! defined $str ) {
			MyX::Generic::Undef::Param->throw(
				usage => 'nuc_str_to_iupac(str)'
			);
			
			return "";
		}
        
        my $sorted = _get_sorted($str);
        
        if ( $sorted eq "A" ) { return "A"; }
        elsif ( $sorted eq "C" ) { return "C"; }
        elsif ( $sorted eq "G" ) { return "G"; }
        elsif ( $sorted eq "T" ) { return "T"; }
        elsif ( $sorted eq "AG" ) { return "R"; }
        elsif ( $sorted eq "CT" ) { return "Y"; }
        elsif ( $sorted eq "CG" ) { return "S"; }
        elsif ( $sorted eq "AT" ) { return "W"; }
        elsif ( $sorted eq "GT" ) { return "K"; }
        elsif ( $sorted eq "AC" ) { return "M"; }
        elsif ( $sorted eq "CGT" ) { return "B"; }
        elsif ( $sorted eq "AGT" ) { return "D"; }
        elsif ( $sorted eq "ACT" ) { return "H"; }
        elsif ( $sorted eq "ACG" ) { return "V"; }
        elsif ( $sorted eq "ACGT" ) { return "N"; }
		elsif ( $sorted =~ m/N/ ) { return "N"; }
		else {
			# the parameter is not a nucleotide base
			MyX::Generic::BadValue->throw(
				error => "Unrecognized nucleotide base",
				value => $str,
			);
			
			return '*';
		}
	}
	
	sub iupac_to_nuc_str {
		my ($iupac) = @_;
		
		if ( ! defined $iupac ) {
			MyX::Generic::Undef::Param->throw(
				usage => 'iupac_to_nuc_str(iupac)'
			);
			
			return "";
		}
		
		if ( $iupac eq "A" ) { return "A"; }
        elsif ( $iupac eq "C" ) { return "C"; }
        elsif ( $iupac eq "G" ) { return "G"; }
        elsif ( $iupac eq "T" ) { return "T"; }
        elsif ( $iupac eq "R" ) { return "AG"; }
        elsif ( $iupac eq "Y" ) { return "CT"; }
        elsif ( $iupac eq "S" ) { return "CG"; }
        elsif ( $iupac eq "W" ) { return "AT"; }
        elsif ( $iupac eq "K" ) { return "GT"; }
        elsif ( $iupac eq "M" ) { return "AC"; }
        elsif ( $iupac eq "B" ) { return "CGT"; }
        elsif ( $iupac eq "D" ) { return "AGT"; }
        elsif ( $iupac eq "H" ) { return "ACT"; }
        elsif ( $iupac eq "V" ) { return "ACG"; }
		elsif ( $iupac eq "N" ) { return "ACGTN"; }
        else {
			# the parameter is not an IUPAC coded letter
			MyX::Generic::BadValue->throw(
				error => "Not an IUPAC coded character",
				value => $iupac,
			);
			
			
			return 0;
		}
	}
	
	sub _get_sorted {
		my ($str) = @_;
		
        $str =~ tr/[a-z]/[A-Z]/;
        
        my @bases = split(//, $str);
        my @sorted = sort {uc($a) cmp uc($b)} @bases;
        
        return (join('', @sorted));
	}
}


1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::Codec::IUPAC - Converts between residues and IUPAC codes


=head1 VERSION

This document describes BioUtils::Codec::IUPAC version 1.0.9


=head1 SYNOPSIS

    use BioUtils::Codec::IUPAC qw(nuc_str_to_iupac iupac_to_nuc_str);
    
    # Convert an nucleotide string to IUPAC
    eval {
        my $iupac = nuc_str_to_iupac('AT');
    };
    if ( my $e = MyX::Generic::BadValue->caught() ) {
        # handle exception
    }
    
    # Convert a IUPAC value to nucleotide string
    eval {
        my $nuc_str = iupac_to_nuc_str('W');
    };
    if ( my $e = MyX::Generic::BadValue->caught() ) {
        # handle exception
    }
    
  
=head1 DESCRIPTION

This perl library is used as a module where subroutines are exported. The
exportable subroutines include nuc_str_to_iupac and iupac_to_nuc_str.

This library has a collection of subroutines to convert between residues and
IUPAC values. 


=head1 METHODS 

	nuc_str_to_iupac
	iupac_to_nuc_str
	_get_sorted

=head1 CONFIGURATION AND ENVIRONMENT
  
BioUtils::Codec::IUPAC requires no configuration files or environment
variables.


=head1 DEPENDENCIES

	Readonly
	Class::Std::Utils
	MyX::Generic 1.0.9
	Exporter qw( import )
	version
	use Readonly;


=head1 INCOMPATIBILITIES

	None reported.


=head1 METHODS DESCRIPTION

=head2 nuc_str_to_iupac
    
    Title: nuc_str_to_iupac
    Usage: my $iupac = nuc_str_to_iupac($nuc_str);
    Function: Converts a string of nucleotide bases to an IUPAC value
    Returns: String
    Args: -nuc_str => a string of concatenated nucleotides.  All nucleotides
					  should be unique.  (e.g. "AT")
    Throws: MyX::Generic::BadValue. MyX::Generic::Undef::Param
    Comments: All nucleotides should be unique.  'ATA' is not an acceptable
			  string
    See Also: NA
	
=head2 iupac_to_nuc_str
    
    Title: iupac_to_nuc_str
    Usage: my $nuc_str = iupac_to_nuc_str($iupac);
    Function: Converts an IUPAC value to a string of nucleotides
    Returns: String
    Args: -iupac => a single IUPAC value (e.g. "W")
    Throws: MyX::Generic::BadValue, MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA

=head2 _get_sorted
    
    Title: _get_sorted
    Usage: _get_sorted($string);
    Function: Sorts the string of nucloetide bases so I can look them up to
			  return he appropriate IUPAC value
    Returns: String
    Args: $string => A string of unique nucloetides
	Throws: NA
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
