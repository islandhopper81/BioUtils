#!/usr/bin/env perl

# converts a string of encoded quality values to ints

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile tempdir);
use Carp;
use Readonly;

use BioUtils::Codec::QualityScores qw(illumina_1_8_to_int);

# Subroutines #


# Variables #
my ($quals, $type, $help, $man);
Readonly my $TYPE => "illumina_1_8";

my $options_okay = GetOptions (
    "quals:s" => \$quals,
    "type:s" => \$type,
    "help" => \$help,                  # flag
    "man" => \$man,                     # flag (print full man page)
);

# check for input errors
if ( $help ) { pod2usage(2) }
if ( $man ) { pod2usage(-verbose => 2) }
if ( ! defined $quals ) { pod2usage(-msg => "No --quals value",
                                    -extval => 2); }
if ( ! defined $type ) { $type = $TYPE; }


### MAIN ###
# process each character in the sting
my $int_quals_aref = illumina_1_8_to_int($quals);
if ( ref $int_quals_aref eq 'ARRAY') {
    print(join(',', @{$int_quals_aref}), "\n");
}
else {
    # must be a scalar
    print $int_quals_aref, "\n";
}


########
# Subs #
########




__END__

# POD

=head1 NAME

quals_to_ints.pl - converts a string of encoded quality values to ints


=head1 VERSION

This documentation refers to quals_to_ints.pl version 0.0.1


=head1 SYNOPSIS

    quals_to_ints.pl
        --quals <encoded quality value string>
        --type [illumina_1_8]
        [--help]
        [--man]

    --quals  = A string of encoded quality values
    --type   = Type of quality value encoding
    --help  = Prints USAGE statement
    --man   = Prints the man page


=head1 ARGUMENTS
    
=head2 --quals

A string of encoded quality values.  REQUIRED
    
=head2 --type

The type of quality value encoding.  Right now the only supported option is
Illumina v1.8 (which is also the Sanger code).  DEFAULT: illumina_1_8

=head2 [--help]
    
An optional parameter to print a usage statement.

=head2 [--man]

An optional parameter to print he entire man page (i.e. all documentation)
    

=head1 DESCRIPTION

This perl script converts a string of encoded quality values into their integer
values.  Currently the only quality value encoding supported is the Illumina v1.8
which also happens to be the Sanger encoding.


=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Pod::Usage
File::Temp qw(tempfile tempdir)
Carp
Readonly
BioUtils::Codec::QualityScores qw(illumin_1_8_to_int);


=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
    
    
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2014, Scott Yourstone
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
