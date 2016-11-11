#! /usr/bin/env perl

# converts all version numbers to a new version number

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Carp;
use BioUtils::FastaIO;

# Subroutines #
sub my_reverse;
sub my_complement;
sub _operations;

# Variables #
my ($str, $file, $rev_flag, $comp_flag, $help);

my $options_okay = GetOptions (
    "str:s" => \$str,
	"file:s" => \$file,
    "rev_flag" => \$rev_flag,
    "comp_flag" => \$comp_flag,
    "help" => \$help,                  # flag
);

# check for input errors
if ( $help ) { pod2usage(2) }
if ( ! defined $str and ! defined $file ) { pod2usage(2) }

# if none of the flags are set automatically do a reverese and compliment
if ( ! defined $rev_flag and ! defined $comp_flag ) {
	$rev_flag = 1;
	$comp_flag = 1;
}

# MAIN
if ( defined $file and -s $file ) {
	# create the bioUtils fasta file object
	my $in = BioUtils::FastaIO->new({stream_type => '<', file => $file});
	while ( my $seq = $in->get_next_seq() ) {
		my $new_seq = _operations($seq->get_seq());
		$seq->set_seq($new_seq);
		print ">", $seq->get_header(), "\n";
		print $seq->get_seq(), "\n";
	}
}
else {
	print _operations($str), "\n";
}






########
# Subs #
########
sub _operations {
	my ($str) = @_;
	
	my $new_str;
	
	if ( $rev_flag and $comp_flag ) {
		# both reverse and complement
		$new_str = my_complement(my_reverse($str));
	}
	elsif ( $rev_flag ) {
		$new_str = my_reverse($str);
	}
	elsif ( $comp_flag ) {
		$new_str = my_complement($str);
	}
	else {
		warn "ERROR";
	}
	
	return $new_str;
}

sub my_reverse {
    my ($str) = @_;
    my $ans = scalar reverse $str;
	$ans =~ tr/[]/][/;
    return $ans;
}

sub my_complement {
    my ($str) = @_;
    
    $str =~ tr/ACGTacgt/TGCAtgca/;
    return $str;
}



__END__

# POD

=head1 NAME

reverse_complement.pl - Simple reverse complement for DNA


=head1 VERSION

This documentation refers to reverse_complement.pl version 0.0.1


=head1 SYNOPSIS

    reverse_complement.pl
					--str ATGTA
					[--file <file.fasta>]
                    [--rev_flag]
                    [--comp_flag]
                    [--help]

    --str  = DNA string to reverse complement
	--file = A fasta file with DNA sequences
    --rev_flag   = ONLY reverse the string
    --comp_flag  = ONLY complement the string
    --help  = Prints USAGE statement


=head1 ARGUMENTS
    
=head2 --str

DNA string to reverse complement

=head2 [--file <file.fasta>]

Alternitively to inputing a DNA sequence directly this options allows users to
pass a fasta file.  The resulting sequecnes will be printed to the screen.
    
=head2 [--rev_flag]

ONLY reverse the string
    
=head2 [--comp_flag]

ONLY complement the string

=head2 [--help]
    
An optional parameter to print a usage statement.
    

=head1 DESCRIPTION

Reverse complements a DNA string


=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

Getopt::Long
Pod::Usage
Carp
BioUtils::FASTAIO


=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
    
    
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
