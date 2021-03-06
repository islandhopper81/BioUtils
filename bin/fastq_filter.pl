#! /usr/bin/env perl

# A perl script to filter fastq sequences

use strict;
use warnings;

use version; our $VERSION = qv('1.2.1');

use BioUtils::QC::FastqFilter 1.2.1;
use Getopt::Long;
use Carp qw(cluck);
use Pod::Usage;


# Variables
my $input_file = undef;
my $input_file2 = undef;
my $output_dir = undef;
my $min_length = undef;  # Default is keep all short sequences
my $min_average_qual = 0;  # Defualt is keep all with low average quality
my $min_base_qual = 0;  # Default is keep all with low base qualities
my $min_c_score = 0;
my $allow_gaps = 0;  # Default is don't allow gaps
my $allow_ambiguous_bases = 0;  # Default is don't allow ambiguous bases

my $help = 0;
my $verbose = 0;

my $options_okay = GetOptions (
    "input_file:s" => \$input_file,
    "input_file2:s" => \$input_file2,
    "output_dir:s" => \$output_dir,
    "min_length:i" => \$min_length,
    "min_average_qual:i" => \$min_average_qual,
    "min_base_qual:i" => \$min_base_qual,
    "min_c_score:i" => \$min_c_score,
    "allow_gaps" => \$allow_gaps,
    "allow_ambiguous_bases" => \$allow_ambiguous_bases,
    
    "help" => \$help,                  # flag
    "verbose" => \$verbose,            # flag
);

# Check the parameters
_check_params();

# Create a BioUtils::QC::FastqFilter object
my $filter_obj = BioUtils::QC::FastqFilter->new({
                        fastq_file => $input_file,
                        fastq_file2 => $input_file2,
                        output_dir => $output_dir,
                        min_len => $min_length,
                        min_avg_qual => $min_average_qual,
                        min_base_qual => $min_base_qual,
                        min_c_score => $min_c_score,
                        allow_gaps => $allow_gaps,
                        allow_ambig_bases => $allow_ambiguous_bases,
                        verbose => $verbose,
                    });

if ( defined $input_file2 ) {
    # run the filter on pairs
    $filter_obj->filter_pairs();
}
else {
    # run the normal filter
    $filter_obj->filter();
}

sub _check_params {
    # Fail if unknown options are encountered
    if ( ! $options_okay ) {
        pod2usage(2);
    }
    
    # Print a help message
    if ( $help ) {
        pod2usage(1);
    }
    
    # input_file is required parameter
    if ( ! $input_file ) {
        pod2usage(-msg => "Undefined --input_file\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # output_dir is a required parameter
    if ( ! $output_dir ) {
        pod2usage(-msg => "Undefined --output_dir\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # Make sure the input file exists
    if ( ! -f $input_file ) {
        pod2usage(-msg => "Input file doesn't exist\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # Make sure the output directory exists
    if ( ! -d $output_dir ) {
        pod2usage(-msg => "Output directory doesn't exist\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # Remember $min_length might be undef which indicates that we keep the
    # entire sequence.
    if ( defined $min_length and $min_length < 0 ) {
        pod2usage(-msg => "--min_length must be a positive integer\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # min_average_qual must be an int greater than 0
    if ( $min_average_qual < 0 ) {
        pod2usage(-msg => "--min_average_qual must be a positive integer\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # min_base_qual must be an int greater than 0
    if ( $min_base_qual < 0 ) {
        pod2usage(-msg => "--min_base_qual must be a positive integer\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # min_c_score must be an int greater than 0
    if ( $min_c_score < 0 ) {
        pod2usage(-msg => "--min_c_score must be a positive digit\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # allow_gaps must be either 0 (NO) or 1 (YES)
    if ( ! ($allow_gaps == 0 or $allow_gaps == 1) ) {
        pos2usage(-msg => "--allow_gaps must be 0 (NO) or 1 (YES)\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
    
    # allow_ambiguous_bases must be either 0 (NO) or 1 (YES)
    if ( ! ( $allow_ambiguous_bases == 0 or $allow_ambiguous_bases == 1) ) {
        pos2usage(-msg => "--allow_ambiguous_bases must be 0 (NO) or 1 (YES)\n",
                  -exitval => 2,
                  -verbose => 0,
                  );
    }
}

__END__

# POD

=head1 NAME

fastq_filter - Filter FASTQ sequences


=head1 VERSION

This documentation refers to fastq_filter.pl version 1.2.1


=head1 SYNOPSIS

    fastq_filter --output_dir --input_file [--input_file2] --min_length
                 --min_average_qual --min_base_qual --allow_gaps
                 --allwo_ambiguous_bases [--verbose] [--help]

    --input_file            = Path to the FASTQ input file
    --input_file2           = Path to the FASTQ input file of matched pairs
    --output_dir            = Path output directory
    --min_length            = A number for the minimum sequence length.
                              DEFAULT = none
    --min_average_qual      = Minimum average quality score for an entire read
                              DEFAULT = 0
    --min_base_qual         = Minimum single base quality score
                              DEFAULT = 0
    --min_c_score           = Minimum c-score value if given in header
                              DEFUALT = 0
    --allow_gaps            = Boolean signalling to retain sequences with gaps
    --allow_ambiguous_bases = Boolean signalling to retain IUPAC coded bases
    --verbose               = Runs all test for each sequence
    --help                  = Prints USAGE statement


=head1 ARGUMENTS
    
=head2 --input_file

The path to the FASTQ formated input file.
    
=head2 [--input_file2]

The patht to the FASTQ formated input file of matched pair.  OPTIONAL.
    
=head2 --output_dir

The path to the output directory.

=head2 --min_length

The minimum sequence length of sequences to keep.
    
=head2 --min_average_qual

The minimum average quality value for the enitre read given as an integer.
    
=head2 --min_base_qual
    
The minimum base quality value.  Any read that has at least one base below
the --min_base_qual value will be filtered out of the final set.  IMPORTANT:
It is important to not that most gaps are scored as a 0.  Therefore you may
set --min_base_qual to some number greater than 0 and that will effectively
remove the gaps even if --allow_gaps is set.

=head2 --min_c_score

A c-score is a metric for a consensus sequence.  It signifies the similarity of
the sequences that were used to create the consensus sequence.  If a sequences
has a c-score it should be declared in the header as c_score=40.0.
    
=head2 --allow_gaps

A flag to keep sequences that have gap characters in them.  If this flag is
not set on the commandline sequences with gap characters are automatically
filtered out of the final set.  IMPORTANT: if --min_base_qual is greater
than 0 the sequences with gaps (which generally have a quality score of 0)
are removed even if --allow_gaps is set.
    
=head2 --allow_ambiguous_bases
    
A flag to keep sequences that have ambiguous IUPAC bases.  if this flag is
not set on the commandline sequences with ambiguous bases are automatically
filtered out fo the final set.
    
=head2 [--verbose]
    
An optional flag too allow for extra processing and extra output files.  One
may want to know why the sequences are being filtered out or what parameters
are causing most sequences to be filtered out.  By turing on this flag each
sequence is tested by each parameter.  A table is output that can be loaded
into gnuplot or other software to diagnose what is filtering out most reads.
Also a file with all the filtered sequences is output for downstream analysis
as needed.  

=head2 [--help]
    
An optional parameter to print a usage statement.
    

=head1 DESCRIPTION

This perl script filters out sequences as specified by the user input parameters.

=head2 OUTPUTS

If the --verbose flag is not used there will be one output file containing
all the sequneces that pass ALL the tests specified by the parameters.
These sequences are called 'unfiltered' sequences.  This output file will be
in FASTQ format, and will be called unfiltered.fastq.

If the --verbose flag is included in the commmand there will be three output    
files.  The first is the typical FASTQ file contain sequences that pass ALL
the tests specified by the parameters.  The second is a FASTQ file that
contains all sequences that do NOT pass at least ONE test (all remaining
sequences that are not in the first file).  This file will be called
filtered.fastq.  The third is a table with a line for each sequence that
failed a test showing which tests it failed.  This table can be loaded into
a gnuplot scipt or other script to show why what tests most of the sequenes
failed.  Parameter can be adjusted accordingly.


=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Carp qw(cluck)
Pod::Usage
BioUtils::QC::FastqFilter 1.2.1


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
