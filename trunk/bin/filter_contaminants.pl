#! /usr/bin/env perl

# This is the driver script to filter contaminants from a seq and/or otu table 
#
# Scott Yourstone
# scott.yourstone81@gmail.com

use strict;
use warnings;

use BioUtils::QC::ContaminantFilter 1.0.6;
use Getopt::Long;
use Config::Std;

print "Start\n\n";

# USAGE
my $usage = "$0 --params_file file [--help]\n";

# Variables to be set by GetOpts and their defaults
my $params_file = undef;
my $help = 0;

my $options_okay = GetOptions (
    "params_file:s" => \$params_file,  # string
    "help" => \$help,                  # flag
);

# Fail if unknown options are encountered
if ( ! $options_okay ) {
    die $usage;
}

# Print a help message
if ( $help ) {
    die $usage;
}

# Fail if no params file is given
if ( ! defined $params_file ) {
    die $usage;
}

print "Read Config File\n";
read_config( $params_file => my %hash);

print "Building Filter\n\n";
my $filter = BioUtils::QC::ContaminantFilter->new({
                blast_db => $hash{'BLAST Parameters'}{'db'},
                query_file => $hash{'BLAST Parameters'}{'query'},
                output_dir => $hash{'Output'}{'dir'},
                output_prefix => $hash{'Output'}{'prefix'},
                eval => $hash{'BLAST Parameters'}{'evalue'},
                perc_iden => $hash{'BLAST Parameters'}{'perc_identity'},
                output_fmt => $hash{'BLAST Parameters'}{'output_fmt'},
                max_targets => $hash{'BLAST Parameters'}{'max_targets'},
                otu_table => $hash{'OTU Table'}{'otu_table_file'},
                });

print "Running Filter\n\n";
$filter->run_filter();

print "Finished\n\n";


__END__

# POD

=head1 NAME

filter_contaminants.pl - Filters out contaminant sequences using blastn


=head1 VERSION

This documentation refers to filter_contaminants.pl version 0.0.2


=head1 USAGE

perl filter_contaminants.pl --params_file <file>


=head1 REQUIRED ARGUMENTS
    
    --params_file
        A file path to a file which contains all necessary and optional parameters


=head1 DESCRIPTION

In sequencing projects filtering out contaminant sequences is almost always a
necessary step.  Contaminant sequences are NOT low quality sequences.
Contaminant sequences come from organisms or samples that should be ignored.
For example, in plant metagenome 16S profiling sequence projects we would like
to ignore any 16S reads that match very closely to the host plant 16S sequence.

To identify contaminants a set of sequences is blasted against a known set of
contaminant sequences.

In 16S profiling it is typical to build OTU sequences and an OTU table before
filtering contaminants.  An optional OTU table can be supplied to this script
from which contaminant OTUs will be seperated from non-contaminant OTUs.


=head1 INPUTS

    The only required parameter is the file path to the paramters file.  The
    parameters file is described below.
    
    
=head1 PARAMETERS FILE
    
    The parameters file is formated as described in Damian Conway's perl package
    
    L<http://search.cpan.org/~bricker/Config-Std-0.900/lib/Config/Std.pm>
    
    There are three important headers in the params file:
    
    1. BLAST Paramters
    2. OTU Table
    3. Output
    
    The parameters for each section are included below their respective headers.
    An example of a correctly formated parameters file is below:
    
    [BLAST Parameters]
    db: /home/me/blast_db
    query: /home/me/my_seqs_file.fasta
    evalue: 0.00001
    perc_identity: 80
    
    [OTU Table]
    otu_table_file: /home/me/my_otu_table_file.txt
    
    [Output]
    dir: /home/me/my_output/
    prefix: my_results


=head1 OUTPUTS
    
    Output file are printed to the output directory specified by the user in the
    parameters file. There are two main output files:  1) a fasta file
    containing all the contaminant sequences and 2) a fasta file containing
    all the non-contaminant sequences.
    
    If an OTU table is provided to filter out OTU contaminants then there are 2
    more output files: 1) a tab delimated OTU table text file with contaminant
    OTUs and 2) a tab delimated OTU table text file with non-contaminant OTUs.
    

=head1 CONFIGURATION AND ENVIRONMENT
    
    BLAST 2.2.25+
    Tested using perl5.8.9 and perl5.12.3.
    
    
    
=head1 DEPENDANCIES

    BioUtils::QC::ContaminantFilter 1.0.6
    Getopt::Long
    Config::Std


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to the author


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

