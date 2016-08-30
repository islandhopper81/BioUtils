#! /usr/bin/env perl

# This is the driver script to filter contaminants from a seq and/or otu table 
#
# Scott Yourstone
# scott.yourstone81@gmail.com

use strict;
use warnings;

use version; our $VERSION = qv('1.2.0');

use BioUtils::QC::ContaminantFilter 1.2.0;
use BioUtils::FastqIO 1.2.0;
use BioUtils::FastaIO 1.2.0;
use Getopt::Long;
use Config::Std;
use Pod::Usage;
use Readonly;
use Cwd;

# Global Variables
Readonly::Scalar my $OUTPUT_DIR => getcwd;
Readonly::Scalar my $OUTPUT_PREFIX => "out";
Readonly::Scalar my $KEEP_TMP => "no";
Readonly::Scalar my $EVAL => 0.00001;
Readonly::Scalar my $PERC_IDEN => 94;
Readonly::Scalar my $OUTPUT_FMT => 6;
Readonly::Scalar my $MAX_TARGETS => 1;
Readonly::Scalar my $PARALLEL => "OFF";
Readonly::Scalar my $SEQS_PER_FILE => 10000;
Readonly::Scalar my $QUEUE => "bigmem";
Readonly::Scalar my $MEM => 100;

print "Start\n\n";

# Variables to be set by GetOpts and their defaults
my $params_file = undef;
my ($blast_db, $query_file, $fastq, $output_dir, $output_prefix, $keep_tmp,
    $eval, $perc_iden, $output_fmt, $max_targets, $otu_table, $parallel,
    $seqs_per_file, $queue, $mem, $help, $man);

my $options_okay = GetOptions (
    "params_file:s" => \$params_file,  # string
    "blast_db:s" => \$blast_db,
    "query_file:s" => \$query_file,
    "fastq" => \$fastq, # flag
    "output_dir:s" => \$output_dir,
    "output_prefix:s" => \$output_prefix,
    "keep_tmp:s" => \$keep_tmp,
    "eval:s" => \$eval,
    "perc_iden:s" => \$perc_iden,
    "output_fmt:s" => \$output_fmt,
    "max_targets:s" => \$max_targets,
    "otu_table:s" => \$otu_table,
    "parallel:s" => \$parallel,
    "seqs_per_file:s" => \$seqs_per_file,
    "queue:s" => \$queue,
    "mem:s" => \$mem,
    "help" => \$help,                  # flag
    "man" => \$man,  #flag
);

# Fail if unknown options are encountered
if ( ! $options_okay ) {
    pod2usage(2);
}

if ( $help ) {
    pod2usage(2);
}

if ( $man ) {
    pod2usage(-verbose => 2);
}

if ( ! defined $blast_db and ! defined $params_file ) {
    pod2usage(-message => "Supply a blast_db or params_file",
              -exit_status => 2,
              -verbose => 0)
}

if ( ! defined $query_file and ! defined $params_file ) {
    pod2usage(-message => "Supply a query_file or params_file",
              -exit_status => 2,
              -verbose => 0)
}

if ( ! defined $output_dir ) { $output_dir = $OUTPUT_DIR; }
if ( ! defined $output_prefix ) { $output_prefix = $OUTPUT_PREFIX; }
if ( ! defined $keep_tmp ) { $keep_tmp = $KEEP_TMP; }
if ( ! defined $eval ) { $eval = $EVAL; }
if ( ! defined $perc_iden ) { $perc_iden = $PERC_IDEN; }
if ( ! defined $output_fmt ) { $output_fmt = $OUTPUT_FMT; }
if ( ! defined $max_targets ) { $max_targets = $MAX_TARGETS; }
if ( ! defined $parallel ) { $parallel = $PARALLEL; }
if ( ! defined $seqs_per_file ) { $seqs_per_file = $SEQS_PER_FILE; }
if ( ! defined $queue ) { $queue = $QUEUE; }
if ( ! defined $mem ) { $mem = $MEM; }





if ( defined $params_file ) {
    print "Read Config File\n";
    read_config( $params_file => my %hash);
    
    $blast_db = $hash{'BLAST Parameters'}{'db'};
    $query_file = $hash{'BLAST Parameters'}{'query'};
    $output_dir = $hash{'Output'}{'dir'};
    $output_prefix = $hash{'Output'}{'prefix'};
    $keep_tmp = $hash{'Output'}{'keep_tmp'};
    $eval = $hash{'BLAST Parameters'}{'evalue'};
    $perc_iden = $hash{'BLAST Parameters'}{'perc_identity'};
    $output_fmt = $hash{'BLAST Parameters'}{'output_fmt'};
    $max_targets = $hash{'BLAST Parameters'}{'max_targets'};
    $otu_table = $hash{'OTU Table'}{'otu_table_file'};
    $parallel = $hash{'BLAST Parameters'}{'parallel'};
    $seqs_per_file = $hash{'BLAST Parameters'}{'seqs_per_file'};
    $queue = $hash{'BLAST Parameters'}{'queue'};
    $mem = $hash{'BLAST Parameters'}{'mem'};
}

# make a temp fasta file if necessary
if ( $fastq ) {
    my $tmp_fasta = $query_file . ".fasta.tmp";
    
    print "Making temporary fasta file: $tmp_fasta";
    my $in_fastq = BioUtils::FastqIO->new({
                            stream_type => '<',
                            file => "$query_file"
                        });
    my $out_fasta = BioUtils::FastaIO->new({
                            stream_type => '>',
                            file => "$tmp_fasta"
                        });
    
    while ( my $fastq_seq = $in_fastq->get_next_seq() ) {
        $out_fasta->write_seq($fastq_seq->to_FastaSeq);
    }
    
    $query_file = $tmp_fasta;
}
    
print "Building Filter\n\n";
my $filter = BioUtils::QC::ContaminantFilter->new({
                blast_db => $blast_db,
                query_file => $query_file,
                output_dir => $output_dir,
                output_prefix => $output_prefix,
                keep_tmp => $keep_tmp,
                eval => $eval,
                perc_iden => $perc_iden,
                output_fmt => $output_fmt,
                max_targets => $max_targets,
                otu_table => $otu_table
                });

if ( $parallel =~ m/ON/i ) {
    print "Running Parallelized Filter\n\n";
    $filter->run_filter(1, $seqs_per_file, $queue, $mem);
}
else {
    print "Running Filter\n\n";
    $filter->run_filter();
}

# remove temp fasta if neccessary
if ( $fastq ) {
    `rm $query_file`;
}

print "Finished\n\n";


__END__

# POD

=head1 NAME

filter_contaminants.pl - Filters out contaminant sequences using blastn


=head1 VERSION

This documentation refers to filter_contaminants.pl version 1.2.0


=head1 SYNOPSIS

perl filter_contaminants.pl
--params_file <file path>
--blast_db <blast db path>
--query_file <file path>
--fastq
--output_dir <dir path>
--output_prefix <prefix str>
--keep_tmp <on | off>
--eval <decimanl>
--perc_iden <int>
--output_fmt <int>
--max_targets <int>
--otu_table <file path>
--parallel <on | off>
--seqs_per_file <int>
--queue <bsub queue>
--mem <bsub -M argument>
--help
--man


=head1 ARGUMENTS

=head2 --params_file

A file path to a file which contains all necessary and optional parameters

=head2 --blast_db

Path to a blast formated database.  REQUIRED

=head2 --query_file

Path to file with query sequecnes in FASTA format.  REQUIRED

=head2 --fastq

Flag indicating that input query files are in FASTQ format.  A temporary query
FASTA file will be created for downstream analysis.  The temporary file will be
removed when analysis is complete.

=head2 --output_dir

Path to directory where output will be printed.  DEFAULT: working directory

=head2 --output_prefix

String to prepend output files.  DEFAULT:  out

=head2 --keep_tmp

"On" or "Off" indicating if the temporary files should be kept.  DEFAULT:  off

=head2 --eval

The minimum evalue to be considered a contaminant.  DEFAULT:  0.00001

=head2 --perc_iden

The minimum percent identity to be considered a contaminant.  DEFAULT:  94

=head2 --output_fmt

The blast output format.  DEFAULT:  6

=head2 --max_targets

The number of blast targets to keep for each query sequence.  DEFAULT:  1

=head2 --otu_table

If these sequecnes are associtated and OTU table any sequences that are classified
as contaminants will be seperated into a seperate OTU table.  DEFAULT: NA

=head2 --parallel

"On" or "Off" indicating a parallel run on an LSF cluster.  The source code will
have to be modified if you are not using UNC's Kure or Killdevil cluster.
DEFAULT:  Off

=head2 --seqs_per_file

For parallel runs, the number of sequence to put in each file for the parallel
analysis.  DEFAULT: 10000

=head2 --queue

The bsub queue to which BLAST jobs should be submitted.  DEFAULT: bigmem

=head2 --mem

The bsub -M argument--amount of memory allocated for each BLAST job.  DEFAUL: 100


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
    parallel: ON
    seqs_per_file: 10000
    queue: bigmem
    mem: 100
    
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

    BioUtils::QC::ContaminantFilter 1.2.0
    Getopt::Long
    Config::Std
    Pod::Usage
    Readonly
    Cwd


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

