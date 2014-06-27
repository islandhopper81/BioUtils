#!/usr/bin/env perl

use strict;
use warnings;

use Benchmark qw(:all);

use BioUtils::Align::Pairwise::SW;

my $sw = BioUtils::Align::Pairwise::SW->new();

# create seqs to align
my $h1 = 'Seq1';
my $s1 = "CGTGAATTCAT";
my $seq1 = BioUtils::FastaSeq->new({header => $h1, seq => $s1});

my $h2 = 'Seq2';
my $s2 = "GACTTAC";
my $seq2 = BioUtils::FastaSeq->new({header => $h2, seq => $s2});


# benchmark the alignment
cmpthese(100000, {
	'New' => sub{ $sw->align($seq1, $seq2) },
	'Old' => sub{ $sw->align_old($seq1, $seq2) },
});
