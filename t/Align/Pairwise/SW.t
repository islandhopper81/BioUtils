
use strict;
use warnings;

use BioUtils::Align::Pairwise::SW;
use BioUtils::FastaSeq;
use Test::More tests => 12;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::Align::Pairwise::SW' ); }

### Create 2 FastqSeq objects to use for testing
# these seqs come from http://www.bioinformatics.wsu.edu/bioinfo_course/notes/Lecture6.pdf
my $h1 = 'Seq1';
my $s1 = "AATGT";
my $seq1 = BioUtils::FastaSeq->new({header => $h1, seq => $s1});

my $h2 = 'Seq2';
my $s2 = "ATGAC";
my $seq2 = BioUtils::FastaSeq->new({header => $h2, seq => $s2});

# Test SW object creation
my $sw;
{
    lives_ok( sub { $sw = BioUtils::Align::Pairwise::SW->new() },
             "creating new SW obj");
    lives_ok( sub { $sw = BioUtils::Align::Pairwise::SW->new({
                            match_score => 1,
                            mismatch_score => -1,
                            gap_score => -2,
                            })
                    },
             "creating new SW obj with parameters"
            );
}

# Test Align method
{
    my $aln;
    lives_ok( sub{ $aln = $sw->align($seq1, $seq2) },
             "expected to live" );
    is( $aln->get_score(), 3, "get_score() = 3" );
    is( $aln->get_perc_iden(), 1, "get_perc_iden() = 1" );
    is( $aln->get_seq1()->get_start(), 1, "seq1->get_start()" );
    is( $aln->get_seq1()->get_end(), 3, "seq1->get_end()" );
    is( $aln->get_seq2()->get_start(), 0, "seq2->get_start()" );
    is( $aln->get_seq2()->get_end(), 2, "seq2->get_end()" );
    is( $aln->get_mismatch_count(), 0, "get_mismatch_count()" );
    is( $aln->get_indel_count(), 0, "get_indel_count()" );
}