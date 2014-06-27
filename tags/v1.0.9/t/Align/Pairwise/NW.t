
use strict;
use warnings;

use BioUtils::Align::Pairwise::NW;
use BioUtils::FastaSeq;
use Test::More tests => 12;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::Align::Pairwise::NW' ); }

### Create 2 FastqSeq objects to use for testing
# These seqs come from http://amrita.vlab.co.in/?sub=3&brch=274&sim=1431&cnt=1
my $h1 = 'Seq1';
my $s1 = "CGTGAATTCAT";
my $seq1 = BioUtils::FastaSeq->new({header => $h1, seq => $s1});

my $h2 = 'Seq2';
my $s2 = "GACTTAC";
my $seq2 = BioUtils::FastaSeq->new({header => $h2, seq => $s2});

# Test NW object creation
my $nw;
{
    lives_ok( sub { $nw = BioUtils::Align::Pairwise::NW->new() },
             "creating new NW obj");
    lives_ok( sub { $nw = BioUtils::Align::Pairwise::NW->new({
                            match_score => 1,
                            mismatch_score => -1,
                            gap_score => -1,
                            })
                    },
             "creating new NW obj with parameters"
            );
}

# Test Align method
{
    my $aln;
    lives_ok( sub{ $aln = $nw->align($seq1, $seq2) },
             "expected to live" );
    is( $aln->get_score(), -1, "get_score() = -1" );
    is( $aln->get_perc_iden(), (5/11), "get_perc_iden() = (5/11)" );
    is( $aln->get_seq1()->get_start(), 0, "seq1->get_start()" );
    is( $aln->get_seq1()->get_end(), 10, "seq1->get_end()" );
    is( $aln->get_seq2()->get_start(), 0, "seq2->get_start()" );
    is( $aln->get_seq2()->get_end(), 6, "seq2->get_end()" );
    is( $aln->get_mismatch_count(), 2, "get_mismatch_count() - 2" );
    is( $aln->get_indel_count(), 4, "get_indel_count() - 4" );
}