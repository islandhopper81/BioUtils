
use strict;
use warnings;

use BioUtils::Align::Pairwise::Alignment;
use BioUtils::Align::FastaSeq;
use Test::More tests => 38;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::Align::Pairwise::Alignment' ); }

# create some seqs for testing
my $seq1 = BioUtils::Align::FastaSeq->new({
                header => 'seq1',
                seq => "ATCG",
                start => 10,
                end => 14});
my $seq2 = BioUtils::Align::FastaSeq->new({
                header => 'seq2',
                seq => "ATCG",
                start => 10,
                end => 14});

# test constructor
my $aln;
{
    throws_ok( sub{ BioUtils::Align::Pairwise::Alignment->new() },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::Pairwise::Alignment->new({
                        seq1 => $seq1}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::Pairwise::Alignment->new({
                        seq1 => $seq1,
                        seq2 => $seq2,}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::Pairwise::Alignment->new({
                        seq1 => $seq1,
                        seq2 => $seq2,
                        score => 10}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::Pairwise::Alignment->new({
                        seq1 => $seq1,
                        seq2 => $seq2,
                        score => 10,
                        mismatch_count => 0}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::Pairwise::Alignment->new({
                        seq1 => $seq1,
                        seq2 => $seq2,
                        score => 10,
                        mismatch_count => 0,
                        indel_count => 0}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    lives_ok( sub{ $aln = BioUtils::Align::Pairwise::Alignment->new({
                    seq1 => $seq1,
                    seq2 => $seq2,
                    score => 10,
                    perc_iden => 1,
                    mismatch_count => 1,
                    indel_count => 0}) },
             "expected to live" );
}

# test get_score
{
    is( $aln->get_score(), 10, "get_score()" );
}

# test get_seq1
{
    is( $aln->get_seq1()->get_header(), 'seq1', "get_seq1()" );
}

# test get_seq2
{
    is( $aln->get_seq2()->get_header(), 'seq2', "get_seq2()" );
}

# test get_seq
{
    is( $aln->get_seq("seq1")->get_header(), 'seq1', "get_seq(seq1)" );
}

# test get_perc_iden
{
    is( $aln->get_perc_iden(), 1, "get_perc_iden() - 1" );
}

# test get_mismatch_count
{
    is( $aln->get_mismatch_count(), 1, "get_mismatch_count() - 1" );
}

# test get_indel_count
{
    is( $aln->get_indel_count(), 0, "get_indel_count() - 0" );
}

# test seq_set1
{
    throws_ok( sub{ $aln->set_seq1() },
              'MyX::Generic::Undef::Param', "set_seq1() - caught" );
    lives_ok( sub{ $aln->set_seq1($seq2) },
              "expected to live" );
    is ( $aln->get_seq1()->get_header(), 'seq2', "set_seq1 test" );
}

# test seq_set2
{
    throws_ok( sub{ $aln->set_seq2() },
              'MyX::Generic::Undef::Param', "set_seq2() - caught" );
    lives_ok( sub{ $aln->set_seq2($seq1) },
              "expected to live" );
    is ( $aln->get_seq2()->get_header(), 'seq1', "set_seq2 test" );
}

# test set_perc_iden
{
    throws_ok( sub{ $aln->set_perc_iden() },
              'MyX::Generic::Undef::Param', "set_perc_iden() - caught" );
    throws_ok( sub{ $aln->set_perc_iden('a') },
              'MyX::Generic::Digit::MustBeDigit',
              "set_perc_iden(a) - caught" );
    throws_ok( sub{ $aln->set_perc_iden(-1) },
              'MyX::Generic::Digit::OOB',
              "set_perc_iden(-1) - caught" );
    lives_ok( sub{ $aln->set_perc_iden(.3) },
              "expected to live" );
    is ( $aln->get_perc_iden(), .3, "set_indel_count(.3) test" );
}

# test set_mismatch_count
{
    throws_ok( sub{ $aln->set_mismatch_count() },
              'MyX::Generic::Undef::Param', "set_mismatch_count() - caught" );
    throws_ok( sub{ $aln->set_mismatch_count('a') },
              'MyX::Generic::Digit::MustBeDigit',
              "set_mismatch_count(a) - caught" );
    throws_ok( sub{ $aln->set_mismatch_count(-1) },
              'MyX::Generic::Digit::TooSmall',
              "set_mismatch_count(-1) - caught" );
    lives_ok( sub{ $aln->set_mismatch_count(3) },
              "expected to live" );
    is ( $aln->get_mismatch_count(), 3, "set_mismatch_count(3) test" );
}

# test set_indel_count
{
    throws_ok( sub{ $aln->set_indel_count() },
              'MyX::Generic::Undef::Param', "set_indel_count() - caught" );
    throws_ok( sub{ $aln->set_indel_count('a') },
              'MyX::Generic::Digit::MustBeDigit',
              "set_indel_count(a) - caught" );
    throws_ok( sub{ $aln->set_indel_count(-1) },
              'MyX::Generic::Digit::TooSmall',
              "set_indel_count(-1) - caught" );
    lives_ok( sub{ $aln->set_indel_count(3) },
              "expected to live" );
    is ( $aln->get_indel_count(), 3, "set_indel_count(3) test" );
}

# test get_str
{
    my $str = "######\n";
    $str .= "+seq1\n";
    $str .= "*seq2\n";
    $str .= "\tScore: 10\n";
    $str .= "\tPercent Identity: 0.9\n";
    $str .= "\tMismatch Count: 0\n";
    $str .= "\tIndel Count: 0\n";
    $str .= "\t+Start: 10\n";
    $str .= "\t+End: 14\n";
    $str .= "\t*Start: 10\n";
    $str .= "\t*End: 14\n";
    $str .= "\n";
    $str .= "+" . $seq1->get_seq() . "\n";
    $str .= "*" . $seq2->get_seq() . "\n";
    $str .= "______\n";
    
    my $seq1 = BioUtils::Align::FastaSeq->new({
                header => 'seq1',
                seq => "ATCG",
                start => 10,
                end => 14});
    my $seq2 = BioUtils::Align::FastaSeq->new({
                header => 'seq2',
                seq => "ATCG",
                start => 10,
                end => 14});
    
    my $aln;
    lives_ok( sub{ $aln = BioUtils::Align::Pairwise::Alignment->new({
                    seq1 => $seq1,
                    seq2 => $seq2,
                    score => 10,
                    perc_iden => 0.9,
                    mismatch_count => 0,
                    indel_count => 0}) },
             "expected to live" );
    
    is( $aln->get_str(), $str, "get_str()" );
}
