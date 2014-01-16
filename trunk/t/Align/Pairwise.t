
use strict;
use warnings;

use BioUtils::Align::Pairwise;
use BioUtils::FastaSeq;
use Test::More tests => 18;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::Align::Pairwise' ); }

### Create 2 FastqSeq objects to use for testing
my $h1 = 'Seq1';
my $s1 = "CGTGAATTCAT";
my $seq1 = BioUtils::FastaSeq->new({header => $h1, seq => $s1});

my $h2 = 'Seq2';
my $s2 = "GACTTAC";
my $seq2 = BioUtils::FastaSeq->new({header => $h2, seq => $s2});

# Test NW object creation
my $obj;
{
    lives_ok( sub { $obj = BioUtils::Align::Pairwise->new() },
             "creating new Pairwise obj");
    lives_ok( sub { $obj = BioUtils::Align::Pairwise->new({
                            match_score => 5,
                            mismatch_score => -1,
                            gap_score => -2,
                            })
                    },
             "creating new Pairwise obj with parameters"
            );
}

# Test Getters
{
    is( $obj->get_match_score(), 5, "get_match_score" );
    is( $obj->get_mismatch_score(), -1, "get_mismatch-score" );
    is( $obj->get_gap_score(), -2, "get_gap_score" );
}

# Test Setters
{
    lives_ok( sub {$obj->set_match_score(3)}, "set_match_score - lives" );
    is( $obj->get_match_score(), 3, "get_match_score(3)" );
    lives_ok( sub {$obj->set_mismatch_score(0)}, "set_mismatch_score - lives" );
    is( $obj->get_mismatch_score(), 0, "get_mismatch_score(0)" );
    lives_ok( sub {$obj->set_gap_score(-1)}, "set_gap_score - lives" );
    is( $obj->get_gap_score(), -1, "get_gap_score(-1)" );
    
    throws_ok( sub{ $obj->set_match_score('a') },
              "MyX::Generic::Digit::MustBeDigit",
              "set_match_score(a) - caught" );
    throws_ok( sub{ $obj->set_match_score('a') },
              qr/Match Score must be a digit/,
              "set_match_score(a) - caught" );
    throws_ok( sub{ $obj->set_mismatch_score('a') },
              "MyX::Generic::Digit::MustBeDigit",
              "set_mismatch_score(a) - caught" );
    throws_ok( sub{ $obj->set_mismatch_score('a') },
              qr/Mismatch Score must be a digit/,
              "set_mismatch_score(a) - caught" );
    throws_ok( sub{ $obj->set_gap_score('a') },
              "MyX::Generic::Digit::MustBeDigit",
              "set_gap_score(a) - caught" );
    throws_ok( sub{ $obj->set_gap_score('a') },
              qr/Gap Score must be a digit/,
              "set_gap_score(a) - caught" );
}
