
use strict;
use warnings;

use BioUtils::ConsensusBuilder::FastqConsensus;
use Test::More tests => 8;
use Test::Exception;


BEGIN { use_ok( 'BioUtils::ConsensusBuilder::FastqConsensus'); }

### Create a FastqColumn object to use for testing
my $seq = "ATCG";
my $quals_str = "AAAA";
my $fc;

# test constructor with only required parameters
{
    lives_ok( sub { $fc = BioUtils::ConsensusBuilder::FastqConsensus->new({
                                seq => $seq,
                                quals_str => $quals_str,
                            } ); },
             "expected to live" );
}

# test with the optional c_score parameter
{
    lives_ok( sub { $fc = BioUtils::ConsensusBuilder::FastqConsensus->new({
                                seq => $seq,
                                quals_str => $quals_str,
                                c_score => 40,
                            } ); },
             "expected to live with c_score parameter" );
}

# test with WRONG optional c_score parameter
{
    throws_ok( sub { $fc = BioUtils::ConsensusBuilder::FastqConsensus->new({
                                seq => $seq,
                                quals_str => $quals_str,
                                c_score => 'a',
                            } ); },
             'MyX::Generic::Digit::MustBeDigit', "new(c_score => a) - caught" );
}

# test get_c_score and set_c_score
{
    $fc = BioUtils::ConsensusBuilder::FastqConsensus->new({
                                seq => $seq,
                                quals_str => $quals_str,
                            } );
    throws_ok( sub{ $fc->get_c_score() }, 'MyX::Generic::Undef::Param',
              "get_c_score(undef) - caught" );
    
    throws_ok( sub{ $fc->set_c_score('c') }, 'MyX::Generic::Digit::MustBeDigit',
              "set_c_score(c) - caught" );
    
    lives_ok( sub{ $fc->set_c_score(40) }, "set_c_score(40) - expected to live");
    is( $fc->get_c_score(), 40, "get_c_score() == 40" );
}



