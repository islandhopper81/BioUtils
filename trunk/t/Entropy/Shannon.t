use strict;
use warnings;

use BioUtils::Entropy::Shannon;
use BioUtils::FastaIO;
use BioUtils::FastaSeq;
use Test::More tests => 6;
use Test::Exception;
use File::Temp qw{tempfile tempdir};

BEGIN { use_ok( 'BioUtils::Entropy::Shannon' ); }

# I need to have a temporary file to use for testing ouput and input of FastqParser
my($fh, $filename) = tempfile();

# Test constructor
my $shannon = undef;
lives_ok( sub { $shannon = BioUtils::Entropy::Shannon->new() },
         "expected to live" );


# test compute_str
{
    is( $shannon->compute_str("ATCG"), 2, "compute_str(ATCG)" );
    is( $shannon->compute_str("AAAA"), 0, "compute_str(AAAA)" );
    is( $shannon->compute_str("ATGGTA"), 1.58496250072116,
       "compute_str(ATGGTA)" );
}


# test compute_seq
{
    my $seq = BioUtils::FastaSeq->new({
        header => "seq1",
        seq => "ATCG" });
    
    is( $shannon->compute_seq($seq), 2, "compute_seq(ATCG)" );
}

# test compute_file -- not tested because it prints to screen.
# I can fix this in the future