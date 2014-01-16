
use strict;
use warnings;

use BioUtils::Align::FastaSeq;
use Test::More tests => 12;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::Align::FastaSeq' ); }

# test constructor
my $aln_seq;
{
    throws_ok( sub{ BioUtils::Align::FastaSeq->new() },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::FastaSeq->new({
                        header => "seq1"}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::FastaSeq->new({
                        header => "seq1",
                        seq => "ATCT",}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub{ BioUtils::Align::FastaSeq->new({
                        header => "seq1",
                        seq => "ATCT",
                        start => 10}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    lives_ok( sub{ $aln_seq = BioUtils::Align::FastaSeq->new({
                    header => "seq1",
                    seq => "ATCT",
                    start => 10,
                    end => 14}) },
             "expected to live" );
}

# test get_start
{
    is( $aln_seq->get_start(), 10, "get_start()" );
}

# test get_end
{
    is( $aln_seq->get_end(), 14, "get_end()" );
}

# test set_start
{
    lives_ok( sub{ $aln_seq->set_start(11) },
              "expected to live" );
    is ( $aln_seq->get_start(), 11, "get_start() test" );
}

# test set_end
{
    lives_ok( sub{ $aln_seq->set_end(15) },
              "expected to live" );
    is ( $aln_seq->get_end(), 15, "get_end() test" );
}
