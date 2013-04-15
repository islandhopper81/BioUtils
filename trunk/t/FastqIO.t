use strict;
use warnings;

use BioUtils::FastqIO 0.0.2;
use BioUtils::FastqSeq;
use Test::More tests => 18;
use Test::Exception;
use File::Temp qw{tempfile tempdir};

BEGIN { use_ok( 'BioUtils::FastqIO' ); }

# I need to have a temporary file to use for testing ouput and input of BioUtils::FastqIO
my($fh, $filename) = tempfile();

# Test constructor
my $out_parser = undef;
throws_ok( sub { BioUtils::FastqIO->new({file => $filename}) },
          'MyX::Generic::Undef::Param', "new() - caught" );
throws_ok( sub { BioUtils::FastqIO->new({stream_type => 'out'}) },
          'MyX::Generic::Undef::Param', "new() - caught" );
lives_ok( sub { $out_parser = BioUtils::FastqIO->new({
                                stream_type => '>', file => $filename}) },
         "expected to live" );


# First test write_seq method on following BioUtils::FastqSeq objects
my $fastq_seq1 = BioUtils::FastqSeq->new({
                    header => 'Seq1', seq => 'ATCG', quals_str => '9876'});
lives_ok( sub { $out_parser->write_seq($fastq_seq1) }, "expected to live" );
my $fastq_seq2 = BioUtils::FastqSeq->new({
                    header => 'Seq2', seq => 'GGTA', quals_str => '1234'});
lives_ok( sub { $out_parser->write_seq($fastq_seq2) }, "expected to live" );

# try printing a sequence without a header.  It should rethrow an error
my $fastq_seq3 = BioUtils::FastqSeq->new({seq => 'TTGA', quals_str => '2938'});
throws_ok( sub { $out_parser->write_seq($fastq_seq3); },
          'MyX::Generic::Undef::Attribute', 'write_seq() - caught' );
dies_ok( sub { $out_parser->write_seq($fastq_seq3) }, 'expected to die' );

# Detele the first parser so that the file handle to that file will be free to
# read with the second parser
lives_ok( sub { $out_parser->DESTROY }, "expected to live");


# Now read the temp.fastq file checking its contents
my $in_parser = undef;
lives_ok( sub { $in_parser = BioUtils::FastqIO->new({
                                stream_type => "<", file => $filename}) },
         "expected to live" );
my $got_seq = undef;
lives_ok( sub { $got_seq = $in_parser->get_next_seq() }, "expected to live" );
is( $got_seq->get_header(), 'Seq1', "get_next_seq()->get_header()" );
is( $got_seq->get_seq(), 'ATCG', "get_next_seq()->get_seq()" );
is( $got_seq->get_quals_str(), '9876', "get_next_seq()->get_quals_str()" );

# Now try on the second sequence
lives_ok( sub { $got_seq = $in_parser->get_next_seq() }, "expected to live" );
is( $got_seq->get_header(), 'Seq2', "get_next_seq()->get_header()" );
is( $got_seq->get_seq(), 'GGTA', "get_next_seq()->get_seq()" );
is( $got_seq->get_quals_str(), '1234', "get_next_seq()->get_quals_str()" );

# Now try on a third non-existing sequence
$in_parser->get_next_seq();

