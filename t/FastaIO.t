use strict;
use warnings;

use BioUtils::FastaIO;
use BioUtils::FastaSeq;
use Test::More tests => 17;
use Test::Exception;
use File::Temp qw{tempfile tempdir};

BEGIN { use_ok( 'BioUtils::FastaIO' ); }

# I need to have a temporary file to use for testing ouput and input of FastqParser
my($fh, $filename) = tempfile();

# Test constructor
my $fasta_out = undef;
throws_ok( sub { BioUtils::FastaIO->new({file => $filename}) },
          'MyX::Generic::Undef::Param', "new() - caught" );
throws_ok( sub { BioUtils::FastaIO->new({stream_type => 'out'}) },
          'MyX::Generic::Undef::Param', "new() - caught" );
lives_ok( sub { $fasta_out = BioUtils::FastaIO->new({
                    stream_type => '>', file => $filename}) },
         "expected to live" );


# First test write_seq method on following FastqSeq objects
my $fasta_seq1 = BioUtils::FastaSeq->new({header => 'Seq1', seq => 'ATCG'});
lives_ok( sub { $fasta_out->write_seq($fasta_seq1) }, "expected to live" );
my $fasta_seq2 = BioUtils::FastaSeq->new({header => 'Seq2', seq => 'GGTA', quals_str => '1234'});
lives_ok( sub { $fasta_out->write_seq($fasta_seq2) }, "expected to live" );

# try printing a sequence without a header.  It should rethrow an error
my $fasta_seq3 = BioUtils::FastaSeq->new({seq => 'TTGA'});
throws_ok( sub { $fasta_out->write_seq($fasta_seq3); },
          'MyX::Generic::Undef::Attribute', 'write_seq() - caught' );
dies_ok( sub { $fasta_out->write_seq($fasta_seq3) },
        'expected to die' );

# Detele the first parser so that the file handle to that file will be free to
# read with the second parser
lives_ok( sub { $fasta_out->DESTROY }, "expected to live");


# Now read the temp.fastq file checking its contents
my $fasta_in = undef;
lives_ok( sub { $fasta_in = BioUtils::FastaIO->new({
                                stream_type => "<", file => $filename}) },
         "expected to live" );
my $got_seq = undef;
lives_ok( sub { $got_seq = $fasta_in->get_next_seq() },
         "expected to live" );
is( $got_seq->get_header(), 'Seq1',
   "get_next_seq()->get_header()" );
is( $got_seq->get_seq(), 'ATCG',
   "get_next_seq()->get_seq()" );

# Now try on the second sequence
lives_ok( sub { $got_seq = $fasta_in->get_next_seq() }, "expected to live" );
is( $got_seq->get_header(), 'Seq2', "get_next_seq()->get_header()" );
is( $got_seq->get_seq(), 'GGTA', "get_next_seq()->get_seq()" );

# Now try on a third non-existing sequence
lives_ok( sub{ $fasta_in->get_next_seq() }, 'expected to live' );

