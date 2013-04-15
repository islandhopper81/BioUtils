
use strict;
use warnings;

use BioUtils::FastaSeq;
use Test::More tests => 16;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::FastaSeq' ); }

### Create a FastqSeq object to use for testing
my $header = 'Seq1 Homo Sapien';
my $seq = "ATCG";
my $fasta_seq = BioUtils::FastaSeq->new({header => $header, seq => $seq});

# Test the simple getter methods
lives_ok( sub { $fasta_seq->get_header() }, "expected to live" );
is ($fasta_seq->get_header(), $header, "get_header()" );
lives_ok( sub { $fasta_seq->get_id() }, "expected to live" );
is( $fasta_seq->get_id(), 'Seq1', "get_id()" );
is( $fasta_seq->get_seq(), $seq, "get_seq()" );

# Test the simple setter methods
my $new_seq = "TCGA";
my $new_header = "my_seq";
lives_ok( sub{ $fasta_seq->set_header($new_header) }, "expected to live" );
lives_ok( sub{ $fasta_seq->set_seq($new_seq) }, "expected to live" );
is( $fasta_seq->get_seq(), $new_seq, "set_seq()" );
is( $fasta_seq->get_header(), $new_header, "set_header()" );


# Create a sequence with an empty header 
lives_ok( sub { $fasta_seq = BioUtils::FastaSeq->new({ seq => $seq}) },
          "expected to live" );
dies_ok( sub { $fasta_seq->get_header() }, "expected to die" );
throws_ok( sub { $fasta_seq->get_header() }, 'MyX::Generic::Undef::Attribute', "get_header() - caught" );
throws_ok( sub { $fasta_seq->get_header() }, qr/Undefined header/, "get_header() - caught" );
throws_ok( sub { $fasta_seq->get_id() }, 'MyX::Generic::Undef::Attribute', "get_id() - caught" );
throws_ok( sub { $fasta_seq->get_id() }, qr/Undefined header/, "get_id() - caught" );

