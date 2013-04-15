
use strict;
use warnings;

use FastqSeq;
use Test::More tests => 35;
use Test::Exception;

BEGIN { use_ok( 'FastqSeq' ); }

### Create a FastqSeq object to use for testing
my $header = '@Seq1 Homo Sapien';
my $seq = "ATCG";
my $quals_str = "1234";
my $fastq_seq = FastqSeq->new({header => $header, seq => $seq, quals_str => $quals_str});

# Test the simple getter methods
lives_ok( sub { $fastq_seq->get_header() }, "expected to live" );
is ($fastq_seq->get_header(), $header, "get_header()" );
lives_ok( sub { $fastq_seq->get_id() }, "expected to live" );
is( $fastq_seq->get_id(), '@Seq1', "get_id()" );
is( $fastq_seq->get_seq(), $seq, "get_seq()" );
is( $fastq_seq->get_quals_str(), $quals_str, "get_quals_str()" );
my @quals = split //, $quals_str;
is_deeply( $fastq_seq->get_quals_aref(), \@quals, "get_quals_aref()" );

# Test the simple setter methods
my $new_seq = "TCGA";
my $new_quals_str = "9876";
$fastq_seq->set_seq($new_seq);
$fastq_seq->set_quals_str($new_quals_str);
is( $fastq_seq->get_seq(), $new_seq, "set_seq()" );
is( $fastq_seq->get_quals_str(), $new_quals_str, "set_quals_str()" );

# test to_FastaSeq
{
    my $fastq_seq = FastqSeq->new({header => $header,
                                   seq => $seq,
                                   quals_str => $quals_str});
    
    my $expected = FastaSeq->new({header => "Seq1 Homo Sapien",
                                  seq => "ATCG"});
    my $fasta_seq;
    lives_ok( sub{ $fasta_seq = $fastq_seq->to_FastaSeq() },
             "to_FastqSeq - lives" );
    is( $fasta_seq->get_header(),
       "Seq1 Homo Sapien",
       "fasta_seq->get_header()" );
    is( $fasta_seq->get_seq(),
       "ATCG",
       "fasta_seq->get_seq()" );
}

# Test the quality indexer method.  Remember, I am using the new_seq values
lives_ok (sub { $fastq_seq->get_qual_at(0) }, "expected to live" );
is( $fastq_seq->get_qual_at(0), 24, "get_qual_at(0)");  # encoded 9 == 51
is( $fastq_seq->get_qual_at(1), 23, "get_qual_at(1)");  # encoded 8 == 50
is( $fastq_seq->get_qual_at(2), 22, "get_qual_at(2)");  # encoded 7 == 49
is( $fastq_seq->get_qual_at(3), 21, "get_qual_at(3)");  # encoded 6 == 48
dies_ok( sub { $fastq_seq->get_qual_at(-1) }, "expected to die" );
throws_ok( sub { $fastq_seq->get_qual_at(-1) }, 'MyX::Generic::OutOfBounds', "get_qual_at(-1) - caught" );
throws_ok( sub { $fastq_seq->get_qual_at(-1) }, qr/Index out of bounds/, "get_qual_at(-1) - caught" );
throws_ok( sub { $fastq_seq->get_qual_at(4) }, 'MyX::Generic::OutOfBounds', "get_qual_at(4) - caught" );
throws_ok( sub { $fastq_seq->get_qual_at(4) }, qr/Index out of bounds/, "get_qual_at(4) - caught" );

# some of my own experimentation with exceptions
#eval { $fastq_seq->get_qual_at(-1); };
#if ( my $e = MyX::Generic::OutOfBounds->caught() ) {
#    print $e->error(), "\n";
#    print $e->index(), "\n";
#}
#eval { $fastq_seq->get_qual_at(4); };
#if ( my $e = MyX::Generic::OutOfBounds->caught() ) {
#    print "error: ", $e->error(), "\n";
#    print "message: ", $e->message(), "\n";
#    print "index: ", $e->index(), "\n";
#    print "file: ", $e->file(), "\n";
#    print "line: ", $e->line(), "\n";
#    print "trace: ", $e->trace(), "\n";
#    print "description: ", $e->description(), "\n";
#    print "as_string: ", $e->as_string(), "\n";
#    print "full_message: ", $e->full_message(), "\n";
#}


# Create a sequence with an empty header 
lives_ok( sub { $fastq_seq = FastqSeq->new({ seq => $seq, quals_str => $quals_str }) },
          "expected to live" );
dies_ok( sub { $fastq_seq->get_header() }, "expected to die" );
throws_ok( sub { $fastq_seq->get_header() }, 'MyX::Generic::Undef::Attribute', "get_header() - caught" );
throws_ok( sub { $fastq_seq->get_header() }, qr/Undefined header/, "get_header() - caught" );
throws_ok( sub { $fastq_seq->get_id() }, 'MyX::Generic::Undef::Attribute', "get_id() - caught" );
throws_ok( sub { $fastq_seq->get_id() }, qr/Undefined header/, "get_id() - caught" );


# Test the _encoding_to_dec method
is( FastqSeq::_encoding_to_dec('I'), 40, "_encoding_to_dec(I)" );
is( FastqSeq::_encoding_to_dec('J'), 41, "_encoding_to_dec(J)" );
is( FastqSeq::_encoding_to_dec('!'), 0, "_encoding_to_dec(!)" );
# I still need to test out of range query errors

# Test the _dec_to_encoding method
is( FastqSeq::_dec_to_encoding(40), 'I', "_dec_to_encoding(40)" );
is( FastqSeq::_dec_to_encoding(41), 'J', "_dec_to_encoding(41)" );
is( FastqSeq::_dec_to_encoding(0), '!', "_dec_to_encoding(0)" );
# I still need to test out of range query errors


