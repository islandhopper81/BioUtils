
use strict;
use warnings;

use BioUtils::FastqSeq;
use Test::More tests => 91;
use Test::Exception;

BEGIN { use_ok( 'BioUtils::FastqSeq' ); }

### Create a FastqSeq object to use for testing
my $header = '@Seq1 Homo Sapien';
my $seq = "ATCG";
my $quals_str = "1234";
my $fastq_seq = BioUtils::FastqSeq->new({
    header => $header,
    seq => $seq,
    quals_str => $quals_str
});

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
my $new_header = "seq2";
is( $fastq_seq->set_header($new_header), 1, "set_header->(seq2)" );
is( $fastq_seq->get_header(), $new_header, "get_header == seq2" );

# test to_FastaSeq
{
    my $fastq_seq = BioUtils::FastqSeq->new({header => $header,
                                             seq => $seq,
                                             quals_str => $quals_str});
    
    my $expected = BioUtils::FastaSeq->new({header => "Seq1 Homo Sapien",
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
lives_ok( sub { $fastq_seq = BioUtils::FastqSeq->new({ seq => $seq,
                                                       quals_str => $quals_str
                                                       }) },
          "expected to live" );
dies_ok( sub { $fastq_seq->get_header() }, "expected to die" );
throws_ok( sub { $fastq_seq->get_header() }, 'MyX::Generic::Undef::Attribute', "get_header() - caught" );
throws_ok( sub { $fastq_seq->get_header() }, qr/Undefined header/, "get_header() - caught" );
throws_ok( sub { $fastq_seq->get_id() }, 'MyX::Generic::Undef::Attribute', "get_id() - caught" );
throws_ok( sub { $fastq_seq->get_id() }, qr/Undefined header/, "get_id() - caught" );


# Test the _encoding_to_dec method
is( BioUtils::FastqSeq::_encoding_to_dec('I'), 40, "_encoding_to_dec(I)" );
is( BioUtils::FastqSeq::_encoding_to_dec('J'), 41, "_encoding_to_dec(J)" );
is( BioUtils::FastqSeq::_encoding_to_dec('!'), 0, "_encoding_to_dec(!)" );
# I still need to test out of range query errors

# Test the _dec_to_encoding method
is( BioUtils::FastqSeq::_dec_to_encoding(40), 'I', "_dec_to_encoding(40)" );
is( BioUtils::FastqSeq::_dec_to_encoding(41), 'J', "_dec_to_encoding(41)" );
is( BioUtils::FastqSeq::_dec_to_encoding(0), '!', "_dec_to_encoding(0)" );
# I still need to test out of range query errors


# test trim_front method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $qual = "AAAABBBB";
    
    my $seq_obj = BioUtils::FastqSeq->new({
        header => $header,
        seq => $seq,
        quals_str => $qual,
    });
    
    my $trimmed_obj;
    
    # test when no parameter is given
    lives_ok( sub{ $trimmed_obj = $seq_obj->trim_front() },
              "trim_front() - lives" );
    is( $trimmed_obj->get_seq(), "", "trim_front() - check trimmed seq" );
    is( $trimmed_obj->get_quals_str(), "", "trim_front() - check qual seq" );
    is( $seq_obj->get_seq(), "ATCGATCG", "trim_front() - check kept seq" );
    is( $seq_obj->get_quals_str(), "AAAABBBB", "trim_front() - check kept quals" );
    
    # test when a string is given as a param
    throws_ok( sub{ $seq_obj->trim_front('a') },
              'MyX::Generic::Digit::MustBeDigit',
              "trim_front(a) - throws" );
    throws_ok( sub{ $seq_obj->trim_front('a') },
              qr/trim_front requires digit > 0/,
              "trim_front(a) - throws" );
    
    # test when a negative number is given as a param
    throws_ok( sub{ $seq_obj->trim_front(-1) },
              'MyX::Generic::Digit::TooSmall',
              "trim_front(-1) - throws" );
    throws_ok( sub{ $seq_obj->trim_front(-1) },
              qr/trim_front requires digit > 0/,
              "trim_front(-1) - throws" );
    
    # test a good trim run
    lives_ok( sub{ $trimmed_obj = $seq_obj->trim_front(2) },
             "trim_front(2) - lives" );
    is( $trimmed_obj->get_seq(), "AT", "trim_front(2) - check trimmed seq" );
    is( $trimmed_obj->get_quals_str(), "AA", "trim_front(2) - check qual seq" );
    is( $seq_obj->get_seq(), "CGATCG", "trim_front(2) - check kept seq" );
    is( $seq_obj->get_quals_str(), "AABBBB", "trim_front(2) - check kept quals" );
}

# test trim_back method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $qual = "AAAABBBB";
    
    my $seq_obj = BioUtils::FastqSeq->new({
        header => $header,
        seq => $seq,
        quals_str => $qual,
    });
    
    my $trimmed_obj;
    
    # test when no parameter is given
    lives_ok( sub{ $trimmed_obj = $seq_obj->trim_back() },
              "trim_back() - throws" );
    is( $trimmed_obj->get_seq(), "", "trim_back(2) - check trimmed seq" );
    is( $trimmed_obj->get_quals_str(), "", "trim_back(2) - check qual seq" );
    is( $seq_obj->get_seq(), "ATCGATCG", "trim_back(2) - check kept seq" );
    is( $seq_obj->get_quals_str(), "AAAABBBB", "trim_back(2) - check kept quals" );
    
    # test when a string is given as a param
    throws_ok( sub{ $seq_obj->trim_back('a') },
              'MyX::Generic::Digit::MustBeDigit',
              "trim_back(a) - throws" );
    throws_ok( sub{ $seq_obj->trim_back('a') },
              qr/trim_back requires digit > 0/,
              "trim_back(a) - throws" );
    
    # test when a negative number is given as a param
    throws_ok( sub{ $seq_obj->trim_back(-1) },
              'MyX::Generic::Digit::TooSmall',
              "trim_back(-1) - throws" );
    throws_ok( sub{ $seq_obj->trim_back(-1) },
              qr/trim_back requires digit > 0/,
              "trim_back(-1) - throws" );
    
    # test a good run of trim_back
    lives_ok( sub{ $trimmed_obj = $seq_obj->trim_back(2) },
             "trim_back(2) - lives" );
    is( $trimmed_obj->get_seq(), "CG", "trim_back(2) - check trimmed seq" );
    is( $trimmed_obj->get_quals_str(), "BB", "trim_back(2) - check qual seq" );
    is( $seq_obj->get_seq(), "ATCGAT", "trim_back(2) - check kept seq" );
    is( $seq_obj->get_quals_str(), "AAAABB", "trim_back(2) - check kept quals" );
}


# test substr method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $qual = "AAAABBBB";
    
    my $seq_obj = BioUtils::FastqSeq->new({
        header => $header,
        seq => $seq,
        quals_str => $qual,
    });
    
    my $substr_obj;
    
    # test when no parameters are given
    throws_ok( sub{ $seq_obj->substr() },
              'MyX::Generic::Undef::Param',
              "substr() - throws" );
    throws_ok( sub{ $seq_obj->substr() },
              qr/start parameter not defined/,
              "substr() - throws" );
    
    # test when no end parameter is given
    throws_ok( sub{ $seq_obj->substr(1) },
              'MyX::Generic::Undef::Param',
              "substr(1) - throws" );
    throws_ok( sub{ $seq_obj->substr(1) },
              qr/end parameter not defined/,
              "substr(1) - throws" );
    
    # test when a non-number is given as start
    throws_ok( sub{ $seq_obj->substr("a", 1) },
              'MyX::Generic::Digit::MustBeDigit',
              "substr(a) - throws" );
    throws_ok( sub{ $seq_obj->substr("a", 1) },
              qr/substr start requires digit > 0/,
              "substr(a) - throws" );
    
    # test when a non-number is given as end
    throws_ok( sub{ $seq_obj->substr(1, "a") },
              'MyX::Generic::Digit::MustBeDigit',
              "substr(1,a) - throws" );
    throws_ok( sub{ $seq_obj->substr(1, "a") },
              qr/substr end requires digit > 0/,
              "substr(1,a) - throws" );
    
    # test when a negative number is given as the start param
    throws_ok( sub{ $seq_obj->substr(-1, 2) },
              'MyX::Generic::Digit::TooSmall',
              "substr(-1,2) - throws" );
    throws_ok( sub{ $seq_obj->substr(-1, 2) },
              qr/start requires digit > 0/,
              "substr(-1,2) - throws" );
    
    # test when a negative number is given as the end param
    throws_ok( sub{ $seq_obj->substr(1, -2) },
              'MyX::Generic::Digit::TooSmall',
              "substr(1,-2) - throws" );
    throws_ok( sub{ $seq_obj->substr(1, -2) },
              qr/end requires digit > 0/,
              "substr(1,-2) - throws" );
    
    # test when the end is smaller than the start
    throws_ok( sub{ $seq_obj->substr(3,2) },
               'MyX::Generic::Digit::TooBig',
               "substr(3,2) - throws" );
    throws_ok( sub{ $seq_obj->substr(3,2) },
              qr/end is larger than start/,
              "substr(3,2) - throws" );
    
    # test a good run of trim_back
    lives_ok( sub{ $substr_obj = $seq_obj->substr(2,4) },
             "substr(2, 4) - lives" );
    is( $substr_obj->get_seq(), "CGA", "substr(2,4) - check substr seq" );
    is( $substr_obj->get_quals_str(), "AAB", "substr(2,4) - check substr quals" );
}


# test rev method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $qual = "AAAABBBB";
    
    my $seq_obj = BioUtils::FastqSeq->new({
        header => $header,
        seq => $seq,
        quals_str => $qual,
    });
    
    # test when no parameter is given
    lives_ok( sub{ $seq_obj->rev() },
              "rev() - lives" );
    is( $seq_obj->get_seq(), (scalar reverse $seq), "rev() - get_seq" );
    is( $seq_obj->get_quals_str(), (scalar reverse $qual),
       "rev() - get_quals_str" );
}

# test comp method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $seq_comp = "TAGCTAGC";
    my $qual = "AAAABBBB";
    
    my $seq_obj = BioUtils::FastqSeq->new({
        header => $header,
        seq => $seq,
        quals_str => $qual,
    });
    
    # test when no parameter is given
    lives_ok( sub{ $seq_obj->comp() },
              "comp() - lives" );
    is( $seq_obj->get_seq(), $seq_comp, "comp() - get_seq" );
    is( $seq_obj->get_quals_str(), $qual,
       "comp() - get_quals_str" );
}

# test rev_comp method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $seq_rev_comp = "CGATCGAT";
    my $qual = "AAAABBBB";
    
    my $seq_obj = BioUtils::FastqSeq->new({
        header => $header,
        seq => $seq,
        quals_str => $qual,
    });
    
    # test when no parameter is given
    lives_ok( sub{ $seq_obj->rev_comp() },
              "rev_comp() - lives" );
    is( $seq_obj->get_seq(), $seq_rev_comp, "rev_comp() - get_seq" );
    is( $seq_obj->get_quals_str(), (scalar reverse $qual),
       "rev_comp() - get_quals_str" );
}


