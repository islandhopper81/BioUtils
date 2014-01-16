
use strict;
use warnings;

use BioUtils::ConsensusBuilder::ConsensusBuilder;
use BioUtils::FastqSeq;
use BioUtils::FastaSeq;
use BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int );
use Data::Dumper qw( Dumper );
use File::Temp qw{tempfile};

use Test::More tests => 34;
use Test::Exception;

# The use_ok does not work with exported methods
BEGIN { use_ok( 'BioUtils::ConsensusBuilder::ConsensusBuilder'); }

# form some reason use_ok doesn't work with exported methods so I added the
# following line.
use BioUtils::ConsensusBuilder::ConsensusBuilder qw(/.*/);

### Create a temporary file representing a clustalw output file
my ($fh1, $file1) = tempfile();
print $fh1 "#Seq1
AATTCCGG
#Seq2
AATTCCGG
#Seq3
AATTCCGG";
close($fh1);  # I have to close the fh because the Consensus::Builder also opens a fh.

my %qualsHash = ();
my @seq1 = map ( int_to_illumina_1_8($_), (40,40,40,40,40,40,40,40) );
my @seq2 = map ( int_to_illumina_1_8($_), (40,40,40,40,40,40,40,40) );
my @seq3 = map ( int_to_illumina_1_8($_), (40,40,40,40,40,40,40,40) );

$qualsHash{'Seq1'} = \@seq1;
$qualsHash{'Seq2'} = \@seq2;
$qualsHash{'Seq3'} = \@seq3;

my $consensus = buildFromClustalwFile($file1, \%qualsHash);
is( $consensus->get_seq(), "AATTCCGG", "buildConFromFile() -- AATTCCGG" );
is( $consensus->get_quals_aref->[1], int_to_illumina_1_8(40), "buildConFromFile() -- conQuals[0]" );
is( $consensus->get_c_score(), (40*8)/8, "buildConFromFile() -- get_c_score()" );

### Create a second temporary file representing a clustalw output file
my ($fh2, $file2) = tempfile();
print $fh2 "#Seq1
TATTCCGG
#Seq2
AAT-CCGG
#Seq3
AATTCCGT
#Seq4
AATTCCGT";
close($fh2);

%qualsHash = ();
@seq1 = map ( int_to_illumina_1_8($_), (40,40,40,40,40,40,40,40) );
@seq2 = map ( int_to_illumina_1_8($_), (40,40,40,40,40,40,40,40) );
@seq3 = map ( int_to_illumina_1_8($_), (30,40,40,40,40,40,40,40) );
my @seq4 = map ( int_to_illumina_1_8($_), (40,40,40,40,40,40,40,40) );

$qualsHash{'Seq1'} = \@seq1;
$qualsHash{'Seq2'} = \@seq2;
$qualsHash{'Seq3'} = \@seq3;
$qualsHash{'Seq4'} = \@seq4;

$consensus = buildFromClustalwFile($file2, \%qualsHash);
is( $consensus->get_seq(), "AATTCCGK", "buildConFromFile() -- AATTCCGG" );
is( $consensus->get_quals_aref->[0], int_to_illumina_1_8(36), "buildConFromFile() -- conQuals[0]" );

### Test this pair of full length sequences using the clustalw approach.
# They were giving me problems.
# Hmmm... it seems to work here.  The problem must be somewhere in MT_Toolbox
my ($fh3, $file3) = tempfile();

my $seq1 = 'GTGCCAGCAGCCGCGGTAA' .
           'TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA' .
           'GGACTACCGGGGTTTCTAAT';
my $seq2 = "GTGCCAGCAGCCGCGGTAA" . 
           "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA" .
           "GGACTACCGGGGTTTCTAAT";

print $fh3 "#1_1_13352_2073
$seq1
#1_1_13352_2074
$seq2";
close($fh3);

my $quals_str1 = "IJJIJGIIGGHIGGIIHID" .
                 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#' .
                 'DDFDHHHBHIIIHIIIIIII';
my $quals_str2 = "IJJIJGIIGGHIGGIIHID" .
                 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#' .
                 'DDFDHHHBHIIIHIIIIIII';

my @quals_arr1 = split //, $quals_str1;
my @quals_arr2 = split //, $quals_str2;
%qualsHash = ();
$qualsHash{'1_1_13352_2073'} = \@quals_arr1;
$qualsHash{'1_1_13352_2074'} = \@quals_arr2;

$consensus = buildFromClustalwFile($file3, \%qualsHash);
is( $consensus->get_seq(), "$seq1", "buildConFromFile() -- long sequence" );
#is( $consensus->get_quals_aref->[0], int_to_illumina_1_8(36), "buildConFromFile() -- conQuals[0]" );

###################
# I have written some new methods to build a consensus sequence from FastQ objs.
# The tests are below:

# build 2 fastqSeq objects for testing
my $s1_obj = BioUtils::FastqSeq->new({
                header => "seq1",
                seq => "ATCG",
                quals_str => "DDDD",
                });

my $s2_obj = BioUtils::FastqSeq->new({
                header => "seq2",
                seq => "ATCT",
                quals_str => "DDDC",
                });
my @seqs_arr = ($s1_obj, $s2_obj);

# test build_consensus_from_aref
{
    my @seq_arr = ($s1_obj, $s2_obj);
    my $con;
    lives_ok( sub{ $con = build_consensus( \@seqs_arr ) },
             "build_consensus(aref) -- lives" );
    is( $con->get_seq(), "ATCG", "testing con base accuracy" );
    is( $con->get_quals_str(), "DDDD", "testing con quals accuracy" );
}

# test when seqs are not square
{
    $s1_obj->set_seq("A");
    $s1_obj->set_quals_str("D");
    throws_ok( sub{ build_consensus(\@seqs_arr) },
              "BioUtils::MyX::ConsensusBuilder::SeqsNotSqr",
              "build_consensus(SeqsNotSqr) - caught" );
    throws_ok( sub{ build_consensus(\@seqs_arr) },
              qr/Alignment seqs are not square/,
              "build_consensus(SeqsNotSqr) - caught" );
}

# test when quals are not square
{
    $s1_obj->set_seq("ATCG");
    $s1_obj->set_quals_str("DD");
    throws_ok( sub{ build_consensus(\@seqs_arr) },
              "BioUtils::MyX::ConsensusBuilder::QualsNotSqr",
              "build_consensus(QualsNotSqr) - caught" );
    throws_ok( sub{ build_consensus(\@seqs_arr) },
              qr/Alignment quals are not square/,
              "build_consensus(QualsNotSqr) - caught" );
}

# test when a hash ref is used
{
    # reset seq1 to something that works
    $s1_obj->set_seq("ATCG");
    $s1_obj->set_quals_str("DDDD");
    my %seq_hash = ("seq1" => $s1_obj, "seq2" => $s2_obj);
    
    my $con;
    lives_ok( sub{ $con = build_consensus( \%seq_hash ) },
             "build_consensus(href) -- lives" );
    is( $con->get_seq(), "ATCG", "testing con base accuracy" );
    is( $con->get_quals_str(), "DDDD", "testing con quals accuracy" );
}

# test _check_seq_ref
{
    # create some object for testing
    my $fastq_obj = BioUtils::FastqSeq->new({
                header => "seq1",
                seq => "ATCG",
                quals_str => "DDDD",
                });
    my $fasta_obj = BioUtils::FastaSeq->new({
                header => "seq1",
                seq => "ATCG",
                });
    my $int = 0;
    my $str = "a";
    my $href = {};
    
    is( BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($fastq_obj),
       1, "_check_seq_ref(fastq_obj) -- lives" );
    is( BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($fasta_obj),
       1, "_check_seq_ref(fasta_obj) -- lives" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($int) },
       "MyX::Generic::Ref::UnsupportedType",
       "_check_seq_ref(int) -- throws" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($int) },
       qr/Reference must be BioUtils::FastaSeq or BioUtils::FastqSeq/,
       "_check_seq_ref(int) -- throws" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($str) },
       "MyX::Generic::Ref::UnsupportedType",
       "_check_seq_ref(str) -- throws" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($str) },
       qr/Reference must be BioUtils::FastaSeq or BioUtils::FastqSeq/,
       "_check_seq_ref(str) -- throws" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($href) },
       "MyX::Generic::Ref::UnsupportedType",
       "_check_seq_ref(href) -- throws" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_ref($href) },
       qr/Reference must be BioUtils::FastaSeq or BioUtils::FastqSeq/,
       "_check_seq_ref(href) -- throws" );
}

# test when I pass a bad reference type to build_consensus
{
    # some bad types to test.  Good ones worked above -- href, aref
    my $int = 0;
    my $str = "a";
    my $href = {};
    
    throws_ok( sub{ build_consensus($int) },
       "MyX::Generic::Ref::UnsupportedType",
       "build_consensus(int) -- throws" );
    throws_ok( sub{ build_consensus($int) },
       qr/build_consensus expects HASH or ARRAY reference/,
       "build_consensus(int) -- throws" );
    throws_ok( sub{ build_consensus($str) },
       "MyX::Generic::Ref::UnsupportedType",
       "build_consensus(str) -- throws" );
    throws_ok( sub{ build_consensus($str) },
       qr/build_consensus expects HASH or ARRAY reference/,
       "build_consensus(str) -- throws" );
}

# test _check_seq_count
{
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_count(0) },
       "BioUtils::MyX::ConsensusBuilder::NoSeqs",
       "_check_seq_count(0) - throws" );
    throws_ok( sub{ BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_count(0) },
       qr/No seqs to build a consensus from/,
       "_check_seq_count(0) - throws" );
    
    is( BioUtils::ConsensusBuilder::ConsensusBuilder::_check_seq_count(1),
       1, "_check_seq_count(1)" );
    
    my $href = {};
    throws_ok( sub{ build_consensus($href) },
       "BioUtils::MyX::ConsensusBuilder::NoSeqs",
       "build_consensus(empty href) -- throws" );
    throws_ok( sub{ build_consensus($href) },
       qr/No seqs to build a consensus from/,
       "build_consensus(empty href) -- throws" );
}

