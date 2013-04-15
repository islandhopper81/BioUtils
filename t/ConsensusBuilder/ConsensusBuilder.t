
use strict;
use warnings;

use BioUtils::ConsensusBuilder::ConsensusBuilder;
use BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int );
use Data::Dumper qw( Dumper );
use File::Temp qw{tempfile};

use Test::More tests => 6;

BEGIN { use_ok( 'BioUtils::ConsensusBuilder::ConsensusBuilder'); }

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

my $consensus = BioUtils::ConsensusBuilder::ConsensusBuilder::buildFromClustalwFile($file1, \%qualsHash);
is( $consensus->get_seq(), "AATTCCGG", "buildConFromFile() -- AATTCCGG" );
is( $consensus->get_quals_aref->[1], int_to_illumina_1_8(40), "buildConFromFile() -- conQuals[0]" );


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

$consensus = BioUtils::ConsensusBuilder::ConsensusBuilder::buildFromClustalwFile($file2, \%qualsHash);
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

$consensus = BioUtils::ConsensusBuilder::ConsensusBuilder::buildFromClustalwFile($file3, \%qualsHash);
is( $consensus->get_seq(), "$seq1", "buildConFromFile() -- long sequence" );
#is( $consensus->get_quals_aref->[0], int_to_illumina_1_8(36), "buildConFromFile() -- conQuals[0]" );

