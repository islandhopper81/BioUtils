
use BioUtils::ConsensusBuilder::FastqColumn;
use Test::More tests => 59;
use BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int);

BEGIN { use_ok( 'BioUtils::ConsensusBuilder::FastqColumn'); }

### Create a FastqColumn object to use for testing
my $fc = BioUtils::ConsensusBuilder::FastqColumn->new();
is( $fc->getBaseCount('A'), 0, "init A count == 0" );
is( $fc->getBaseCount('T'), 0, "init T count == 0" );
is( $fc->getBaseCount('C'), 0, "init C count == 0" );
is( $fc->getBaseCount('G'), 0, "init G count == 0" );
is( $fc->getBaseCount('N'), 0, "init N count == 0" );
is( $fc->getBaseCount('-'), 0, "init - count == 0" );
is( $fc->getBaseQuals("A")->[0], undef, "init A quals == ''" );

# Test the int_to_illumina_1_8 and illumina_1_8_to_int methods so I can use them later
is( int_to_illumina_1_8(40), 'I', "int_to_illumina_1_8(40)" );
is( illumina_1_8_to_int('I'), 40, "illumina_1_8_to_int('I')" );

# Test the addBase function
# I think it is easier to think in number so I am using the encoding methods
$fc->addBase("A", int_to_illumina_1_8(10));
is( $fc->getBaseCount("A"), 1, "addBase(A, +) -- count goes up" );
is( $fc->getBaseQuals("A")->[0], int_to_illumina_1_8(10), "addBase(A, +) -- quals are stored" );

$fc->addBase('-', int_to_illumina_1_8(0));
is( $fc->getBaseCount("-"), 1, "addBase(-, !) -- count goes up" );
is( $fc->getBaseQuals("-")->[0], int_to_illumina_1_8(0), "addBase(-, !) -- quals are stored" );


# Test the _getMeanQual function on the A base
is( $fc->_getMeanQual("A"), int_to_illumina_1_8(10), "_getMeanQual(A) -- one element (10)" );
$fc->addBase("A", int_to_illumina_1_8(20));
is( $fc->_getMeanQual("A"), int_to_illumina_1_8(15), "_getMeanQual(A) -- two elements (10, 20)" );
$fc->addBase("A", int_to_illumina_1_8(20));
is( $fc->_getMeanQual("A"), int_to_illumina_1_8(16), "_getMeanQual(A) -- two elements (10, 20, 20)" );
is( $fc->_getMeanQual("-"), int_to_illumina_1_8(0), "_getMeanQual(-) -- one element (0)" );
is( $fc->_getMeanQual("A-"), int_to_illumina_1_8(12), "_getMeanQual(A-) -- two bases (10, 20, 20, 0)" );


# Test the _getSortedTopBases function
is( $fc->_getSortedTopBases("TG"), "GT", "_getSortedTopBases(TG)" );
is( $fc->_getSortedTopBases("AG"), "AG", "_getSortedTopBases(AG)" );
is( $fc->_getSortedTopBases("ca"), "AC", "_getSortedTopBases(CA)" );
is( $fc->_getSortedTopBases("TGCA"), "ACGT", "_getSortedTopBases(TGCA)" );


# Test the getIupacCode function
is( $fc->getIupacCode("A"), "A", "getIupacCode(A)" );
is( $fc->getIupacCode("T"), "T", "getIupacCode(T)" );
is( $fc->getIupacCode("C"), "C", "getIupacCode(C)" );
is( $fc->getIupacCode("G"), "G", "getIupacCode(G)" );
is( $fc->getIupacCode("AG"), "R", "getIupacCode(AG)" );
is( $fc->getIupacCode("CT"), "Y", "getIupacCode(CT)" );
is( $fc->getIupacCode("GC"), "S", "getIupacCode(GC)" );
is( $fc->getIupacCode("AT"), "W", "getIupacCode(AT)" );
is( $fc->getIupacCode("GT"), "K", "getIupacCode(GT)" );
is( $fc->getIupacCode("AC"), "M", "getIupacCode(AC)" );
is( $fc->getIupacCode("AGT"), "D", "getIupacCode(AGT)" );
is( $fc->getIupacCode("ACT"), "H", "getIupacCode(ACT)" );
is( $fc->getIupacCode("ACG"), "V", "getIupacCode(ACG)" );
is( $fc->getIupacCode("CGT"), "B", "getIupacCode(CGT)" );


# Test getConBaseAndQual() function
my ($topBases, $meanQual) = $fc->getConBaseAndQual();
is( $topBases, "A", "getConBaseAndQual() -- topBases" );
is( $meanQual, int_to_illumina_1_8(16), "getConBaseAndQual() -- meanQual" );

# Now I am going to create a tie to resolve
$fc->addBase("T", int_to_illumina_1_8(10));
$fc->addBase("T", int_to_illumina_1_8(20));
$fc->addBase("T", int_to_illumina_1_8(10));
# Now I have a tie between the total counts of A and T BUT A has higher quals
is( $fc->_resolveTie("TA"), "A", "_resolveTie(TA)" );
($topBases, $meanQual) = $fc->getConBaseAndQual();
is( $topBases, "A", "getConBaseAndQual() -- topBases -- tie with unequal qual means" );
is( $meanQual, int_to_illumina_1_8(16), "getConBaseAndQual() -- meanQual -- tie with unequal qual menas" );

$fc->addBase("C", int_to_illumina_1_8(10));
$fc->addBase("C", int_to_illumina_1_8(20));
$fc->addBase("C", int_to_illumina_1_8(20));
# Now I have a three way tie between A, T, and C.  A and C have equal mean quals which are also greater than T
is( $fc->_resolveTie("TCA"), "AC", "_resolveTie(TCA)" );
($topBases, $meanQual) = $fc->getConBaseAndQual();
is( $topBases, "M", "getConBaseAndQual() -- topBases -- three way tie with top two equal qual means" );
is( $meanQual, int_to_illumina_1_8(16), "getConBaseAndQual() -- meanQual -- three way tie with top two equal qual means" );


# Now I am adding in a four way tie with a dash.  I should get the dash back as the topBase.
$fc->addBase("-", int_to_illumina_1_8(0));
$fc->addBase("-", int_to_illumina_1_8(0));
is( $fc->_resolveTie("TCA-"), "-", "_resolveTie(TCA-)" );
($topBases, $meanQual) = $fc->getConBaseAndQual();
is( $topBases, '-', "getConBaseAndQual() -- topBases -- four way tie with a dash" );
is( $meanQual, int_to_illumina_1_8(0), "getConBaseAndQual() -- meanQual -- four way tie with a dash" );



#### Create a new FastqColumn for testing.  I specifically want to test a dot base value here.
$fc = BioUtils::ConsensusBuilder::FastqColumn->new();
$fc->addBase("A", int_to_illumina_1_8(40));
$fc->addBase("A", int_to_illumina_1_8(30));
$fc->addBase("T", int_to_illumina_1_8(40));
$fc->addBase(".", int_to_illumina_1_8(0));
( $topBase, $meanQual ) = $fc->getConBaseAndQual();
is( $topBase, "A", "getConBaseAndQual() -- topBases -- with a dot");
is( $meanQual, int_to_illumina_1_8(35), "getConBaseAndQual() -- meanQual -- with a dot");



# test _qual_to_prob
{
    is( BioUtils::ConsensusBuilder::FastqColumn::_qual_to_prob(40), .0001, "_qual_to_prob(40)" );
    is( BioUtils::ConsensusBuilder::FastqColumn::_qual_to_prob(30), .001, "_qual_to_prob(30)" );
    is( BioUtils::ConsensusBuilder::FastqColumn::_qual_to_prob(20), .01, "_qual_to_prob(20)" );
}

# test _base_probability
{
    is( BioUtils::ConsensusBuilder::FastqColumn::_base_probability("T", "T", 40), .9999,
       "_base_probability(T, T, 40)" );
    is( BioUtils::ConsensusBuilder::FastqColumn::_base_probability("A", "T", 40), 3.33333333333333e-05,
       "_base_probability(A, T, 40" );
    is( BioUtils::ConsensusBuilder::FastqColumn::_base_probability("C", "T", 40), 3.33333333333333e-05,
       "_base_probability(A, T, 40" );
    is( BioUtils::ConsensusBuilder::FastqColumn::_base_probability("G", "T", 40), 3.33333333333333e-05,
       "_base_probability(A, T, 40" );
}

# test _by_probability
{
    $fc = BioUtils::ConsensusBuilder::FastqColumn->new();
    $fc->addBase("A", int_to_illumina_1_8(40));
    ( $topBase, $qual ) = $fc->_by_probability();
    is( $topBase, "A", "_by_probability base => A - 40" );
    is( $qual, .9999, "_by_probability qual => A - 40" );
}
