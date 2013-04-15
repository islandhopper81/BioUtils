use Test::More tests => 11;

BEGIN {
use_ok( 'BioUtils::ConsensusBuilder::ConsensusBuilder' );
use_ok( 'BioUtils::ConsensusBuilder::FastqColumn' );
use_ok( 'BioUtils::ConsensusBuilder::FastqConsensus' );
use_ok( 'BioUtils::Codec::QualityScores', qw(int_to_illumina_1_8 illumina_1_8_to_int) );
use_ok( 'BioUtils::FastaSeq' );
use_ok( 'BioUtils::FastaIO' );
use_ok( 'BioUtils::FastqSeq' );
use_ok( 'BioUtils::FastqIO' );
use_ok( 'BioUtils::QC::FastqFilter' );
use_ok( 'BioUtils::QC::ContaminantFilter' );
use_ok( 'BioUtils::SeqSet::Diagnostics' );
}

diag( "\n" );
diag( "Testing BioUtils::ConsensusBuilder::ConsensusBuilder $BioUtils::ConsensusBuilder::ConsensusBuilder::VERSION" );
diag( "Testing BioUtils::ConsensusBuilder::FastqColumn $BioUtils::ConsensusBuilder::FastqColumn::VERSION" );
diag( "Testing BioUtils::ConsensusBuilder::FastqConsensus $BioUtils::ConsensusBuilder::FastqConsensus::VERSION" );
diag( "Testing BioUtils::Codec::QualityScores $BioUtils::Codec::QualityScores::VERSION" );
diag( "Testing BioUtils::FastaSeq $BioUtils::FastaSeq::VERSION" );
diag( "Testing BioUtils::FastaIO $BioUtils::FastaIO::VERSION" );
diag( "Testing BioUtils::FastqSeq $BioUtils::FastqSeq::VERSION" );
diag( "Testing BioUtils::FastqIO $BioUtils::FastqIO::VERSION" );
diag( "Testing BioUtils::QC::FastqFilter $BioUtils::QC::FastqFilter::VERSION" );
diag( "Testing BioUtils::QC::ContaminantFilter $BioUtils::QC::ContaminantFilter::VERSION" );
diag( "Testing BioUtils::SeqSet::Diagnostics $BioUtils::SeqSet::Diagnostics::VERSION" );
