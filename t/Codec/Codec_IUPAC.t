
use strict;
use warnings;

use Test::More tests => 17;
use Test::Exception;

# Test Module loading
BEGIN { use_ok( 'BioUtils::Codec::IUPAC',
               qw(nuc_str_to_iupac iupac_to_nuc_str) ); }
can_ok( __PACKAGE__, 'nuc_str_to_iupac');
can_ok( __PACKAGE__, 'iupac_to_nuc_str');

# Test module functionality
{
    # Test going from nucleotide string to IUPAC
    lives_ok( sub { nuc_str_to_iupac('A') }, 'Expected to live' );
    is( nuc_str_to_iupac('A'), 'A', 'A == A' );
    is( nuc_str_to_iupac("ACG"), 'V', 'ACG == V' );
    dies_ok( sub { nuc_str_to_iupac("X") }, 'Expect to die' );
    throws_ok( sub { nuc_str_to_iupac("X") }, 'MyX::Generic::BadValue',
              'Caught MyX::Generic::BadValue' );
    throws_ok( sub { nuc_str_to_iupac("X") }, qr/Unrecignized nucleotide base/,
              'Unrecignized nucleotide base' );
    throws_ok( sub { nuc_str_to_iupac() }, 'MyX::Generic::Undef::Param',
              'Caught MyX::Generic::Undef::Param' );
    
    # Test going from IUPAC to nucleotide string
    lives_ok( sub { iupac_to_nuc_str('T') }, 'Expected to live' );
    is( iupac_to_nuc_str('T'), 'T', 'T == T' );
    is( iupac_to_nuc_str('H'), 'ACT', 'ACT == ACT' );
    dies_ok( sub { iupac_to_nuc_str('Z') }, 'Expect to die' );
    throws_ok( sub { iupac_to_nuc_str('Z') }, 'MyX::Generic::BadValue',
              'Caught MyX::Generic::BadValue' );
    throws_ok( sub { iupac_to_nuc_str('Z') }, qr/Not an IUPAC coded character/,
              'Not an IUPAC coded character' );
    throws_ok( sub { iupac_to_nuc_str() }, 'MyX::Generic::Undef::Param',
              'Caught MyX::Generic::Undef::Param' );
}

