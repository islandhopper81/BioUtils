
use strict;
use warnings;

use Test::More tests => 46;
use Test::Exception;

# Test Module loading
BEGIN { use_ok( 'BioUtils::Codec::QualityScores',
               qw(int_to_illumina_1_8 illumina_1_8_to_int) ); }
can_ok( __PACKAGE__, 'int_to_illumina_1_8');
can_ok( __PACKAGE__, 'illumina_1_8_to_int');

# Test module functionality
{
    # Test going from int to illumina-1.8 encoding
    lives_ok( sub { int_to_illumina_1_8(0) }, 'Expected to live' );
    is( int_to_illumina_1_8(0), '!', '0 == !' );
    is( int_to_illumina_1_8(41), 'J', '41 == J' );
    dies_ok( sub { int_to_illumina_1_8(-1) }, 'Expect to die' );
    throws_ok( sub { int_to_illumina_1_8(-1) }, 'MyX::Generic::OutOfBounds',
              'Caught MyX::Generic::OutOfBounds' );
    throws_ok( sub { int_to_illumina_1_8(-1) }, qr/Int out of ASCII bounds/,
              'Int out of ASCII bounds' );
    dies_ok( sub { int_to_illumina_1_8(94) }, 'Expected to die' );
    throws_ok( sub { int_to_illumina_1_8(94) }, 'MyX::Generic::OutOfBounds',
              'Caught MyX::Generic::OutOfBounds' );
    throws_ok( sub { int_to_illumina_1_8(94) }, qr/Int out of ASCII bounds/,
              'Int out of ASCII bounds' );
    is( int_to_illumina_1_8([40, 39, 38]), 'IHG', '[40,39,38] == IHG' );
    
    # Test going from illumina-1.8 to int
    lives_ok( sub { illumina_1_8_to_int('I') }, 'Expected to live' );
    is( illumina_1_8_to_int('!'), 0, '! == 0' );
    is( illumina_1_8_to_int('J'), 41, 'J == 41' );
    dies_ok( sub { illumina_1_8_to_int(' ') }, 'Expect to die' );
    throws_ok( sub { illumina_1_8_to_int(' ') }, 'MyX::Generic::OutOfBounds',
              'Caught MyX::Generic::OutOfBounds' );
    throws_ok( sub { illumina_1_8_to_int(' ') }, qr/Char out of ASCII bounds/,
              'Char out of ASCII bounds' );
    is_deeply( illumina_1_8_to_int('IHG'), [40,39,38], 'IHG == [40,39,38]' );
}

# Test the Object functionality
{
    # test object creation (ie new)
    my $codec;
    lives_ok( sub{ $codec = BioUtils::Codec::QualityScores->new() }, 'Expected to live' );
    
    # test _set_int_to_illumina_1_8_lookup
    my $lookup_href;
    lives_ok( sub{ $lookup_href = $codec->_set_int_to_illumina_1_8_lookup() },
             'Expected to live' );
    is( $lookup_href->{0}, '!', '_set_int_to_illumina_1_8_lookup: 0 => !' );
    is( $lookup_href->{93}, '~', '_set_int_to_illumina_1_8_lookup: 93 => ~' );
    is( $lookup_href->{-1}, undef, '_set_int_to_illumina_1_8_lookup: -1 => undef' );
    is( $lookup_href->{94}, undef, '_set_int_to_illumina_1_8_lookup: 94 => undef' );
    
    # test _set_illumina_1_8_to_int_lookup
    $lookup_href = ();
    lives_ok( sub{ $lookup_href = $codec->_set_illumina_1_8_to_int_lookup() },
             'Expected to live' );
    is( $lookup_href->{'!'}, 0, '_set_int_to_illumina_1_8_lookup: ! => 0' );
    is( $lookup_href->{'~'}, 93, '_set_int_to_illumina_1_8_lookup: ~ => 93' );
    is( $lookup_href->{' '}, undef, '_set_int_to_illumina_1_8_lookup:  => undef' );
    
    # test int_to_illumina_1_8_lookup
    is( $codec->int_to_illumina_1_8_lookup(0), '!', '0 == !' );
    is( $codec->int_to_illumina_1_8_lookup(41), 'J', '41 == J' );
    dies_ok( sub { $codec->int_to_illumina_1_8_lookup(-1) },
            'Expect to die' );
    throws_ok( sub { $codec->int_to_illumina_1_8_lookup(-1) },
              'MyX::Generic::OutOfBounds',
              'Caught MyX::Generic::OutOfBounds' );
    throws_ok( sub { $codec->int_to_illumina_1_8_lookup(-1) },
              qr/Int out of ASCII bounds/,
              'Int out of ASCII bounds' );
    dies_ok( sub { $codec->int_to_illumina_1_8_lookup(94) },
            'Expected to die' );
    throws_ok( sub { $codec->int_to_illumina_1_8_lookup(94) },
              'MyX::Generic::OutOfBounds',
              'Caught MyX::Generic::OutOfBounds' );
    throws_ok( sub { $codec->int_to_illumina_1_8_lookup(94) },
              qr/Int out of ASCII bounds/,
              'Int out of ASCII bounds' );
    is( $codec->int_to_illumina_1_8_lookup([40, 39, 38]),
       'IHG',
       '[40,39,38] == IHG' );
    
    # test illumina_1_8_to_int_lookup
    lives_ok( sub { $codec->illumina_1_8_to_int_lookup('I') },
             'Expected to live' );
    is( $codec->illumina_1_8_to_int_lookup('!'), 0, '! == 0' );
    is( $codec->illumina_1_8_to_int_lookup('J'), 41, 'J == 41' );
    dies_ok( sub { $codec->illumina_1_8_to_int_lookup(' ') },
            'Expect to die' );
    throws_ok( sub { $codec->illumina_1_8_to_int_lookup(' ') },
              'MyX::Generic::OutOfBounds',
              'Caught MyX::Generic::OutOfBounds' );
    throws_ok( sub { $codec->illumina_1_8_to_int_lookup(' ') },
              qr/Char out of ASCII bounds/,
              'Char out of ASCII bounds' );
    is_deeply( $codec->illumina_1_8_to_int_lookup('IHG'),
              [40,39,38],
              'IHG == [40,39,38]' );
}
