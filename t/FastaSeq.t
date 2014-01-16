
use strict;
use warnings;

use BioUtils::FastaSeq;
use Test::More tests => 42;
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





# test trim_front method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    
    my $seq_obj = BioUtils::FastaSeq->new({
        header => $header,
        seq => $seq,
    });
    
    my $trimmed_obj;
    
    # test when no parameter is given
    lives_ok( sub{ $trimmed_obj = $seq_obj->trim_front() },
              "trim_front() - throws" );
    is( $trimmed_obj->get_seq(), "", "trim_front(2) - check trimmed seq" );
    is( $seq_obj->get_seq(), "ATCGATCG", "trim_front(2) - check kept seq" );
    
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
    is( $seq_obj->get_seq(), "CGATCG", "trim_front(2) - check kept seq" );
}

# test trim_back method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    
    my $seq_obj = BioUtils::FastaSeq->new({
        header => $header,
        seq => $seq,
    });
    
    my $trimmed_obj;
    
     # test when no parameter is given
    lives_ok( sub{ $trimmed_obj = $seq_obj->trim_back() },
              "trim_back() - throws" );
    is( $trimmed_obj->get_seq(), "", "trim_back(2) - check trimmed seq" );
    is( $seq_obj->get_seq(), "ATCGATCG", "trim_back(2) - check kept seq" );
    
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
    is( $seq_obj->get_seq(), "ATCGAT", "trim_back(2) - check kept seq" );
}


# test rev method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    
    my $seq_obj = BioUtils::FastaSeq->new({
        header => $header,
        seq => $seq,
    });
    
    # test when no parameter is given
    lives_ok( sub{ $seq_obj->rev() },
              "rev() - lives" );
    is( $seq_obj->get_seq(), (scalar reverse $seq), "rev() - get_seq" );
}

# test comp method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $seq_comp = "TAGCTAGC";
    
    my $seq_obj = BioUtils::FastaSeq->new({
        header => $header,
        seq => $seq,
    });
    
    # test when no parameter is given
    lives_ok( sub{ $seq_obj->comp() },
              "comp() - lives" );
    is( $seq_obj->get_seq(), $seq_comp, "comp() - get_seq" );
}

# test rev_comp method
{
    my $header = "seq1";
    my $seq = "ATCGATCG";
    my $seq_rev_comp = "CGATCGAT";
    
    my $seq_obj = BioUtils::FastaSeq->new({
        header => $header,
        seq => $seq,
    });
    
    # test when no parameter is given
    lives_ok( sub{ $seq_obj->rev_comp() },
              "rev_comp() - lives" );
    is( $seq_obj->get_seq(), $seq_rev_comp, "rev_comp() - get_seq" );
}



