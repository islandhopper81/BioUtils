
use BioUtils::SeqSet::Diagnostics;
use Test::More tests => 70;
use Test::Exception;
use Test::Warn;

BEGIN { use_ok( 'BioUtils::SeqSet::Diagnostics'); }

# Create a SeqienceDiagnostics object to test
my $my_seq_diag = BioUtils::SeqSet::Diagnostics->new();

# Test the gap character stuff
{
    # Test defualt settings
    is( $my_seq_diag->has_gap_character("-"), 1, "has_gap_character(-) - double quotes" );
    is( $my_seq_diag->has_gap_character('-'), 1, "has_gap_character(-) - single quotes" );
    is( $my_seq_diag->has_gap_character('#'), 0, "has_gap_character(#) - no" );
    is( $my_seq_diag->has_gap_character('.'), 0, "has_gap_character(.) - no" );
    
    # Test has_gap_character -- some of the above tests apply here too.
    dies_ok( sub { $my_seq_diag->has_gap_character("12") },
            "expected to die -- has_gap_character(12)" );
    dies_ok( sub { $my_seq_diag->has_gap_character(" ") },
            "expected to die -- has_gap_character( )" );
    
    # Test add_gap_character
    is( $my_seq_diag->has_gap_character('#'), 0, "has_gap_character(#) - no");
    is( $my_seq_diag->add_gap_character('#'), 1, "add_gap_character(#)" );
    is( $my_seq_diag->has_gap_character('#'), 1, "has_gap_character(#) - yes" );
    dies_ok( sub { $my_seq_diag->add_gap_character("123") },
            "expected to die -- add_gap_character(123)" );
    dies_ok( sub { $my_seq_diag->add_gap_character(" ") },
            "expected to die -- add_gap_character( )" );
    
    # Test remove_gap_character
    is( $my_seq_diag->has_gap_character('-'), 1, "has_gap_character(-) - yes" );
    is( $my_seq_diag->remove_gap_character('-'), 1, "remove_gap_character(-)" );
    is( $my_seq_diag->has_gap_character('-'), 0, "has_gap_character(-) - no" );
    dies_ok( sub { $my_seq_diag->remove_gap_character("123") },
            "expected to die -- remove_gap_character(123)" );
    dies_ok( sub { $my_seq_diag->remove_gap_character(" ") },
            "expected to die -- remove_gap_character( )" );
    
    # Remove all gap characters
    warning_is( sub { $my_seq_diag->remove_gap_character('#') },
               "All gap characters have been removed",
               "remove_gap_character(#)"
               );
    
    # Add back a gap character
    warning_is( sub { $my_seq_diag->add_gap_character('-') },
               "Gap character no longer empty after adding: -",
               "add_gap_character(-)"
               );
}

# Test the initialized attributes
is( $my_seq_diag->get_gap_count(), 0, "Initial gap count" );
is_deeply( $my_seq_diag->get_gap_pos_dist(), {}, "Initial gap pos dist" );
is_deeply( $my_seq_diag->get_gaps_per_seq(), {}, "Initial gaps per seq" );
is( $my_seq_diag->get_ambig_base_count() , 0, "Initial ambiguous base count" );
is_deeply( $my_seq_diag->get_ambig_base_pos_dist(), {},
   "Initial ambiguous base pos dist" );
is_deeply( $my_seq_diag->get_ambig_base_comp_dist(), {},
   "Initial ambigous base composition dist" );
is_deeply( $my_seq_diag->get_ambig_bases_per_seq(), {},
    "Initial gaps per seq" );


# Test _add_gap_at_pos
{
    $my_seq_diag->_init(); # reset the attribute values
    # add one gap
    $my_seq_diag->_add_gap_at_pos(4);
    is( $my_seq_diag->get_gap_count(), 1, "_add_gap_at_pos(4) - gap_count increment" );
    is_deeply( $my_seq_diag->get_gap_pos_dist(), {4 => 1}, "_add_gap_at_pos(4) - {4 => 1}" );
    
    # add a second gap at the same position
    $my_seq_diag->_add_gap_at_pos(4);
    is( $my_seq_diag->get_gap_count(), 2, "_add_gap_at_pos(4) - gap_count increment" );
    is_deeply( $my_seq_diag->get_gap_pos_dist(), {4 => 2}, "_add_gap_at_pos(4) - {4 => 1}" );
    
    # add a third gap at a different position
    $my_seq_diag->_add_gap_at_pos(5);
    is( $my_seq_diag->get_gap_count(), 3, "_add_gap_at_pos(4) - gap_count increment" );
    is_deeply( $my_seq_diag->get_gap_pos_dist(), {4 => 2, 5 => 1}, "_add_gap_at_pos(4) - {4 => 1}" );
}


# Test _add_ambig_base_at_pos
{
    $my_seq_diag->_init();  # reset the attribute ablues
    
    # add one ambiguous base - N
    $my_seq_diag->_add_ambig_base_at_pos("N", 0);
    is( $my_seq_diag->get_ambig_base_count(), 1,
       "_add_ambig_base_at_pos() - ambig_base_count increment" );
    is_deeply( $my_seq_diag->get_ambig_base_pos_dist(), {0 => 1},
              "_add_ambig_base_at_pos(N, 0) - {0 => 1}" );
    is_deeply( $my_seq_diag->get_ambig_base_comp_dist(), {N => 1},
              "_add_ambig_base_at_pos(N, 0) - {N => 1}" );
    
    # add a second ambiguous base - Y
    $my_seq_diag->_add_ambig_base_at_pos("Y", 2);
    is( $my_seq_diag->get_ambig_base_count(), 2,
       "_add_ambig_base_at_pos() - ambig_base_count increment" );
    is_deeply( $my_seq_diag->get_ambig_base_pos_dist(), {0 => 1, 2 => 1},
              "_add_ambig_base_at_pos(N, 0) - {0 => 1, 2 => 1}" );
    is_deeply( $my_seq_diag->get_ambig_base_comp_dist(), {N => 1, Y => 1},
              "_add_ambig_base_at_pos(N, 0) - {N => 1, Y => 1}" );
    
    # add a third ambiguous base - Y (again)
    $my_seq_diag->_add_ambig_base_at_pos("Y", 2);
    is( $my_seq_diag->get_ambig_base_count(), 3,
       "_add_ambig_base_at_pos() - ambig_base_count increment" );
    is_deeply( $my_seq_diag->get_ambig_base_pos_dist(), {0 => 1, 2 => 2},
              "_add_ambig_base_at_pos(N, 0) - {0 => 1, 2 => 2}" );
    is_deeply( $my_seq_diag->get_ambig_base_comp_dist(), {N => 1, Y => 2},
              "_add_ambig_base_at_pos(N, 0) - {N => 1, Y => 2}" );
}


# Test diagnose_sequence
{
    $my_seq_diag->_init(); # reset the attribute values
    
    # Test with only one gap and no ambiguous bases
    $my_seq_diag->diagnose_sequence("ATCG-");
    is( $my_seq_diag->get_gap_count(), 1,
       "diagnose_sequence() - gap_count increment" );
    is_deeply( $my_seq_diag->get_gap_pos_dist(), {4 => 1},
              "diagnose_sequence() - {4 => 1}" );
    is_deeply( $my_seq_diag->get_gaps_per_seq(), {1 => 1},
              "diagnose_sequence() - {1 => 1}" );
    
    # Test a sequence with multiple gaps
    $my_seq_diag->diagnose_sequence("AT--T");
    is( $my_seq_diag->get_gap_count(), 3,
       "diagnose_sequence() - gap_count increment by 2" );
    is_deeply( $my_seq_diag->get_gap_pos_dist(), {2 => 1, 3 => 1, 4 => 1},
              "diagnose_sequence() - {2 => 1, 3 => 1, 4 => 1}" );
    is_deeply( $my_seq_diag->get_gaps_per_seq(), {1 => 1, 2 => 1},
              "diagnose_sequence() - {1 => 1, 2 => 1}" );
    
    # Test a sequence with an ambiguous base
    $my_seq_diag->diagnose_sequence("ATCGN");
    is ( $my_seq_diag->get_ambig_base_count(), 1,
        "diagnose_sequence(ATCGN) - ambig_base_count increment" );
    is_deeply( $my_seq_diag->get_ambig_base_pos_dist(), {4 => 1},
              "diagnose_sequence(ATCGN) - {4 => 1}" );
    is_deeply( $my_seq_diag->get_ambig_base_comp_dist(), {N => 1},
              "diagnose_sequence(ATCGN) - {N => 1}" );
    is_deeply( $my_seq_diag->get_ambig_bases_per_seq(), {1 => 1},
              "diagnose_sequence(ATCGN) - {1 => 1}" );
    
    # Test a sequnce with multiple anbiguous bases
    $my_seq_diag->diagnose_sequence("ATCBD");
    is ( $my_seq_diag->get_ambig_base_count(), 3,
        "diagnose_sequence(ATCBD) - ambig_base_count increment by 2" );
    is_deeply( $my_seq_diag->get_ambig_base_pos_dist(), {4 => 2, 3 => 1},
              "diagnose_sequence(ATCBD) - {4 => 1, 3 => 1}" );
    is_deeply( $my_seq_diag->get_ambig_base_comp_dist(), {N => 1, B => 1, D => 1},
              "diagnose_sequence(ATCBD) - {N => 1, B => 1, D => 1}" );
    is_deeply( $my_seq_diag->get_ambig_bases_per_seq(), {1 => 1, 2 => 1},
              "diagnose_sequence(ATCBD) - {1 => 1, 2 => 1}" );
}


# Test the distribution string getter methods
{
    $my_seq_diag->_init(); # reset the attribute values
    
    # Test get_gap_pos_dist_str by diagnosing the following 3 sequences
    $my_seq_diag->diagnose_sequence("AT-TG");
    $my_seq_diag->diagnose_sequence("-ATTG");
    $my_seq_diag->diagnose_sequence("AT-AT");
    my $dist_str = "0\t1\n2\t2\n";
    is( $my_seq_diag->get_gap_pos_dist_str(), $dist_str,
       "get_gap_pos_dist_str()" );
    
    # Test get_gaps_per_seq_str using the above sequences
    my $gaps_per_seq_str = "1\t3\n";
    is( $my_seq_diag->get_gaps_per_seq_str(), $gaps_per_seq_str,
       "get_gaps_per_seq_str()" );
    
    # Test get_ambig_base_pos_dist_str by diagnosing the following 3 sequences
    $my_seq_diag->diagnose_sequence("ATTRS");
    $my_seq_diag->diagnose_sequence("HAGGG");
    $my_seq_diag->diagnose_sequence("ATTCW");
    my $pos_dist_str = "0\t1\n3\t1\n4\t2\n";
    is( $my_seq_diag->get_ambig_base_pos_dist_str(), $pos_dist_str,
       "get_ambig_base_pos_dist_str()" );
    
    # Test get_ambig_bases_per_seq_str using the above sequences
    my $ambig_bases_per_seq_str = "1\t2\n2\t1\n";
    is( $my_seq_diag->get_ambig_bases_per_seq_str(), $ambig_bases_per_seq_str,
       "get_ambig_bases_per_seq_str()" );
    
    # Test get_ambig_base_comp_dist_str using the above sequences
    my $comp_dist_str = "H\t1\nR\t1\nS\t1\nW\t1\n";
    is( $my_seq_diag->get_ambig_base_comp_dist_str(), $comp_dist_str,
       "get_ambig_base_comp_dist_str()" );
}


# Test the graphing methods
{
    $my_seq_diag->_init(); # reset the attribute values
    
    # Diagnose some sequences
    $my_seq_diag->diagnose_sequence("AT-TG");
    $my_seq_diag->diagnose_sequence("-ATTG");
    $my_seq_diag->diagnose_sequence("AT-AT");
    $my_seq_diag->diagnose_sequence("ATTRS");
    $my_seq_diag->diagnose_sequence("HAGGG");
    $my_seq_diag->diagnose_sequence("ATTCW");
    
    lives_ok( sub { $my_seq_diag->graph_gap_pos_dist("temp_gap_pos_graph") },
             "Expected to live");
    lives_ok( sub { $my_seq_diag->graph_gaps_per_seq("temp_gaps_per_seq_graph") },
             "Expected to live");
    lives_ok( sub { $my_seq_diag->graph_ambig_base_pos_dist("temp_ambig_base_pos_graph") },
             "Expected to live");
    lives_ok( sub { $my_seq_diag->graph_ambig_base_comp_dist("temp_ambig_base_comp_graph") },
             "Expected to live");
    lives_ok( sub { $my_seq_diag->graph_ambig_bases_per_seq("temp_ambig_bases_per_seq_graph") },
             "Expected to live");
    
    is( -f "temp_gap_pos_graph.png", 1, "temp_gap_pos_graph made" );
    is( -f "temp_gaps_per_seq_graph.png", 1, "temp_gaps_per_seq_graph made" );
    is( -f "temp_ambig_base_pos_graph.png", 1, "temp_ambig_base_pos_graph made" );
    is( -f "temp_ambig_base_comp_graph.png", 1, "temp_ambig_base_comp_graph made" );
    is( -f "temp_ambig_bases_per_seq_graph.png", 1, "temp_ambig_bases_per_seq_graph made");
    
    # remove the temporary graphs that I made
    `rm temp_gap_pos_graph.png`;
    `rm temp_gaps_per_seq_graph.png`;
    `rm temp_ambig_base_pos_graph.png`;
    `rm temp_ambig_base_comp_graph.png`;
    `rm temp_ambig_bases_per_seq_graph.png`;
}

