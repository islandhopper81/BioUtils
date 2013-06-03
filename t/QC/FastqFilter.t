
use BioUtils::QC::FastqFilter 1.0.2;
use Test::More tests => 127;
use Test::Exception;
use Test::Warn;
use File::Temp qw/ tempfile tempdir /;

sub _build_fastq_file;
sub _build_fastq_file2;

BEGIN { use_ok( 'BioUtils::QC::FastqFilter'); }

# make a temp fastq file
my ($fh, $filename) = tempfile();
_build_fastq_file($fh);

# make a temp output dir
my $output_dir = tempdir();

# make a global filter
my $filter_g = BioUtils::QC::FastqFilter->new({
                    fastq_file => $filename,
                    output_dir => $output_dir
                });

# test set_fastq_file
{
    # test undef file
    throws_ok( sub{$filter_g->set_fastq_file(undef)},
              qr/Undefined fastq_file/, "set_fastq_file undef - throws");
    
    # test empty file
    my ($fh_new, $file_new) = tempfile();
    close($fh_new);
    warnings_exist {$filter_g->set_fastq_file($file_new)}
              [qr/Empty fastq_file/], "set_fastq_file empty - throws";
    
    # test good file
    lives_ok( sub{$filter_g->set_fastq_file($filename)},
             "set_fastq_file - lives");
}

# test set_output_dir
{
    # an undef output_dir
   throws_ok( sub{$filter_g->set_output_dir(undef)},
              qr/Undefined output dir/, "set_output_dir undef - throws");
              
    # a file for an output_dir
    throws_ok( sub{$filter_g->set_output_dir($filename)},
              qr/Output dir doesn't exist/,
              "set_output_dir as file - throws");
    
    # all good
    lives_ok( sub{$filter_g->set_output_dir($output_dir)},
             "set_output_dir - lives");
}

# test set_min_len
{
    # negative
    my $min_len = -1;
    throws_ok( sub{$filter_g->set_min_len($min_len)},
              qr/min_len must be > 0/,
              "set_min_len(-1) - throws");
    
    # all good
    $min_len = 100;
    lives_ok( sub{$filter_g->set_min_len($min_len)},
             "set_min_len(100) - lives");
}

# test set_min_avg_qual
{
    # negative
    my $min_avg_qual = -1;
    throws_ok( sub{$filter_g->set_min_avg_qual($min_avg_qual)},
              qr/min_avg_qual must be > 0/,
              "set_min_avg_qual(-1) - throws");
    
    # all good
    $min_avg_qual = 100;
    lives_ok( sub{$filter_g->set_min_avg_qual($min_avg_qual)},
             "set_min_avg_qual(100) - lives");
}

# test set_min_base_qual
{
    # negative
    my $min_base_qual = -1;
    throws_ok( sub{$filter_g->set_min_base_qual($min_base_qual)},
              qr/min_base_qual must be > 0/,
              "set_min_base_qual(-1) - throws");
    
    # all good
    $min_base_qual = 100;
    lives_ok( sub{$filter_g->set_min_base_qual($min_base_qual)},
             "set_min_base_qual(100) - lives");
}

# test _Y_N_translate -- I test this here because I use it in set_allow_gaps
# which is tested next
{
    is( BioUtils::QC::FastqFilter::_Y_N_translate('Y'), 1, "_Y_N_translate(Y)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('Yes'), 1, "_Y_N_translate(Yes)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('YES'), 1, "_Y_N_translate(YES)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('y'), 1, "_Y_N_translate(y)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('N'), 0, "_Y_N_translate(N)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('No'), 0, "_Y_N_translate(No)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('no'), 0, "_Y_N_translate(no)" );
    is( BioUtils::QC::FastqFilter::_Y_N_translate('n'), 0, "_Y_N_translate(n)" );
    
    throws_ok( sub{BioUtils::QC::FastqFilter::_Y_N_translate('Blah')},
              qr/Bad Yes\/No value/, "_Y_N_translate(Blah) - throws" );
}

# test set_allow_gaps
{
    # too big
    my $allow_gaps = 10;
    throws_ok( sub{$filter_g->set_allow_gaps($allow_gaps)},
              qr/allow_gaps must be 0 or 1/,
              "set_allow_gaps(10) - throws");
    
    # too small
    $allow_gaps = -10;
    throws_ok( sub{$filter_g->set_allow_gaps($allow_gaps)},
              qr/allow_gaps must be 0 or 1/,
              "set_allow_gaps(-10) - throws");
    
    # all good
    $allow_gaps = 1;
    lives_ok( sub{$filter_g->set_allow_gaps($allow_gaps)},
             "set_allow_gaps(1) - lives");
    
    # all good
    $allow_gaps = 0;
    lives_ok( sub{$filter_g->set_allow_gaps($allow_gaps)},
             "set_allow_gaps(0) - lives");
    
    # all good
    $allow_gaps = 'Y';
    lives_ok( sub{$filter_g->set_allow_gaps($allow_gaps)},
             "set_allow_gaps(Y) - lives");
    
    # all good
    $allow_gaps = 'N';
    lives_ok( sub{$filter_g->set_allow_gaps($allow_gaps)},
             "set_allow_gaps(N) - lives");
}

# test allow_ambig_bases
{
    # too big
    my $allow_ambig = 10;
    throws_ok( sub{$filter_g->set_allow_ambig_bases($allow_ambig)},
              qr/allow_ambig_bases must be 0 or 1/,
              "set_allow_ambig_bases(10) - throws");
    
    # too small
    $allow_ambig = -10;
    throws_ok( sub{$filter_g->set_allow_ambig_bases($allow_ambig)},
              qr/allow_ambig_bases must be 0 or 1/,
              "set_allow_ambig_bases(-10) - throws");
    
    # all good
    $allow_ambig = 1;
    lives_ok( sub{$filter_g->set_allow_ambig_bases($allow_ambig)},
             "set_allow_ambig_bases(1) - lives");
    
    # all good
    $allow_ambig = 0;
    lives_ok( sub{$filter_g->set_allow_ambig_bases($allow_ambig)},
             "set_allow_ambig_bases(0) - lives");
    
    # all good
    $allow_ambig = 'Y';
    lives_ok( sub{$filter_g->set_allow_ambig_bases($allow_ambig)},
             "set_allow_ambig_bases(Y) - lives");
    
    # all good
    $allow_ambig = 'N';
    lives_ok( sub{$filter_g->set_allow_ambig_bases($allow_ambig)},
             "set_allow_ambig_bases(N) - lives");
}

# test _init
{
    my $filter;
    
    # no fastq file given
    throws_ok( sub{$filter = BioUtils::QC::FastqFilter->new()},
              qr/Undefined fastq_file/, "new no fastq file - throws");
    
    # empty fastq file
    my ($fh_2, $filename_2) = tempfile();
    close($fh_2);
    warnings_exist {$filter = BioUtils::QC::FastqFilter->new({
                                                fastq_file => $filename_2,
                                                output_dir => $output_dir
                                                })}
        [qr/Empty fastq_file/,],
        "new empty fastq file - warnings";
    
    # good fastq file
    lives_ok( sub{$filter = BioUtils::QC::FastqFilter->new({fastq_file => $filename,
                                              output_dir => $output_dir})},
             "new - lives" );
    
    # no output dir
    throws_ok( sub{$filter = BioUtils::QC::FastqFilter->new({fastq_file => $filename})},
              qr/Undefined output dir/, "new no output dir- throws");
    
    # good fastq file and output dir
    lives_ok( sub{$filter = BioUtils::QC::FastqFilter->new({fastq_file => $filename,
                                              output_dir => $output_dir,
                                              })},
             "new - lives" );
}

# test all getters
{
    my $filter_l = BioUtils::QC::FastqFilter->new({fastq_file => $filename,
                                     output_dir => $output_dir,
                                     min_len => 100,
                                     min_avg_qual => 30,
                                     min_base_qual => 5,
                                     allow_gaps => 1,
                                     allow_ambig_bases => 1,
                                     });
    
    is( $filter_l->get_fastq_file(), $filename, "get_fastq_file()" );
    is( $filter_l->get_output_dir(), $output_dir, "get_output_dir()" );
    is( $filter_l->get_min_len(), 100, "get_min_len()" );
    is( $filter_l->get_min_avg_qual(), 30, "get_min_avg_qual()" );
    is( $filter_l->get_min_base_qual(), 5, "get_min_base_qual()" );
    is( $filter_l->get_allow_gaps(), 1, "get_allow_gaps()" );
    is( $filter_l->get_allow_ambig_bases(), 1, "get_allow_ambig_bases()" );
}

# test _too_short
{
    my $seq = "ATCG";
    is( BioUtils::QC::FastqFilter::_too_short(10,$seq), 1,
       "_too_short(10,ATCG)" );
    is( BioUtils::QC::FastqFilter::_too_short(undef, $seq), 0,
       "_too_short(undef, ATCG)" );
    is( BioUtils::QC::FastqFilter::_too_short(2, $seq), 0,
       "_too_short(2, ATCG)" );
    is( BioUtils::QC::FastqFilter::_too_short(4, $seq), 0,
       "_too_short(4, ATCG)" );
}

# test _below_avg_qual
{
    is( BioUtils::QC::FastqFilter::_below_avg_qual(30, 30), 0,
       "_below_avg_qual(30,30)" );
    is( BioUtils::QC::FastqFilter::_below_avg_qual(30, 29), 1,
       "_below_avg_qual(30,29)" );
    is( BioUtils::QC::FastqFilter::_below_avg_qual(30, 31), 0,
       "_below_avg_qual(30,31)" );
}

# test _below_base_qual
{
    is( BioUtils::QC::FastqFilter::_below_min_base_qual(20, 20), 0,
       "_below_base_qual(20,20)" );
    is( BioUtils::QC::FastqFilter::_below_min_base_qual(20, 19), 1,
       "_below_base_qual(20,19)" );
    is( BioUtils::QC::FastqFilter::_below_min_base_qual(20, 21), 0,
       "_below_base_qual(20,21)" );
}

# test _has_gaps
{
    is( BioUtils::QC::FastqFilter::_has_gaps("ATGC"), 0,
       "_has_gaps(ATCG)" );
    is( BioUtils::QC::FastqFilter::_has_gaps("ATC-"), 1,
       "_has_gaps(ATC-)" );
    is( BioUtils::QC::FastqFilter::_has_gaps("ATC."), 1,
       "_has_gaps(ATC.)" );
    is( BioUtils::QC::FastqFilter::_has_gaps("-TC."), 1,
       "_has_gaps(-TC.)" );
}

# test _has_ambiguous_bases 
{
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("ATCG"), 0,
       "_has_ambiguous_bases(ATCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("ANCG"), 1,
       "_has_ambiguous_bases(ANCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("ARCG"), 1,
       "_has_ambiguous_bases(ARCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AYCG"), 1,
       "_has_ambiguous_bases(AYCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("ASCG"), 1,
       "_has_ambiguous_bases(ASCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AWCG"), 1,
       "_has_ambiguous_bases(AWCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AKCG"), 1,
       "_has_ambiguous_bases(AKCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AMCG"), 1,
       "_has_ambiguous_bases(AMCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("ABCG"), 1,
       "_has_ambiguous_bases(ABCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("ADCG"), 1,
       "_has_ambiguous_bases(ADCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AHCG"), 1,
       "_has_ambiguous_bases(AHCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AVCG"), 1,
       "_has_ambiguous_bases(AVCG)" );
    is( BioUtils::QC::FastqFilter::_has_ambiguous_bases("AnCG"), 1,
       "_has_ambiguous_bases(AnCG)" );
}

# test _get_file_prefix
{
    is( BioUtils::QC::FastqFilter::_get_file_prefix("my_prefix.txt"),
       "my_prefix", "_get_file_prefix(my_prefix.txt)" );
    is( BioUtils::QC::FastqFilter::_get_file_prefix("my_prefix.tmp.txt"),
       "my_prefix.tmp", "_get_file_prefix(my_prefix.tmp.txt)" );
    is( BioUtils::QC::FastqFilter::_get_file_prefix("something"),
       "something", "_get_file_prefix(something)" );
    is( BioUtils::QC::FastqFilter::_get_file_prefix("/Users/Scott/my_prefix.txt"),
       "my_prefix", "_get_file_prefix(/Users/Scott/my_prefix.txt)" );
}

# test _print_daignostics
{
    my @seq1 = ("seq1",0,0,0,0,0);
    my @seq2 = ("seq2",1,1,1,1,1);
    my @diag_table = (\@seq1, \@seq2);
    
    my $dir = tempdir();
    
    lives_ok( sub{BioUtils::QC::FastqFilter::_print_diagnostics(\@diag_table, $dir, "seqs")},
        "_print_diagnostics - lives");
    
    close($fh);
    
    cmp_ok( -s $filename, ">", 0, "_print_diagnostics - file created" );
}

# test _test_seq
{
    # make a local filter
    my $filter_l = BioUtils::QC::FastqFilter->new({
                        fastq_file => $filename,
                        output_dir => $output_dir,
                        min_len => 4,
                        min_avg_qual => 20,
                        min_base_qual => 3,
                        allow_gaps => 0,
                        allow_ambig_bases => 0,
                        verbose => 1
                    });
    
    # make an array ref as a place holder for the diagnositcs aref required to
    # run _test_seq
    my $aref = ();
    
    # test a good sequence
    is( $filter_l->_test_seq($aref, "seq1", "ATCG", "IIII", 1),
       0,
       "_test_seq -- good seq");
    
    # test a sequence with below avg qual score
    is( $filter_l->_test_seq($aref, "seq1", "ATCG", "++++", 1),
       1,
       "_test_seq -- below avg qual");
    
    # test sequences with below min qual score
    is( $filter_l->_test_seq($aref, "seq1", "ATCG", "III#", 1),
       1,
       "_test_seq -- below min qual");
    
    # test sequences with below min len
    is( $filter_l->_test_seq($aref, "seq1", "ATC", "III", 1),
       1,
       "_test_seq -- below min len");

}

# test filter -- this is the main subroutine
{
    my $prefix = BioUtils::QC::FastqFilter::_get_file_prefix($filename);
    #print "filename: $filename\n";
    #print "output dir: $output_dir\n";
    
    # test with set min_len
    my $filter_l = BioUtils::QC::FastqFilter->new({fastq_file => $filename,
                                    output_dir => $output_dir,
                                    min_len => 8,
                                    min_avg_qual => 20,
                                    min_base_qual => 10,
                                    allow_gaps => 0,
                                    allow_ambig_bases => 0,
                                    verbose => 1,
                                    });
    lives_ok( sub{ $filter_l->filter() }, "filter - lives" );
    
    cmp_ok( length `grep seq1 $output_dir/$prefix\_unfiltered.fastq`, '>', 0,
       "filter(min_len) - seq1 is in unfiltered file" );
    cmp_ok( length `grep seq1 $output_dir/$prefix\_filtered.fastq`, '==', 0,
       "filter(min_len) - seq1 NOT in filtered file" );
    
    cmp_ok( length `grep seq2 $output_dir/$prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq2 is in filtered file" );
    cmp_ok( length `grep seq2 $output_dir/$prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq2 NOT in unfiltered file" );
    
    cmp_ok( length `grep seq3 $output_dir/$prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq3 is in filtered file" );
    cmp_ok( length `grep seq3 $output_dir/$prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq3 NOT in unfiltered file" );
    
    cmp_ok( length `grep seq4 $output_dir/$prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq4 is in filtered file" );
    cmp_ok( length `grep seq4 $output_dir/$prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq4 NOT in unfiltered file" );
    
    cmp_ok( length `grep seq5 $output_dir/$prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq5 is in filtered file" );
    cmp_ok( length `grep seq5 $output_dir/$prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq5 NOT in unfiltered file" );
    
    cmp_ok( length `grep seq6 $output_dir/$prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq6 is in filtered file" );
    cmp_ok( length `grep seq6 $output_dir/$prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq6 NOT in unfiltered file" );
    
    # make sure diagnostics file was create
    # the diagnostic file requires more testing!!
    cmp_ok( -s "$output_dir/$prefix\_diagnostics.txt", '>', 0,
           "filter - diagnostics file created" );
}

# test with pairs (eg given a second fastq file to filter based on paired seqs)
{
    # build a second fastq file
    my ($fh2, $file2name) = tempfile();
    
    # get the output files prefixes
    my $fwd_prefix = BioUtils::QC::FastqFilter::_get_file_prefix($filename);
    
    my $filter_l; # a FastqFilter object for this block
    
    # check when a file is given but is empty
    warnings_exist {my $filter_l = BioUtils::QC::FastqFilter->new({
                                    fastq_file => $filename,
                                    fastq_file2 => file2name,
                                    output_dir => $output_dir})}
        [qr/Empty fastq_file2/,],
        "new empty fastq file2 - warnings";
    
    # populate the new second fastq file
    # NOTE: this method will create a file with a mismatch number of seqs.
    _build_fastq_file3($fh2);
    
    # now it should be created without warnings
    # Here I just call set_fastq_file2 directly
    lives_ok(sub{$filter_l = BioUtils::QC::FastqFilter->new({
                    fastq_file => $filename,
                    fastq_file2 => $file2name,
                    output_dir => $output_dir,
                    min_len => 16,
                    min_avg_qual => 20,
                    min_base_qual => 10,
                    allow_gaps => 0,
                    allow_ambig_bases => 0,
                    verbose => 1,
                    })},
        "lives FastqFilter->new()");
    
    # test the mismatch number of sequences crash with fewer rev seqs
    throws_ok( sub{$filter_l->filter_pairs()},
              qr/Mismatched Pairs: fewer rev seqs/,
              "throws - Mismatched Pairs: fewer rev seqs");
    
    # test the mismatch number of sequences crash with fewer fwd seqs
    throws_ok( sub{
        ($fh2, $file2name) = tempfile();
        _build_fastq_file4($fh2);
        $filter_l->set_fastq_file2($file2name);
        $filter_l->filter_pairs()},
              qr/Mismatched Pairs: fewer fwd seqs/,
              "throws - Mismatched Pairs: fewer fwd seqs");
    
    # test the mismatch of names warnings
    warnings_exist{my ($fh2, $file2name) = tempfile();
        _build_fastq_file5($fh2);
        $filter_l->set_fastq_file2($file2name);
        $filter_l->filter_pairs()}
        [qr/Fwd and Rev IDs do NOT match:/],
        "Warn - Fwd and Rev IDs do NOT match";
    
    # for clearity and testing I create a new output dir
    my $output_dir = tempdir();
    $filter_l->set_output_dir($output_dir);
    
    lives_ok( sub{
                ($fh2, $file2name) = tempfile();
                _build_fastq_file2($fh2);
                $filter_l->set_fastq_file2($file2name);
                $filter_l->filter_pairs()},
             "lives_ok - filter_pairs");
    
    $rev_prefix = BioUtils::QC::FastqFilter::_get_file_prefix($file2name);

    cmp_ok( length `grep seq1 $output_dir/$fwd_prefix\_unfiltered.fastq`, '>', 0,
       "filter(min_len) - seq1 is in fwd unfiltered file" );
    cmp_ok( length `grep seq1 $output_dir/$fwd_prefix\_filtered.fastq`, '==', 0,
       "filter(min_len) - seq1 NOT in fwd filtered file" );
    cmp_ok( length `grep seq1 $output_dir/$rev_prefix\_unfiltered.fastq`, '>', 0,
       "filter(min_len) - seq1 is in rev unfiltered file" );
    cmp_ok( length `grep seq1 $output_dir/$rev_prefix\_filtered.fastq`, '==', 0,
       "filter(min_len) - seq1 NOT in rev filtered file" );
    
    cmp_ok( length `grep seq2 $output_dir/$fwd_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq2 is in fwd filtered file" );
    cmp_ok( length `grep seq2 $output_dir/$fwd_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq2 NOT in fwd unfiltered file" );
    cmp_ok( length `grep seq2 $output_dir/$rev_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq2 is in rev filtered file" );
    cmp_ok( length `grep seq2 $output_dir/$rev_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq2 NOT in rev unfiltered file" );
    
    cmp_ok( length `grep seq3 $output_dir/$fwd_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq3 is in fwd filtered file" );
    cmp_ok( length `grep seq3 $output_dir/$fwd_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq3 NOT in fwd unfiltered file" );
    cmp_ok( length `grep seq3 $output_dir/$rev_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq3 is in rev filtered file" );
    cmp_ok( length `grep seq3 $output_dir/$rev_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq3 NOT in rev unfiltered file" );
    
    cmp_ok( length `grep seq4 $output_dir/$fwd_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq4 is in fwd filtered file" );
    cmp_ok( length `grep seq4 $output_dir/$fwd_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq4 NOT in fwd unfiltered file" );
    cmp_ok( length `grep seq4 $output_dir/$rev_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq4 is in rev filtered file" );
    cmp_ok( length `grep seq4 $output_dir/$rev_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq4 NOT in rev unfiltered file" );
    
    cmp_ok( length `grep seq5 $output_dir/$fwd_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq5 is in fwd filtered file" );
    cmp_ok( length `grep seq5 $output_dir/$fwd_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq5 NOT in fwd unfiltered file" );
    cmp_ok( length `grep seq5 $output_dir/$rev_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq5 is in rev filtered file" );
    cmp_ok( length `grep seq5 $output_dir/$rev_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq5 NOT in rev unfiltered file" );
    
    cmp_ok( length `grep seq6 $output_dir/$fwd_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq6 is in fwd filtered file" );
    cmp_ok( length `grep seq6 $output_dir/$fwd_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq6 NOT in fwd unfiltered file" );
    cmp_ok( length `grep seq6 $output_dir/$rev_prefix\_filtered.fastq`, '>', 0,
       "filter(min_len) - seq6 is in rev filtered file" );
    cmp_ok( length `grep seq6 $output_dir/$rev_prefix\_unfiltered.fastq`, '==', 0,
       "filter(min_len) - seq6 NOT in rev unfiltered file" );
}



### Subroutines ###

sub _build_fastq_file {
    my ($fh_l) = @_;
    
    my $fastq_str = '@seq1' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq1\n" .
                    "IIIIIIII\n" .
                    '@seq2' . "\n" .
                    "AAAATTT\n" .
                    "+seq2\n" .
                    "IIIIIII\n" .
                    '@seq3' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq3\n" .
                    "44444444\n" .
                    '@seq4' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq4\n" .
                    "4444444#\n" .
                    '@seq5' . "\n" . 
                    "AAA-TTTT\n" .
                    "+seq5\n" .
                    "4444444#\n" .
                    '@seq6' . "\n" . 
                    "AAANTTTT\n" .
                    "+seq6\n" .
                    "4444444#\n";
                    
    
    print $fh_l $fastq_str;
    
    close($fh_l);
}

sub _build_fastq_file2 {
    my ($fh_l) = @_;
    
    my $fastq_str = '@seq1' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq1\n" .
                    "IIIIIIII\n" .
                    '@seq2' . "\n" .
                    "AAAATTT\n" .
                    "+seq2\n" .
                    "IIIIIII\n" .
                    '@seq3' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq3\n" .
                    "44444444\n" .
                    '@seq4' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq4\n" .
                    "4444444#\n" .
                    '@seq5' . "\n" . 
                    "AAA-TTTT\n" .
                    "+seq5\n" .
                    "4444444#\n" .
                    '@seq6' . "\n" . 
                    "AAANTTTT\n" .
                    "+seq6\n" .
                    "4444444#\n";
                    
    
    print $fh_l $fastq_str;
    
    close($fh_l);
}

sub _build_fastq_file3 {
    my ($fh_l) = @_;
    
    # This data has a mismatching number of sequences (fewer rev seqs)
    
    my $fastq_str = '@seq1' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq1\n" .
                    "IIIIIIII\n" .
                    '@seq2' . "\n" .
                    "AAAATTT\n" .
                    "+seq2\n" .
                    "IIIIIII\n" .
                    '@seq3' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq3\n" .
                    "44444444\n" .
                    '@seq4' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq4\n" .
                    "4444444#\n" .
                    '@seq5' . "\n" . 
                    "AAA-TTTT\n" .
                    "+seq5\n" .
                    "4444444#\n";
                    
    
    print $fh_l $fastq_str;
    
    close($fh_l);
}

sub _build_fastq_file4 {
    my ($fh_l) = @_;
    
    my $fastq_str = '@seq1' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq1\n" .
                    "IIIIIIII\n" .
                    '@seq2' . "\n" .
                    "AAAATTT\n" .
                    "+seq2\n" .
                    "IIIIIII\n" .
                    '@seq3' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq3\n" .
                    "44444444\n" .
                    '@seq4' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq4\n" .
                    "4444444#\n" .
                    '@seq5' . "\n" . 
                    "AAA-TTTT\n" .
                    "+seq5\n" .
                    "4444444#\n" .
                    '@seq6' . "\n" . 
                    "AAANTTTT\n" .
                    "+seq6\n" .
                    "4444444#\n" . 
                    '@seq7' . "\n" . 
                    "AAANTTTT\n" .
                    "+seq7\n" .
                    "4444444#\n";
                    
    
    print $fh_l $fastq_str;
    
    close($fh_l);
}

sub _build_fastq_file5 {
    my ($fh_l) = @_;
    
    # this data has mimatching sequence names.  It should throw a warning.
    
    my $fastq_str = '@seq1' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq1\n" .
                    "IIIIIIII\n" .
                    '@seq2' . "\n" .
                    "AAAATTT\n" .
                    "+seq2\n" .
                    "IIIIIII\n" .
                    '@seq3' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq3\n" .
                    "44444444\n" .
                    '@seq4' . "\n" . 
                    "AAAATTTT\n" .
                    "+seq4\n" .
                    "4444444#\n" .
                    '@seq5' . "\n" . 
                    "AAA-TTTT\n" .
                    "+seq5\n" .
                    "4444444#\n" .
                    '@seq7' . "\n" . 
                    "AAANTTTT\n" .
                    "+seq7\n" .
                    "4444444#\n";
                    
    
    print $fh_l $fastq_str;
    
    close($fh_l);
}