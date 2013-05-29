
use BioUtils::QC::ContaminantFilter 1.0.2;
use Test::More tests => 43;
use Test::Exception;
use File::Temp qw(tempfile tempdir);

sub _create_seq_file;
sub _create_otu_table_file;
sub _create_blast_db;


# Test use BioUtils::QC::ContaminantFilter
BEGIN { use_ok( 'BioUtils::QC::ContaminantFilter'); }


# Create the temporary files used in testing
my ($seq_fh, $seq_file) = tempfile();
_create_seq_file($seq_fh);

my ($otu_table_fh, $otu_table_file) = tempfile();
_create_otu_table_file($otu_table_fh);

my $blast_db_dir = tempdir();
my $blast_db_name = _create_blast_db($blast_db_dir);

my $output_dir = tempdir();  # temporary output dir

#print "seq_file: $seq_file\n";
#print "otu_table_file: $otu_table_file\n";
#print "blast_dir: $blast_db_dir\n";
#print "output_dir: $output_dir\n";


### Tests ###

# test constructor
my $filter;
{
    # lives with no arguments (args can be set later)
    lives_ok( sub{ $filter = BioUtils::QC::ContaminantFilter->new() },
             'new() - lives'
             );
    
    # lives with args
    lives_ok( sub{ $filter = BioUtils::QC::ContaminantFilter->new({
                                    blast_db => $blast_db_name,
                                    query_file => $seq_file,
                                    output_dir => $output_dir,
                                    output_prefix => "test",
                                    eval => 1,
                                    perc_iden => 94,
                                    output_fmt => 6,
                                    max_targets => 1,
                                    otu_table => $otu_table_file,
                                })},
            "new(arg_href) - lives");
}

# test get/set_blast_db
{
    is( $filter->set_blast_db("$blast_db_name"), 1, "set_blast_db($blast_db_name)" );
    is( $filter->get_blast_db(), $blast_db_name, "get_blast_db()" );
}

# test get/set_query_file
{
    is( $filter->set_query_file("my_query.fasta"), 1,
       "set_query_file(my_query.fasta)" );
    is( $filter->get_query_file(), "my_query.fasta", "get_query_file()" );
}

# test get/set_output_dir
{
    is( $filter->set_output_dir("test_out"), 1,
       "set_output_dir(test_out)" );
    is( $filter->get_output_dir(), "test_out", "get_output_dir()" );
}

# test get/set_output_prefix
{
    is( $filter->set_output_prefix("my_out"), 1,
       "set_output_prefix(my_out)" );
    is( $filter->get_output_prefix(), "my_out", "get_output_prefix()" );
}

# test get/set_eval
{
    is( $filter->set_eval(10), 1, "set_eval(10)" );
    is( $filter->get_eval(), 10, "get_eval()" );
}

# test get/set_perc_iden
{
    is( $filter->set_perc_iden(97), 1, "set_perc_iden(97)" );
    is( $filter->get_perc_iden(), 97, "get_perc_iden()" );
}

# test get/set_output_fmt
{
    is( $filter->set_output_fmt(7), 1, "set_output_fmt(7)" );
    is( $filter->get_output_fmt(), 7, "get_output_fmt()" );
}

# test get/set_max_targets
{
    is( $filter->set_max_targets(1), 1, "set_max_targets(1)" );
    is( $filter->get_max_targets(), 1, "get_max_targets()" );
}

# test get/set_otu_table
{
    is( $filter->set_otu_table("my_table"), 1, "set_otu_table(my_table)" );
    is( $filter->get_otu_table(), "my_table", "get_otu_table()" );
}

# reload the correct filter parameters
$filter = BioUtils::QC::ContaminantFilter->new({  blast_db => $blast_db_name,
                                    query_file => $seq_file,
                                    output_dir => $output_dir,
                                    output_prefix => "test",
                                    eval => 1,
                                    perc_iden => 80,
                                    output_fmt => 6,
                                    max_targets => 1,
                                    otu_table => $otu_table_file,
                                 });

# test _run_blast -- these test are not comprehensive
my $blast_output;
{
    lives_ok( sub{ $blast_output = $filter->_run_blast() },
             '_run_blast() - lives' );
    is( -s $blast_output > 0, 1, 'blast file created');
    #print "blast_output: $blast_output\n";
}


# test _parse_blast_file
my $got; # This is also used in _sequence_printing
{
    my $expected = {'OTU_7'  => 'Arabidopsis_thaliana_chloroplast_rRNA_AP000423.1',
                    'OTU_10' => 'Arabidopsis_thaliana_chloroplast_rRNA_AP000423.1'};

    lives_ok( sub{ $got = $filter->_parse_blast_file($blast_output) },
             '_parse_blast_file - lives');
    is_deeply( $got, $expected, "_parse_blast_file() - deeply" );
}


# test _sequence_printing
{
    lives_ok( sub{ $filter->_sequence_printing($got) },
             "_sequence_printing - lives" );
    is( -s "$output_dir/test_contaminants.fasta" > 0, 1,
       'contaminant fasta file created');
    is( -s "$output_dir/test_non_contaminants.fasta" > 0, 1,
       'non_contaminant fasta file created');
    
    # NOTE: this could be further tested to make sure what ends up in the
    # output files is correct.  Right now I am just testing that there is
    # non-zero sized output files.
}


# test _otu_table_printing
{
    lives_ok( sub{ $filter->_otu_table_printing($got) },
             "_otu_table_printing - lives" );
    is( -s "$output_dir/test_contaminants_otu_table.txt" > 0, 1,
       'contaminants_otu_table file created');
    is( -s "$output_dir/test_non_contaminants_otu_table.txt" > 0, 1,
       'non_contaminants_otu_table file created');
    
    # NOTE: this could be further tested to make sure what ends up in the
    # output files is correct.  Right now I am just testing that there is
    # non-zero sized output files.
}


# test _print_results
{
    # this is mainly tested by _sequence_printing and _otu_table_printing;
    lives_ok( sub{ $filter->_print_results($got) },
             '_print_results() - lives' );
}


# test run_filter
{
    # This is mainly tested by testing all the contained subroutines
    lives_ok( sub{ $filter->run_filter() },
             'run_filter() - lives' );
}

# test for when no otu table is given
{
    my $new_output_dir = tempdir();
    # reload the correct filter parameters
    $filter = BioUtils::QC::ContaminantFilter->new({
                                        blast_db => $blast_db_name,
                                        query_file => $seq_file,
                                        output_dir => $new_output_dir,
                                        output_prefix => "test",
                                        eval => 1,
                                        perc_iden => 80,
                                        output_fmt => 6,
                                        max_targets => 1,
                                     });
    
    lives_ok( sub{ $filter->run_filter() },
             'run_filter() - lives' );
    
    is( -s "$new_output_dir/test_contaminants.fasta" > 0, 1,
       'contaminant fasta file created');
    is( -s "$new_output_dir/test_non_contaminants.fasta" > 0, 1,
       'non_contaminant fasta file created');
    
    is( defined -s "$new_output_dir/test_contaminants_otu_table.txt" > 0, '',
       'contaminants_otu_table file created');
    is( defined -s "$new_output_dir/test_non_contaminants_otu_table.txt" > 0, '',
       'non_contaminants_otu_table file created');
}

# test for when a blast database has not been made but a fasta file is provided
{
    my $new_output_dir = tempdir();
    my $blast_dir = tempdir();
    my $new_blast_db = _create_contam_file($blast_dir);
    # reload the correct filter parameters
    $filter = BioUtils::QC::ContaminantFilter->new({
                                        blast_db => $new_blast_db,
                                        query_file => $seq_file,
                                        output_dir => $new_output_dir,
                                        output_prefix => "test",
                                        eval => 1,
                                        perc_iden => 80,
                                        output_fmt => 6,
                                        max_targets => 1,
                                        otu_table => $otu_table_file,
                                     });
    
    lives_ok( sub{ $filter->run_filter() },
             'run_filter() - lives' );
    
    is( -s "$new_output_dir/test_contaminants.fasta" > 0, 1,
       'contaminant fasta file created');
    is( -s "$new_output_dir/test_non_contaminants.fasta" > 0, 1,
       'non_contaminant fasta file created');
    
    is( -s "$new_output_dir/test_contaminants_otu_table.txt" > 0, 1,
       'contaminants_otu_table file created');
    is( -s "$new_output_dir/test_non_contaminants_otu_table.txt" > 0, 1,
       'non_contaminants_otu_table file created');
}




### Subroutines ###
sub _create_seq_file {
    my ($fh) = @_;
    
    my $str = <<SEQ;
>OTU_1
GTGCCAGCCGCCGCGGTAATACGAGTGCCTCAAGCGTTATCCGGAATCATTGGGCGTAAAGGTTGTGTAGGTGGTTTTAT
TAGTCTTCTGTTAAATTCTTCGGCTTAACCGGGGGCATGCAGAGGAAACGGTAAAACTAGAGGATGCGAGGGGTTAGCGG
AACTCATAGTGTAGCGGTGAAATGCGTTGATATTATGGGGAACACCAAAAGCGAAGGCAGCTAACTGGAGCACTCCTGAC
ACTGAAACAAGAAAGCGTGGGTCGCGAATGGGATTAGATACCCCTGTAGTCC
>OTU_2
GTGCCAGCAGCCGCGGTAAGACGAACCGTGCAAACGTTATTCGGAATCACTGGGCATAAAGGGCGCGTAGGCGGCTCCAA
AAGTCAGGGGTGAAATCCGGCAGCTTAACTGTCGCAGTGCCTTTGATACTGTGGAGCTAGAGGGAGGTAGGGGTCTGTGG
AACTTCCGGTGGAGCGGTGAAATGCGTTGATATCGGAAGGAACGCCGGTGGCGAAAGCGACGGACTGGATCTCTTCTGAC
GCTGAGGCGCGAAAGCCAGGGGAGCAAACGGGATTAGATACCCCGGTAGTCC
>OTU_3
GTGCCAGCCGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACGCTT
AAGTCAGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTTTGATACTGGGTGTCTGGAGGTCGAGAGAGGTGAGTGG
AATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGAC
GCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCGGGTAGTCC
>OTU_4
GTGCCAGCCGCCGCGGTAATACGGAGGGTGCAAACGTTGCTCGGAATCATTGGGCGTAAAGCGCACGTAGGCGGCCCGCT
AAGTCGGATGTGAAATCCCTCGGCTCAACCGAGGACGTGCATTCGATACTGGCAGGCTTGAATATGGAAGAGGGTCGCGG
AATTCCCGGTGTAGAGGTGAAATTCGTAGATATCGGGAGGAACACCAGTGGCGAAGGCGGCGACCTGGGCCAATATTGAC
GCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCGGTAGTCC
>OTU_5
GTGCCAGCCGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTCGCT
GTGTCCGATGTGAAAGCCCCGGGCTTAACCTGGGAATGGCATTGGAAACTGGCGGGCTTGAGTGCGGCAGAGGGGGGTGG
AATTCCGCGTGTAGCAGTGAAATGCGTAGAGATGCGGAGGAACACCGATGGCGAAGGCAGCCCCCTGGGTCGACACTGAC
GCTCAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCCGGTAGTCC
>OTU_6
GTGCCAGCCGCCGCGGTAATACAGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGGGTGCGTAGGCGGTTGTTT
AAGTTAGATGTGAAATCCCCGGGCTTAACCTGGGAATTGCATTTAATACTGAATGACTAGAGTAGAGTAGAGGGAAGTGG
AATTTCCGGTGTAGCGGTGAAATGCGTAGATATCGGAAAGAACATCAGTGGCGAAGGCGGCTTCCTGGACTCATACTGAC
GCTGAGGCACGAAAGCATGGGGAGCAAACAGGATTAGATACCCCAGTAGTCC
>OTU_7
GTGCCAGCCGCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCCGCAGGTGGCGATGT
AAGTCTGCTGTTAAAGAGCAAAGCTTAACTTTGTAAAAGCAGTGGAAACTACATAGCTAGAGTACGTTCGGGGCAGAGGG
AATTCCTGGTGTAGCGGTGAAATGCGTAGAGATCAGGAAGAACACCGGTGGCGAAGGCGCTCTGCTAGGCCGTAACTGAC
ACTGAGGGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCTGTAGTCC
>OTU_8
GTGCCAGCCGCCGCGGTAATACAGAGGGTGCAAGCGTTGCTCGGAATCATTGGGCGTAAAGGGCAAGTAGGTGGTCTCAT
TTGTCTTGTGTGAAATCCTTGGGCTTAACTCAAGAAGTGCGCAAGAAACGGTGAGACTCGAGTTCTGGAGAGGGTCGTGG
AATTCCCGGTGTAGCGGTGAAATGCGTAGAGATCGGGAGGAACACCAGAGGCGAAGGCGGCGACCTGGACAGATACTGAC
ACTCAACTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCGGGTAGTCC
>OTU_9
GTGCCAGCCGCCGCGGTAATACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGTGCGTAGGTGGGTTTGT
AAGTCAGTGGTGAAATCTCCAGGCTTAACCTGGAAACTGCCATTGATACTATAGATCTTGAATTTTCTGGAGGTAAGCGG
AATATGTCATGTAGCGGTGAAATGCTTAGATATGACATAGAACACCAATTGCGAAGGCAGCTTGCTACAGGGATATTGAC
ACTGAGGCACGAAAGCGTGGGGATCAAACAGGATTAGATACCCCTGTAGTCC
>OTU_10
GTGCCAGCAGCCGCGGTAAGACGAACCGCCCGAACGTTGTTCGGATTCACTGGGCTTAAAGGGCGCGTAGGCGGGCGGCC
CGGTCGGGGGTGAAATCTTTCAGCTCAACTGGAAAAGAGCTTCCGATACCGGCCGTCTGGAGGGAGGTAGGGGCATGCGG
AACTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAGGAACGCCGGTGGCGAAAGCGGCGTGCTGGACCTCTTCTGAC
GCTGATGCGCGAAAGCCAGGGGAGCGAACGGGATTAGATACCCCAGTAGTCC
SEQ

    print $fh $str;
    
    close($fh);
    
    return 1;
}

sub _create_otu_table_file {
    my ($fh) = @_;
    
    my $str = <<OTU;
# QIIME v1.3.0 OTU table
#OTU ID P1      P10     P11     P12     P13     P14     P15     P16     P17     P18     P19     P2      P20     P21     P22     P23     P24     P25     P26     P27   
1       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
2       0       0       1       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
3       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
4       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
5       0       0       0       0       0       0       0       1       0       0       0       0       0       0       0       0       0       0       0       1
6       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
7       0       0       0       2       0       1       0       0       0       0       1       0       0       0       0       0       0       0       0       0
8       0       0       0       0       0       1       0       0       0       0       0       0       0       0       0       0       0       0       0       0
9       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
10      0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0
OTU

    print $fh $str;
    
    close($fh);
    
    return 1;
}

sub _create_blast_db {
    my ($dir) = @_;
    
    my $file = _create_contam_file($dir);
    
    
    # make the blast db -- I used qx{} instead of system to supress output
    `makeblastdb -in $file -dbtype nucl 2>/dev/null`;
    
    return "$file";
}

sub _create_contam_file {
    my ($dir) = @_;
    
    my $db_seqs = ">Arabidopsis_thaliana_chloroplast_rRNA_AP000423.1
AAGGAAGCTATAAGTAATGCAACTATGAATCTCATGGAGAGTTCGATCCTGGCTCAGGATGAACGCTGGCGGCATGCTTAACACATGCAAGTCGGACGGGAAGTGGTGTTTCCAGTGGCGGACGGGTGAGTAACGCGTAAGAACCTGCCCTTGGGAGGGGAACAACAGCTGGAAACGGCTGCTAATACCCCGTAGGCTGAGGAGCAAAAGGAGGAATCCGCCCGAGGAGGGGCTCGCGTCTGATTAGCTAGTTGGTGAGGCAATAGCTTACCAAGGCGATGATCAGTAGCTGGTCCGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCCGCAATGGGCGAAAGCCTGACGGAGCAATGCCGCGTGGAGGTAGAAGGCCTACGGGTCCTGAACTTCTTTTCCCAGAGAAGAAGCAATGACGGTATCTGGGGAATAAGCATCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCTAGCCGTAAACGATGGATACTAGGCGCTGTGCGTATCGACCCGTGCAGTGCTGTAGCTAACGCGTTAAGTATCCCGCCTGGGGAGTACGTTCGCAAGAATGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAAAGCGAAGAACCTTACCAGGGCTTGACATGCCGCGAATCCTCTTGAAAGAGAGGGGTGCCTTCGGGAACGCGGACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGCCGTAAGGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTGTTTAGTTGCCACCGTTGAGTTTGGAACCCTGAACAGACTGCCGGTGATAAGCCGGAGGAAGGTGAGGATGACGTCAAGTCATCATGCCCCTTATGCCCTGGGCGACACACGTGCTACAATGGCCGGGACAAAGGGTCGCGATCCCGCGAGGGTGAGCTAACTCCAAAAACCCGTCCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGCCGGAATCGCTAGTAATCGCCGGTCAGCCATACGGCGGTGAATTCGTTCCCGGGCCTTGTACACACCGCCCGTCACACTATGGGAGCTGGCCATGCCCGAAGTCGTTACCTTAACCGCAAGGAGGGGGGTGCCGAAGGCAGGGCTAGTGACTGGAGTGAAGTCGTAACAAGGTAGCCGTACTGGAAGGTGCGGCTGGATCACCTCCTTTTCAGGGAGAGCTAATGCTTCTT";

    open my $fh, ">", "$dir/contaminants.fasta" or
        croak("Cannot open file: $dir/contaminants.fasta");
    
    print $fh $db_seqs;
    
    close($fh);
    
    return "$dir/contaminants.fasta";
}

