use strict;
use warnings;

use Test::More tests => 57;
use Test::Exception;
use Test::Warn;
use BioUtils::FileReformat::GffToFas;
use Log::Log4perl qw(:easy);
use Log::Log4perl::CommandLine qw(:all);
my $logger = get_logger();

# others to include
use File::Temp qw/ tempfile tempdir /;

# helper subroutines
sub _make_test_gff_file;
sub _make_test_genome_file;


# make a temp input gff file for testing
my ($gff_fh, $gff_filename) = tempfile();
_make_test_gff_file($gff_fh);
#print("gff file: $gff_filename\n");

my ($genome_fh, $genome_filename) = tempfile();
_make_test_genome_file($genome_fh);
#print("genome file: $genome_filename\n");

# make a temp output fas file
my ($fas_fh, $fas_filename) = tempfile();
#print("fas temp file: $fas_filename\n");

# make an empty file -- for testing
my ($empty_fh, $empty_filename) = tempfile();


### Begin Testing ###
BEGIN { use_ok( 'BioUtils::FileReformat::GffToFas' ); }

# test constructor
my $obj = undef;
{
    throws_ok( sub { BioUtils::FileReformat::GffToFas->new({
                        gff_file => $gff_filename}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    throws_ok( sub { BioUtils::FileReformat::GffToFas->new({
                        gff_file => $gff_filename}) },
              'MyX::Generic::Undef::Param', "new() - caught" );
    lives_ok( sub { $obj = BioUtils::FileReformat::GffToFas->new({
                                gff_file => $gff_filename,
                                genome_file => $genome_filename,
                                fas_file => $fas_filename,
                                save_feature => "CDS"}) },
             "expected to live" );
    
    # test constructor when save_feature is not provided
        lives_ok( sub { $obj = BioUtils::FileReformat::GffToFas->new({
                                gff_file => $gff_filename,
                                genome_file => $genome_filename,
                                fas_file => $fas_filename}) },
             "expected to live" );
}

# Test the simple getter methods
{
    lives_ok( sub { $obj->get_gff_file() }, "expected to live" );
    is ($obj->get_gff_file(), $gff_filename, "get_gff_file()" );
    
    lives_ok( sub { $obj->get_genome_file() }, "expected to live" );
    is ($obj->get_genome_file(), $genome_filename, "get_genome_file()" );
    
    lives_ok( sub { $obj->get_fas_file() }, "expected to live" );
    is ($obj->get_fas_file(), $fas_filename, "get_fas_file()" );
    
    lives_ok( sub { $obj->get_save_feature() }, "expected to live" );
    is ($obj->get_save_feature(), "CDS", "get_save_feature()" );
    
    lives_ok( sub { $obj->get_print_exons() }, "expected to live" );
    is( $obj->get_print_exons(), "0", "get_print_exons()" );
    
    lives_ok( sub { $obj->get_fasta_type() }, "expected to live" );
    is( $obj->get_fasta_type(), "nucl", "get_fasta_type()" );
}

# Test the simple setter methods
{
    throws_ok( sub { $obj->set_gff_file() },
              'MyX::Generic::Undef::Param', "set_gff_file() - caught" );
    throws_ok( sub { $obj->set_gff_file("blah.gff") },
              'MyX::Generic::DoesNotExist::File', "set_gff_file(blah.txt) - caught" );
    throws_ok( sub { $obj->set_gff_file($empty_filename) },
              'MyX::Generic::File::Empty', "set_gff_file(empty) - caught" );
    lives_ok( sub { $obj->set_gff_file($gff_filename) },
            "expected to live" );
    is( $obj->get_gff_file(), $gff_filename, "set_gff_file()" );
    
    throws_ok( sub { $obj->set_genome_file() },
              'MyX::Generic::Undef::Param', "set_genome_file() - caught" );
    lives_ok( sub { $obj->set_genome_file($genome_filename) },
            "expected to live" );
    is( $obj->get_genome_file(), $genome_filename, "set_genome_file()" );
    
    throws_ok( sub { $obj->set_fas_file() },
              'MyX::Generic::Undef::Param', "set_fas_file() - caught" );
    lives_ok( sub { $obj->set_fas_file($fas_filename) },
            "expected to live" );
    is( $obj->get_fas_file(), $fas_filename, "set_fas_file()" );
    
    lives_ok( sub { $obj->set_save_feature("test") },
            "expected to live" );
    lives_ok( sub { $obj->set_save_feature() },
             "exected to live" );
    is( $obj->get_save_feature(), "CDS", "set_save_feature()" );
    
    lives_ok( sub {$obj->set_print_exons(1) },
             "expected to live" );
    lives_ok( sub { $obj->set_print_exons() },
             "exected to live" );
    is( $obj->get_print_exons(), 0, "set_print_exons()" );
    
    lives_ok( sub {$obj->set_fasta_type("prot") },
             "expected to live" );
    lives_ok( sub { $obj->set_fasta_type() },
             "exected to live" );
    is( $obj->get_fasta_type(), "nucl", "set_fasta_type()" );
    throws_ok( sub { $obj->set_fasta_type("blah") },
              'MyX::Generic::BadValue', "set_fasta_type(blah) - caught" );
}

# test the _compare_plus method
{
    
    # set up two object to compare
    my $seq1 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=1  start=7917 end=8087 strand=+ contig=scaffold_1",
        seq => "ATCG"
    });
    my $seq2 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=2  start=8123 end=8355 strand=+ contig=scaffold_1",
        seq => "GCTA"
    });
    
    cmp_ok(BioUtils::FileReformat::GffToFas->_compare_plus($seq1, $seq2), '<', 0,
           "_compare_plus");
}

# test the _compare_minus method
{
    
    # set up two object to compare
    my $seq1 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=1  start=7917 end=8087 strand=- contig=scaffold_1",
        seq => "ATCG"
    });
    my $seq2 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=2  start=8123 end=8355 strand=- contig=scaffold_1",
        seq => "GCTA"
    });
    
    cmp_ok(BioUtils::FileReformat::GffToFas->_compare_minus($seq1, $seq2), '>', 0,
           "_compare_minus");
}

# test the _get_gene_header method
{
    # set up two object to compare
    my $seq1 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=1  start=7917 end=8087 strand=+ contig=scaffold_1",
        seq => "ATCG"
    });
    my $seq2 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=2  start=8123 end=8355 strand=+ contig=scaffold_1",
        seq => "GCTA"
    });
    
    my @arr = ($seq1, $seq2);
    my $new_header = BioUtils::FileReformat::GffToFas->_get_gene_header(\@arr);
    is($new_header, "695827 NA name='gm1.2_g' exonNumber=2  start=7917 end=8355 strand=+ contig=scaffold_1",
       "_get_gene_header");
}

# test the _get_gene_seq method when fasta_type == nucl
{
    
    # set up two object to compare
    my $seq1 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=1  start=7917 end=8087 strand=+ contig=scaffold_1",
        seq => "ATCG"
    });
    my $seq2 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=2  start=8123 end=8355 strand=+ contig=scaffold_1",
        seq => "GCTAA"
    });
    
    my @arr = ($seq1, $seq2);
    my $new_seq_obj = $obj->_get_gene_seq(\@arr);
    is($new_seq_obj->get_seq(), "ATCGGCTAA", "_get_gene_seq");
    
    # test when there is only one object in the arr
    @arr = ($seq1);
    $obj->set_fasta_type("nucl");
    $new_seq_obj = $obj->_get_gene_seq(\@arr);
    is($new_seq_obj->get_seq(), "ATCG", "_get_gene_seq");
}

# test the _get_gene_seq method when fasta_type == prot
{
    $obj->set_fasta_type("prot");
    $obj->set_print_exons(0);
    
    # set up two object to compare
    my $seq1 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=1  start=7917 end=8087 strand=+ contig=scaffold_1",
        seq => "ATGG"
    });
    my $seq2 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=2  start=8123 end=8355 strand=+ contig=scaffold_1",
        seq => "GCTAA"
    });
    my $bad_seq2 = BioUtils::FastaSeq->new({
        header => "695827 NA name='gm1.2_g' exonNumber=2  start=8123 end=8355 strand=+ contig=scaffold_1",
        seq => "GCTA-"
    });
    my $prot_seq;
    
    # test when there is a bad codon
    my @arr = ($seq1, $bad_seq2);
    throws_ok( sub { $obj->_get_gene_seq(\@arr) },
              'BioUtils::MyX::Fasta::BadCodon', "_get_gene_seq() - caught" );
    
    # test when there is not %3 bases
    @arr = ($seq1);
    warning_is { $prot_seq = $obj->_get_gene_seq(\@arr) } 'Seq (695827) is not divisible by 3.  Trailing bases are ignored.', "Trailing bases warning";
    is( $prot_seq->get_seq(), "M", "_get_gene_seq() - prot = M");
    
    # test when it should work
    @arr = ($seq1, $seq2);
    lives_ok( sub { $prot_seq = $obj->_get_gene_seq(\@arr) }, "_get_gene_seq() - expected to live");
    is( $prot_seq->get_seq(), "MG", "_get_gene_seq() - MG" );
}

# test the _init method
{
    lives_ok( sub { $obj->_init({
                                 gff_file => $gff_filename,
                                 genome_file => $genome_filename,
                                 fas_file => $fas_filename}) },
             "_init() - expected to live" );
    is( $obj->get_gff_file(), $gff_filename, "get_gff_file() after _init()" );
    is( $obj->get_genome_file(), $genome_filename, "get_genome_file() after _init()" );
    is( $obj->get_fas_file(), $fas_filename, "get_fas_file() after _init()" );
    is( $obj->get_save_feature(), "CDS", "get_save_feature() after _init()" );
    is( $obj->get_print_exons(), "0", "get_print_exons() after _init()" );
    is( $obj->get_fasta_type(), "nucl", "get_fasta_type() after _init()" );
}

# test the reformat method -- This is the most important one
# unfortunately, all I can do is test that the output file was created and is
# non-empty
{
    lives_ok( sub { $obj->reformat() },
             "reformat() - expected to live" );
    cmp_ok( -s $obj->get_fas_file(), '>', 0, "fas is non-empty ($fas_filename)" );
}












sub _make_test_gff_file {
    my ($fh) = @_;
    
    my $str = "scaffold_1	JGI	exon	7917	8087	.	+	.	name='gm1.2_g';transcriptId=695887
scaffold_1	JGI	CDS	7917	8087	.	+	0	name='gm1.2_g';proteinId=695827;exonNumber=1
scaffold_1	JGI	start_codon	7917	7919	.	+	0	name='gm1.2_g'
scaffold_1	JGI	exon	8123	8355	.	+	.	name='gm1.2_g';transcriptId=695887
scaffold_1	JGI	CDS	8123	8355	.	+	0	name='gm1.2_g';proteinId=695827;exonNumber=2
scaffold_1	JGI	exon	8494	8737	.	+	.	name='gm1.2_g';transcriptId=695887
scaffold_1	JGI	CDS	8494	8737	.	+	1	name='gm1.2_g';proteinId=695827;exonNumber=3
scaffold_1	JGI	stop_codon	8735	8737	.	+	0	name='gm1.2_g'
scaffold_1	JGI	exon	10204	10506	.	-	.	name='estExt_fgenesh1_pg.C_1_t10001';transcriptId=683617
scaffold_1	JGI	CDS	10204	10506	.	-	0	name='estExt_fgenesh1_pg.C_1_t10001';proteinId=683557;exonNumber=3
scaffold_1	JGI	stop_codon	10204	10206	.	-	0	name='estExt_fgenesh1_pg.C_1_t10001'
scaffold_1	JGI	exon	10581	11442	.	-	.	name='estExt_fgenesh1_pg.C_1_t10001';transcriptId=683617
scaffold_1	JGI	CDS	10581	11442	.	-	1	name='estExt_fgenesh1_pg.C_1_t10001';proteinId=683557;exonNumber=2
scaffold_1	JGI	exon	11634	11736	.	-	.	name='estExt_fgenesh1_pg.C_1_t10001';transcriptId=683617
scaffold_1	JGI	CDS	11634	11647	.	-	0	name='estExt_fgenesh1_pg.C_1_t10001';proteinId=683557;exonNumber=1
scaffold_1	JGI	start_codon	11645	11647	.	-	0	name='estExt_fgenesh1_pg.C_1_t10001'
";

    print $fh $str;
    
    close($fh);
    
    return 1;
}

sub _make_test_genome_file {
    my ($fh) = @_;
    
    my $str = <<'FILE';
>scaffold_1
atataataaaaaagtctttatttatataataataaaaaatattataaagattaaaatatagtatttaata
attcaattatattaaaaatattaaattaaggaaaagagaaaagaagaaaagagaaatcaaaaagagaaaa
gtattttataatcttttaattaataataaatataagttttaaaaatattaatataataaatataataatc
ggaatatttattatataattaatttaattataattatattaGTNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCATTTTTCGCACCCatttttaa
aaaaattaaattatCACGAATAATAGTCTAAATCACGTGTAATTTTTTTCTCTAGTAGTCCATTAACTAG
GCTTAAATTACTAGTCAAACGAAATCCATATAAAGCACTTAATTAATACGGTATTTTAAGTATTAAATAA
GTGATACAAATCGTGAAAACTCTTAAATTTAAGGCCGCCGCACCTATATAACTCCCGTGTATCCATTTCG
TATCCATTTCGCACCCATTTTCGACACCCattttttaaaaaaataaattatCACGAATAACAGTCCGAAT
CACGTGCAAATTTTTTTCCTAACAATCCATTAACTAGGCCCAAATTACTAGTCAAATGAAATCCATATAA
AGCACTTAATTAATACGGCATTTTAAGTATTAAATAAATAGTACAAATCGTGAAAACCCTTAAATTTAAG
GCCGCCGCACCTATACCACTCCTGTACACCCATTTCGCACCCATTTTCGGCACCCatttttttaaaaatt
aaattatCACGAATAATTATCTAAATCACATGTAATTTTTTTCCCTAGTAGCCCATTAATTAGCTCCAAA
TTACTAGTTAAAAGGAATCCTTATAATGCCGTTTATTAAAGCGGTATTTTAAGTCTTAAATAAGTGATGA
AAATCGTGATTTTGGTATTCCTACCTTATACAATAACGTCATTATAACTATACCTCAAAATCTCAATATC
GGGTTCAAATTTTTAAAAAATTACAACTTCCTTAATAATAAGGCTAAGCACGTGATTTTATTTTTTATAA
ATAGGCCTATTTCAAGTCTATAATTTAAATTAGACTAATATCACGTGACTTATCAAGTGCTAGCGCATTA
AATAAGAATGTGATTTTAAATTCGGGAAATTACTCCCTTGGATTTAAGGGCTCTTATTGGCTTAAGGGCC
CACACTAGGCCCAAACAGGAGATTCACATCACTAATCTTGCAAAAACCTTATCTTTTTTTGCCAAGACTT
TAGTGTCCCTTATTAGTAAAAACCCCTATTCTAGAGGAGTGTACTCAGCGAAAATCTAACCtttaaatct
attatttaaagatattatataaatccttttctatttttataatttacttttctattatattattaaatat
aataatagttttaaaaGAATATTCGATTTTCGATTATTTTCGATTTTATTATTTTAAGTAATAAATAATt
ttatatttatatttataatcagattctctattaaatttatattcttatcttctcttatatattatactct
ttttatttatttttatattttatatttTAGAATCAAATTTATACTCTTTTATAAAAACTTTAAAAGTTAC
TTGAACTTGATTTAGAGTATTAAAGAATTGaattaatattaatatatttatagaataattaaagtttatt
tttaataataaaaaatatataataataataatatacttttaaataatagatttaaaGATTAGATTTTCGC
TTGGTTTACCTATACTAGAATAGGGATTTTTACTAATAAGGGACATTAAAGTCTTAGTAAAAAAAGATAA
GGTTTCTGTAAGATTAGTGATATAAATCTCTTGTTTGGACCTAGTATAGGCCCTTGGGCCAATAAGAGCC
CTTAAATTCAAAGGGATAATCTCCCGGAATTTAAAATCACATTCTTATTTAATATACTAGTATTTAATAA
ATCACGTAATATTAGTCTAATTTAAATTATAGACCTGAACTAGGCCTATTTCAGGAAATTAAAATCACGT
GCTTAGTCttattattaaagaagttataattttttaaaaaaaattaaaaccgatattaaagttttgaaat
ataattataataatattattatatataataGAAATACTAAAATCACGATTTCTATAATTTATTTAAGACT
TGAAATGCCGCCTTAATTAACCGTATTATAAGGATTCTATTTAACTAGTAATTTGGAACTAATTGATAGG
TTACTAGaaaaaaaaagtatatataatttagactattatttataataatttaattttttaaaaaatGGGT
GCGGAAAATGAATACGAAATGGATGCGAAATGAGTACACGAGAGTTATATAGGTGCGGCGGCCTTAAATT
TAAGAGTCTTTACAATTTCTATCACttatttaatatttaaaatattactttaattaaGTACTTTACGGAG
ATTTCGTTTAACTAGTAATTTAGGCCTAGTTAATGGGTTATTGAAAAAAAAAGTTTACACGTGATTCGGA
CTATTATTCGTAATGATTTATTTTTTTTAAAAATGAGTGCGGAAAATGGATACGAAATAtataataatcg
ggatacccgttatataattagcttgattgtgattatattagtttcagattaagggaatcgatatttgtaa
tattaatatttcagattaagggtatgagtaagtattgggagagtataagacacttctctcttcttgactt
ctcttctcttctctttaatttagtattcctaatataactggattattgaatattatatcttgacttttat
agtatttcctattattataAAATGGGTATGAAATAATATACGGGaattatataagtataataactttaaa
tttaaaagtctttataatttttattatttatttaatatttaaaataCTGTTTTAATTAAGTACTTTATAA
GGATTTTATTTGATTAATAATTTAGGCCTAGTTAATGGATTATTAAAAAAAAGTTTATATGTGATCCGGA
TTATTATTCGTAATAGtttaatttttttaaaaattaaaacggaaaataaatataataaaGANNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAATTTAAAATTA
ATCAATTAATATCACGTGATTTGGtaaattaaaatactttaatttaatataaaattgaatataattaaat
aaattctttatataaatattaaaatttaaaaaattaaatttaaAATTTCCTTTTAATTAATAAACTAATA
AGCAATCACGTGAGTATTAATATAAGCTTTTAACTATTAAATCACAAAGCTAATCATGtattttaatttt
aattttaatatttaatataaataataattatttattttaaaaatattaattttttttaaaaaaaaagaaa
ttaaaaattaaaattattaaaacttattattaaaGCTCACGTGACTTTGCGCATATAAATATATTAAATT
AACGTGTCTTTAACTTGAATTAAACGAACttttaatatttacttttatttttaaaaatatactcgaattt
tttaatttttattaaatttaaacttttataaaaaagtaaaccttgatttaaagactcttgaattgatcga
aggctcgattatactattaaaaaagaatataataaaaaattaaaagttaatttattaatatattaatata
tatattataaaataatatataaattcatttattatttttagtaaagaaatattactaaaattaaaCGAAA
AATCTAATAAGTACTTAAAAGTATATATACTCGTATTAATAATCTCAAATCTAGTAATCTTTTTAAATAT
TATCTTTATTACTCATTTACTTGAATATTTACTTATTATAAAATAACCTTAATATAAGTATTTATAATCA
GATTATCTAGTAAATAAATATTAGAAATAATTCACTTTTATATATTATAGTATAGTATATAATAGCTTTA
TAAAGACAGAAGTACTTGATATTACTAGGAAACCACTATAAAAAGTAAGTACCTTAGCCTTTCACGAATC
ACGATTTTTCGGATCTTACTATACTGATTTATCTCCGAGCAGAATCGCCCAAAACGGCCATTTTCAACTG
GCCACCTCGATCGGGAAGCGTGCTCGGCCCATAACCCGATTAGGTTAAGCGGGAGGTGAGTCGACTTAAA
GTGAGTCGACTTAGGGAAAGCCGCCCCCGAGGTCTCTCGTTGAACCGCGTTTGGATGCTCCATTTGGGGC
TCTTTCTGTTGGTCACATATGCTGCAACTTGGCGTTCGATGAAGATGATCGAGTAAGCCGAGGGCTGCCC
GTGAAGATGAAGTATGCAACCCCAGATGCGACCCCTCCACGAGAGCCGGGTATTCCCCGTTTGATTCGCT
GGTTGCACAGCGTGGTTCGAGATCTCCAGCACATCAGTATCTTGTTAGCTCCGAGGCTTGTCCCGTCGAT
TCAGGATCCACAATATGCCAACAGAAGGCGCCAGCCTGTATGAGTAAGCCAAGCACCTGCACATCCAATC
AGAAATACTCCGCAATTCTGGATCTGCTCAAGTGTCGTGCAGACCTTGGAAGCGATGCACCTCTAGACAG
AGTTGAAAGCACGGCGTTTGCGGTTAATTCTCCGTACAACTTCTTACCGCACTGTGCACGCAAAGGACTA
CAGCCACTGCCAGCATTTTCCACATCAATCGAGCAAGAATGGTGTCGTCCAGCTCAAAGAGCTGGAAGAA
AAGACCTCTATACCACATTGATGGTGTAAGAGTAACTACCGATGGCAGCAGCTACGCAATTTTGTTTGAT
ATCATAAATATCCAGCCTTCAATACCGTCTTCTTCGTCTTTTCTCTGCAACCACACCTCTGCAGTTTTTG
CTCCGACTCTAAGAACTTCCTAGCCGGCAAATCTACTACTTTGCCCGAGCCGCGCTTACGGCAACCATCT
CAAGCCTCACCTCGAAGCTGCAAACCATCGGTAGAAAACTAACACAATCTGCAATAAGAACTTTAATTCC
CGTTCAACTCTCACCACATGCATCTGCTATCCTCCACAAGATCATCGCCTGCGGATCGAGCCGCTGCGGC
ATTATGCGAACCTGTCTGATTGCAGATGCTAGTTCTGCAGCTCTGAGTGATGGCTAGATAATCATGCTGG
CATCAAGCGATGAAGTGTCAGCGCAGTCGGGCACACCATTTATCCCTAAGTCAGAACTCTATGTTACTTC
TCCAGTGAGGCGAAAATGAGGACAATCTCTCGCACACCCACCTAACGAATTGGTGGTGCATCGATAGTTC
TCTGCAACGAGGAGTCTCTTCTTACTGCAGATACTAGCTGCAGCTTGGAGTGAAGGATATGTGTTCATGC
TCACTTCACATGAGGAAGTATTAGCGCTATTGGCACGCTATTCCATCTCTCAAGTAAGAACTCTATCCTC
ATTCTCCAGTGAGGAGACAGTGAGGACCACCACTCGCACACCATCCAAATGACTTGTTGGTAACATCAAT
AGTTTCTGCAGCAAGGAACCTTTTCTCACTGCGGAACATAGTTCTGCAGCTTTGAGTAAAGAACCGATCC
TCGCATCAAGTGACGAGGTATTAACGCTGTTGGCACGCCATTCAATCTCTGAAGTCAGAACTCTATTCTC
GTTCTCCAGTAAGAGAACAACGAGGATCACCACCTTGATGAAATGGTAGTGAATTCAATACTTTCTCTGC
AGCTAGAAACCTTTTATCACTGCAGAAACAAGTTCTGCAGGTTTGAGTGAAGGAGCAATGCTACTTAATT
CATGTGATGAGATACCAACGCTGTTGGTACACCATTCAATGGCCCAAGTCAGATCTCTATTTTCGTTCTC
CAGCGAGAAGACCATGAGGGCCACTACTCTCACACCACCCAAATGTTTTGTTGGTGATTCAATACTTCTT
TTGCATCGAGGGACCACAGCCTCGCACCACCCTCAGCGAGCAAAGTAAAGACATGTTGGTGAAATCAATA
TTGATTCTGCAATGCATTCTGATTCGGGAAGTTCCCAGGAAAGCCACGACAGTTACTCTGAAATCCTTTT
GTAGATAGGCAGACACAGTACGTCGAGCAAGCTTTGAGATTGGATACTACGCCCTATTTCTGGCAGAAAG
ACATGGAAGTACACACACTTCTGCAACCTTCATCACACTCAACTTTGCAGTTCGATAAATCCCTCAACAT
ATTCTTGACTCCGAGCGCCTGTACCCGACTTCAGGATTACTTTGGATCTTGCTGCCGGCGGCAATGCAAG
TTGATGCTATATTTGTGAAGTACAGAGAGGTACTTAGATGTATCTTTTCTTTGATCATGTTCTTGTGGGT
CAGGAGGAGAGAACGATAGTCATACAAATTGGTGTGGCTGGGTGGTGTAGTTGGTGTTGGTCTTGAGGGT
TACTGGTAGGAAAGATGGGGGGGAAGGGGTTGAGAGCGAGGGGTGCCTGAGACGCGTTTCGGGTTCGGAA
AGTGCGTGGTAGAGATATCACATTCTCTTTACAAACCTTTTTGGAATACTAGTGCTTCAAGAACATGATT
CATTGATCTATTTCGTGTATCCGTCTACTGGCCAGAGAAACTCATTACGCCATCTACGAAATAGTCTTCT
CAACCAAGAAGATTCCATACTTCCTCAACTCCTCCATTAGCGGATCATTAATCTTACTGCTCATCGGCGC
CAAAATACCCTTCTCAAAAATCGTACCATCAAGCACCTGCTTCACAACAACACCACATAAAACACCGACC
AACTTTACCATAGCCAAGTACCCCTTCAGATCACCGTACTCCACCAGCGTACTAGTCCTCGTCTCCTTGC
TCCCATCCTTGTTCTCGATCTCGAACTTGTGCTGCAGCATCACAAAATCGCGCTCGTCTTTCTCGAACTG
CATCTTCTGCTCTAGTGTAGCGCACAGAGTGTCGATTGGTTGCCACGGGGGATGATTTTCTCGTCAGAGA
AGATCCCGACCCAGCGGAGACCGGCGAGGATACGCTCTTTCTCTTCGGCGTCTTTGAAGTAGGCTTTTGA
GGCAATGGAAAATTTGAGGTCTTCCTCTGTGGATGAGCTTGTTGCGAGGATCTTCTGCGTCGCTTCTTTC
CAGGCAATGGGTTCTTTCAGGAAGGGCTGTATTTCGTCGGAGAGGAATCCCATGTCAACTAATACACGCA
CAAACTCGGGGAAGCCTTGGTAGCGGAGGGTGCCGCGGATAATGGTTTGTGCTTCGGGAATATTGTAGCG
CTCTTTGTATGGTGTGGAGTCGCGATTCGGATAGGCGACAAATGCGTAGCCAGGATAGATGAAGTATAGT
TTTGCGGTGCCCATCAGATCTTTTCCGGCAATGTCCACGACCTTGCCATCTTTGTAGTATTTTGCTGCGT
TTCTCAACGCGAGTAGAACTCCGCATGATGACCAGGAGAATTTGTAGCCGAGAGGGTTGTCTGAAGCTTC
CGGAGCTGGGAGTCCACCGCAGTACGACAGGAATGAAAGGATATTTCCGCCGGCCTTGTGCACCTCTTCA
ATTGTCTTGATGGCGTACAAGTGGTCGATTCCAGGGTCGAGACCAATCTCGTTCATGACTGTGATGCCGG
CCTTCTTGGCCTCCTCATCCAGTTCCATCATGGCGTCCGAGACGTAGCTGGTTGTCACCACATTCTTCTT
GTTTCTGATGGCGGATTTTATGGCTGTTGCATGAAAGGTATATGGGATCAAGCTTATGACCAAATCATTT
TTTGCAACTTCGGCATCGAGAGCTGCTGCGTTGGTGACGTCGAGAGAGATGGGGTGGGCGTTCTTAACAC
CGGCAGAGAGCTTCTTTGCTGATTCGAGGGTTCTGCAAGCTGAAATGCTATTAACAGGATGTAACGAAGG
GACTGCAATGTTTGGCCATACCAACAGAGACCTCAATGCCCGCATCACTGAGAATGTCAAGGTTGTCGTG
TCACGAAGCCTGCGCCAAGCATAAGAACTTTCCTAAAGAATTGTTAGACAAAATCGGGGAACTTCTGGAG
TTAATTTTCTTACTTTGATCCCGCCATGGTTGCTATTCTCGGATGAGTAGAATTGAAAGTGGAAATCAAT
CAAATAAGGTGCTAAAAAGAGATCAAGACCCCTCACCGAAAGACTCAGACGGAGTGAAAGATTTTGGTGG
GGTTTTGCATGGCTGCGGGCCGTCATCGGGGCTTTTTTGCGCTAACTCCAATCAGGTTAGCAAGGTACCG
AAAATGCTGACTGATCTAAACCGAGGAAAGAACACGCTTTCAAAACCCTTCGACTTCTGCCTTCTTCACC
TCGCGAGCTTCCCAAAGCTTCGTCGTCACATCACCATCACCCAAACGCAATCATGTCGGCACTGGATATG
TTCTGGGCAGCTCCGCCCATCTCGAGGTTTGTCTTGCCGTACCCCGCACTTCCTCATTGGTGCTGATACT
GATTCAGAGCAGGACGTTGGCTCTCTCGGCCTTCACCCTCTCGATTGCGGTCTATACGGGCATGTTGTCC
GGATATTATGTGATATGGCATACTCCATATATTTGGCAGTTCCCACCACAAATATGGAGATTGGTCACAT
CATTCCTCATAACGGGAAAAGACCTGAGTATCATATTCGATACGTATTTCTGTATGCATCTCCTTGGGAC
TGTGCCGCATGTCAGCAACCATTAACATTTGCAACAGTGTACACATATGGGAGCAAGCTTGAGACCGCTT
CTCCAAGGTTTTCGCAGCCAGGAGATTTCATGACCTATCTTCTCTTTGTCTGCTTCACAATATTGGTAGG
CATATGTTCTCTTCTCGTTCGCTTCATTTCTTTGGCACCTCACATTATCTGCCCGCATAGTGCTATTCAC
CCAACTATTGCGGTTCCTGGAGATGAGGAAGATTACCCTTGCACTTCGGCTGGCCCGTCATTCGCAATCA
ATCCAAAGGGTCAGTTCGCGGTGTAGGCATGGTGGGAAGTTCCTTTGAGAAAACTCGTGTCGATTCATCA
ACTTCACTCGGGTGGTTGTGACTCCCTACATTTCATTGATCTTTGAATGTATTCTGTTTCATCTGGTCGA
TGGATCAATCATAAAAATGTCTAACAAACAGAATAGGGCTTGAACATCTTCGTCACAGGTGGCGTGGTAT
TTACATCTTGCCTTGTGCTGGCCTTTGCATACACATCGACCCAAGATGACCGAGGCATGAAGGCCAACTT
CTTCATCCTCACAATCCCTGCACAGTGGATACCATATGCCATGCTTCTTATGACATTGGTCACAAATGGA
CCAAGCGCAACAAAAGTTCAGGCAACCGGCCTTATCGCAGCTCATCTCCACGACTTCTTGACACGACTGT
GGCCAGAATTTGGTGGTGGCTCAAACCTAGTTCTGACACCAACTTTCATCAGAAACATTTTTGCAAAGCC
GCATGCTACAGTCGAGCACAAATCGTATGGTACTGCATTTGCACAGCCACAGCGAGCGGGCTCAGGCTCT
GATAATGTGCTTCCAGAGTCGTGGAAAAGTAGGGGCTCCGGGCACAGATTGGGGGGGAATTGATTCCTTC
GTTGCATTGCAAATGTATACGACGTAAAGAATTAGAAATCTTGTCTATCCAAATTCGTTTCCCTCTAGCG
ACACAAGATGTGACACCAATAGTAGAGCACCGACAAATGCTTCAGCCAGGTGCTAATGTCAATGAGCACT
AATTGTCCACAATTCAGCAGAGGGAGCACAAGTGATGGATCCTCGGATCTTGTGGCCTGCCCCTTGCATG
CTTAGTCACCTTTTCATATGACATCTCGGGTTCAATCCTCAGATTGGCGATCAAAACAATCGACGCCAAA
ACATACCGCACGGAAGGGCATGCTTCAAAAGAGTCTACCAATTGTCATATTTGAGAGCCCCGGATGCGGA
GAAGATAGCCTGGGAGGGGTTGTAGTATGTACAGCAGCTACCTCTACCTCGGTCACGCTTCAGAACATGC
CAACCCTTGGTTCCGAGCTCCAAAAACAGAAGAATATTTCCGAAGACATTCAAGTGCCGTCGAAAAAAAA
AAATCATGAGACCCTCAATCTGGACCGATACATATTACCTACGTGTCCAATGGGTCTCGGCTTGCTCGAA
ATACAGACAGCCTCCTCGGCCACGCCACCCACCCAGTCGGGTGACATTTACCACTCAACGGCAAGCGTCA
ACATGGACTGAATGGTCTGAAGAGGACCGCACAAGAATCCACATTCTGGGTATAAACAGTTGCGGCATGC
TATATGGGCACTCACTTGGAAACCTCCAAGCACGTCCACCTTTGACATATATAGTGCCTTCGAACTGGCG
TGTCATCCAACTCGCCTCGTTGAACGGCGAGATTGAAGTCATCGAAGGTGATGATCACTACATAACAGGC
GGGTTTGATATTGAGGCCATGGAAAGGCCGCCCATCGAGTCTATTTTCGAGCGTCCATCAAGAAAGCACC
TGTCATTCAATACCTTGACTCTACCAGCGAAAACAAAATTCGACATTAGGACCGAACCTCAACCACAAGA
TTCATACTTACCTTTCAGAGGCTTTTCTCAATACCACTATGGGAAACCCATTAAATTCAATAATGACCCG
TCCGAATTAGGATATGGCTCAGATGGGAAGCTTGTGGTTGAATCCAGGCCTCATGAGCGGAAAATAAACA
CATCTGATCTCTTGGTCCTGGATAATGAGGATGCGAAGAAAACATTGGGTCCTATCAAAAATCTCGTATG
TGCAGTCGACGCACAGGTTGTGGTTCACGCTCTGCAGTCGCAAAAGCACAGACTCAACTGCGATTCGACT
ATTCTTTTCACTCAGAAGGGCATGAGAATTATGGAAATGGTAAATGAGCGAGTCTTCCCAGATCCAAAGA
CTCGACCGACCTACATCCCTGGTATTTTCTCTCATGCTGTTTGGAAACCGGACGATCACGCTACTAGTCC
TGATTTCGCTCATCAAATGGAATCAGGATCATCTGGCCTTGAAGCAAGTGTTCCTCCAAAGAAATTGTCC
GTGAAACACTCACTTTTTGGGAGTCTAATCCTTGGTCCGGTGGTCCTGTTGGAGGAGAAAAAGTCATGCA
GAGAGTGCGCCGCCAGCAGAGTGCAAATTATTTCGTCAGTGCTCTTCTTGCAGCCCCTTCCCTCAGGGCG
CGCTGCTTGCCAGCGCATTTGCTTCTCCGCATCAGACTTAGAGATTTAGCAGTCAATTCTGTTGTAGGCC
CGTCCACCGTTCGCTTTCAGTGTAATAATGGCGGCCTCCTTGCGAATGAAGAGCGAATCGCTCACCTCCA
ATCTCGTCTGAGAGAAGCGGCTAGAATCGTCCAAATATTCGATTCATCTCTTACGTACGAACATCTTGAA
AAGCGACTCGGAGATTATCTTGTCATGTCGGGGGATAGCGTCAACATCCTATTCAAAAACGTTCTCGCAG
GAAGAAAGTCCGATATTGAGGTATACTGATTCGATTTCTCTTCACCAAAGTCTATTTCTACTCTGCATAT
ATAATTTGCATTCCTCTTTTTTTTTAACTGCATGTTGCAGAATTTCTTTAGGCAGAGACTATAGTACTTT
TTTTTTACAGATGTGTATAGATACTAATTGTTTCCGCAGTGGAATAGCGGCTGGCTAGTGGAGCAGGGCA
GAATCCTAGCAAAAGGTGGCGCCAGTTTTTGCCCAGGGCACGAAGAAGACTTAGGCTTAATCCAAGAGAT
GACAGAAGCAGCTGCTATTACATGCTGAGAAGTTCGAGGCAGAGCGCCTGGTTAGGACAAAGCAGCGAAA
AATACAACGGGCCAAAAATTGCAAGAAAGAAGTTCCAAATCAGGTCGAGAGAGATCAGGAGGAACTCGCA
ATTATTTTGCATAATCAAAGAATGAAAGAGGAGAGGAATAAGAGGAAGTAGCAGAAACATCTAGAGGCGA
TACAGGAAGTCATTGTTGGCTCGGATGGCGGATCAAGTTCATCTGGAATGAACACAGTAAAGGAGGTCAT
TCTTGATCCTGGGACGAACGAAAATGACATCTCTGAGGATGCTGCTTTCGAATCGCAGGATTCTCCTCTC
AATAGTGGGTTCAAGATCGGGAAGGATTTGCCTGGCGTCGAAAGCCACATTTCTTTCCCCAACGCTTCGC
ATCCAAAGAGATGAAGGATGTGTCTCTACCTCTTGGGGCCAATGCGGGCTCGAATACTGTCGAGCCTGTC
ATGTGAGGCTGCGCTCAAAACTCGGGTTGACTGAAAGACACAGCTCTCTTCTCCGATTTGGTGGGAATCG
GACTCAGAAATTTCGGTTAGACGGCATCCTCACGAGCTATTCGACTACTAGTGGTAAGTCGAATCTTCTT
TGGCTATTCTTCAAAGATAGCAGAAACTTTGAACCCATCACGGGAGACTCATCTCGTATTTCAAATTGCG
TGAAAGGTGTTGGCCTCTCAGGATTCCCTCTGCCTCTGGGCAGAATCTTCGCTCCTGTCTTAGGCCGTCG
TGTTAGCTGCCAAGCCATGGTCCGGAACAGGTTCCCACCAAATTTGGCCTCAATTCCTCTTATGCCATCA
AGATCATCCCACAACCCTTCATTTCACTCAAAATATATCAAGAAATTGAGAAAAACAACTCGATGATGTA
CACACTTTGAACTAAATTCGTCATTTTATCTTACAACTATATTACAGAACATCCATCATGCCCAGCACAA
AGTATAGTTTACAGAATGGGGTAAGAAGAAGATATTTTTGGATTCCCTTTTCACTTCAACAATCAAAACC
GTAAACTCCCTCCTCTGCTTCATGCTTCCAACTGACCACATTACTGCTGCTTCCCAGGTGGGGTATCGAC
TACCCAATTCTCAATCCCAAGCAGAAATCTAGCCTCTCAGATCCTTTGGCAGACACTCCCTCGACAAGGT
CCACGGCCAACCCCGGGCAGCAAACACGCTCTTCTTGCTCTGGGACGTCGACCTCAGACTCGAAACCAAT
AAATACCAAACTTTCACCAATGACCCCTCCGCCGCACTCGATACCCGTTCCCAGGCATTCCCGTAGTAGC
AACCGTGGCCGTCGATTCCATTTCAGGACAATCTGTACTCGCTGCAGCCTGCAAGCTCACCGATCCCGAC
TCCGTCCTCGTCCCCGTGAGTGTGAACTCCGGTGCCGGAGTAGGTGGTAAAGTTGGAAGTGGTATCAACC
CCGTATTTGTCAACGTGAAGAAAGGGATCTTCGTCAAGGTCGTTGACTTAATATTTGCTGAACCGATTGG
ATGATGGTGATGAGTGTGGTTTCCATATTTTGGATAGCCTGTTGGTTGAGGAAAGCGACTGCCGGTCCCG
GTACCGTGACCAGGGAAAGGGAAGAAGATGCCGGAGGGGATGCCTTTTTCTGGTGCGAAGGATGGATGAG
GGAGAGAGGATGCGGGAGTGAGTGGGACGGGAACCAGGGTGATGGTGCTTATTATGACTACCGGTGCTGG
GGGACTTCCGAGCGTCAATTCTGATGGAGAGGCTAACGTTGCTGTGGAAGTTGAAGTGGACTGGCCGGCT
GAGACGGTCGTTGTGCGTCGTTGGATTCCAGCTGATACCATATTCAATAAGAGAAAGGCTAGCAAGGCAA
ACGACAAAAACGTGCGCATTTTGATTGATTGAAAGGTTGAACTCAGAGTCGAAGGGAATAGAACGAGGGT
AGAATGAGTGTAAAAGTCGAGGCCAGCCTGTACCTGACAAGTCGGGCAATTAAAGGTCTGTTTGGTGCAG
AGAGATCCAGCTGAGGAGTAATTGCAATTCTCTCTGGAAGGGCTGAAGTTTTTGATAAATCGCGAACACC
GAGGGGTTTCTGATGACTCAAGAGTTGAAGAACATGGATAAAAGCGAATGTGAATGGAAACTTGTCCAAA
GAAGGGGATCCAAGCATCGGAAGGGAACGAATGGACAATAACAAATCACAAACCATGAGAGACGACCAGA
AGGAAAAAGAGAGGAGGACAGCATCCCGCCAGTGCACGGCTACAGTCTGCCATATAGCGGTTGAGCATAC
ATGAAACTTTGATGGCAGGAGTTAGGAGAGTGTGACGGTTGTTATTTCTTCATTGCAATTGGCGGGACGA
TAGTTCGTTTCAAGGGACGATGTTTCGAAACACCTCAACACCTCCTGAGAGCCAATCTCTCGCTTAGCCA
AAAGCTTGTGCAGACGCCCCGGAAGCACCTCATCCAATCAAGATGAATGTCCAAGGTTCCGCCAATCATA
GTGCCCAGTACATATGCGCACTCAAGCAGTCGAAATACATGGAATGTTGAGGGGATCCTCTGAGCCTGGA
CGAGATAACGCAGATGAGGGTGTGACTTGAAGAAGCCCCTGACCACCAAGAGTATCGTCAACCAGGTCTT
ATTCACGTGGCCTTTATTTAGTTGACCTACAGGTTAAAACCCGCTGAAGGTATATTTAATCCATCCATTC
AAATGTTGGAAGAATGATTTCTGACTTCCTGTTCAGTGCTTTTGAGATCTCATGCCCTGTCATTCTCCAT
AAACTACATTTTAGAACTTCTCTGCAGCAGGTCTCATATCCAGCTCGTGCGTCCACAGCTCAGGCTCTTG
CTGAGATAGCCATAAATACGTCTTCCCAACTGCATCAGCACTCATCTTCTTTCCTATCCTTTGATCTTCT
CCCTCCTTATCCTCAATGGCTCCATTCGCGATGGTATGAACGACATGCACACCCTTCGCGCTCATCTCTC
TCGCCAGAGTCTGCGCTAATTGCCTAACACTCGCTCTTCCAGCACCATATGCTGCATACTCGGCACTGCA
ACGTAACGCTCCCAGAGTTCCCGTAAAGATCAATGTCCCCTTCTTCTCCCCCGTCTCCGAAAGCGGCGCC
TCGCCATGATGCTCAAAGAACAACTTCAACGCTTCCTGAGCAAAGGTGAATGCACCGCCGACGTATTGTT
CTAAACTTTCCGTGAACTCTTCATACGTTTCATTCATGAATGGCTTCTTGCTGGAATGCTTCACGGAGTA
GACTGCTACTTTGAGCTTTAGATCTTTGAAGGAGCTGTGAGATCTGATGGATTGGAAGGCGCTTTTCAGG
TTTTTGGGTGAAGTATCGGTGGTGAAGGTTTCTATGATGCTGCCTGGGGATGTCCTTGCGAGGTTTGAGG
CTACACTGGAAAGGTTCTCCTGTCGTCGGGCGAGGAGGGCTACGGCCATATTTCCATGCTTGGGGTGGGA
GAGAATGCGAGCGATTCCCGTGCCCTGTGAAAAGAATGAGTTTTGAAATTTGTAGGGTACAGATATTGCA
ATGCCTGGTATTCTTATAACATAACATGACTCTTTACGGACCTTTGAACAAAGAAAATGAACTCACAGTG
TTAGGCCCTGCGCCAATGATGATTGCAAGTCCGTTAGGCGGCATATTATCTATCTTATTTAAATCTTAGT
AAAGAACCGACGATACCAAGTGAAGATATAAATCTGTGAAGTTTAGAAGAAAAATCAAACCACCTACCCC
ATATTTTCACGGATTTGGATATCAAGGCTATGACGTCTGCATACCGAGGTACAGTAGTGGGCTCGCGGTG
GACAAAAAAAGCGGCAGTAATTAGCGATCTCACCCAAACTCAGTGTTTTGCATCGCCACCTACGGTATTC
AACGCAACAAACTTCAAAAAATAAACCAAAAATATATCCGCATAATTCAGAAGCTCTAGATCAATACGTC
TGACTATTAAAAGATTACATGGCGTTGGGCTGTACTGATACGGTCGTAAAGAGGGGATTCATCGTGGTGT
CCTGAACTTTGGTGAGGGGGCCACAGAGTGAGCGGCCGCTAATTACTGCCGCTTTTTTTGTCCACCGCGA
GGCGGGGTGGGGTTGCCCGGAGTAAAATCCGACCGAAGTTTCTCCGGTGGACAACCCGATTTCTATTTGA
ACCACACATCACAGTCCACCAACATTATTTGAGCACAAGTATCATTATCCGACATCATCATAGCATTCTA
TTGCATTTCAGTTTTTGCTCGAAAGTCTATATGGTGATCATCAAATTCTAAGATTCTGAAAGAGCTCCAA
TAATTCCCATTTGGAACTTGAGATCCCGTCTGATCACCGGTCTCCACAATAAGCTGCAGCTGTACAACTC
CACTGATTCATCGTCTTGGTCAAAGTTATCGATACTATAACACCTTTCTCAAGCTTATCATCGACGCAGA
CTGCCTGTTCTCGCACCAAGATGCGTCTCGCTGCATTATCAGATACCACCAATCCGTCCCAAGGTTCCAG
CTTTAAAACGCCAATTCGGGAGCCAACCTTCATCTATAAGCCTCTCCCCGAGTTCCACCATTTTGCCGAT
TTCCCAGCAGAACTGCGTCTCATGATTTGGAGATCCGCAATACCAGTAGATTTTAGCGATAATCCAATCA
GTGTCAAACCTTACTCGTACACCTCTCGTAATCGACGAATGCCACCCAAAACCTCGGGAATCAACCAAGA
AAGCCGCAAGGAAACGCTTAAGTATTTACGTCGTCTGGAGCTTCAAACACCTGAGGTATCTGAACACGAG
GGTCACAAATGGCACCGAAAATACCCTATCCATTTTATCCTCTGGGATATTTCCATCGACGTTTTGACCC
TCGATTGGTGGTCCCTAGTCTCCAACGTTACAGATGTCCTCAAGTATCTGAATCTCGAAACTTTCGTCAA
CGAGAAGAAAGGGACCATTAGCTGTGTTCAAATCCTCGACATCGATACTCAAGCCTGGAGTTCCGGTTAT
GGAACCCTTTCCCTAGGACCGCATGTTATGGAACCCATCGAGAAGAAATTAAAGTATTTTGTGAGTTTAA
GGGAGGTGCGACTTAGTATGTTTTCTCGGAATCGGAGTTTTGATGACGAGATACATGTTAGAAAGTCATT
TGAGGAGTGCTTTGTGAGGCTTGCAGCGCAGAATTGCGGATATTGTGTTCCCAGGGTATATTTGGCGCCG
CAGGAATGGAATGATCCTAGAAGACTGATGAGATGAGCTGGCCTATGCCGAAGATGGCTAAAGAGTTAAG
AGCGAAGTTGCTTGAACCTAGGCTTTGAAATATTGTCATTGGCGAAGCAATTTGGGAACCCTAGATGCTT
TGGTTTACATTGTACAGGACATTGGATGCCTACTGTATTTGGATATTTTTCATGTCAAATGTTCTAACTC
TAACGGGACATTCAAATGCCATTTTGGTTGGATGTTGACGCCCGCATGGATACAATGGAATGCACAATAA
TCTTGAAAGTACTTGAACTGTTGGGCCATTCTTCTCTTCCCTGAAATTATTGGGTTATCTCTATCTTTCT
FILE

    print $fh $str;
    
    close($fh);
    
    return 1;
}
