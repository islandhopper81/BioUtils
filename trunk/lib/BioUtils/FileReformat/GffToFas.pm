package BioUtils::FileReformat::GffToFas;

use warnings;
use strict;
use Carp;
use Log::Log4perl qw(:easy);
use Readonly;
use Class::Std::Utils;
use List::MoreUtils qw(any);
use version; our $VERSION = qv('1.0.11');
use MyX::Generic;
use Bio::SeqIO;
use BioUtils::FastaIO;


my $logger = get_logger();

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new( {
		gff_file => ,
        genome_file => ,
        fas_file => ,
        save_feature => ,
		print_exons => ,
		fasta_type => ,
	})};
    Readonly::Scalar my $SAVE_FEATURE => "CDS";
	Readonly::Scalar my $PRINT_EXON_BOOL => "0";  # set to false
	Readonly::Scalar my $FASTA_TYPE => "nucl";
	Readonly::Hash my %FASTA_TYPES => {"nucl" => 1, "prot" => 1};

	# Attributes #
	my %gff_file_of;
    my %genome_file_of;
    my %fas_file_of;
    my %save_feature_of;
	my %print_exons_of;
	my %fasta_type_of;
	
	# Getters #
	sub get_gff_file;
    sub get_genome_file;
    sub get_fas_file;
    sub get_save_feature;
	sub get_print_exons;
	sub get_fasta_type;
	
	# Setters #
	sub set_gff_file;
    sub set_genome_file;
    sub set_fas_file;
    sub set_save_feature;
	sub set_print_exons;
	sub set_fasta_type;
	
	# Others #
	sub reformat;
	sub _init;
	sub _print_seqs;
	sub _get_seq;
	sub _get_header;
	sub _get_gene_seq;
	sub _get_gene_header;
	sub _get_strand_info;
	sub _format_strand_info;
	sub _compare_plus;
	sub _compare_minus;
	sub _get_start;
	sub _to_bool;

	# Constructor #
	sub new {
		my ($class, $arg_href) = @_;
		
		# Croak if calling new on already blessed reference
        $logger->logdie("Constructor called on existing object instead of class")
            if ref $class;

		# Make sure the required parameters are defined
		if ( any {!defined $_}
			$arg_href->{gff_file},
            $arg_href->{genome_file},
            $arg_href->{fas_file},
			) {
			MyX::Generic::Undef::Param->throw(
				error => 'Undefined parameter value',
				usage => $NEW_USAGE,
			);
		}

		# Bless a scalar to instantiate an object
		my $new_obj = bless \do{my $anon_scalar}, $class;

		# Initialize the parameters
		$new_obj->_init($arg_href);

		return $new_obj;
	}

	# Getters #
	sub get_gff_file {
		my ($self) = @_;
		
		return $gff_file_of{ident $self};
	}
    
    sub get_genome_file {
        my ($self) = @_;
        
        return $genome_file_of{ident $self};
    }
    
    sub get_fas_file {
		my ($self) = @_;
		
		return $fas_file_of{ident $self};
	}
    
    sub get_save_feature {
        my ($self) = @_;
        
        return $save_feature_of{ident $self};
    }
	
    sub get_print_exons {
        my ($self) = @_;
        
        return $print_exons_of{ident $self};
    }
	
	sub get_fasta_type {
		my ($self) = @_;
		
		return $fasta_type_of{ident $self};
	}
	
	# Setters #
	sub set_gff_file {
		my ($self, $gff_file) = @_;
		
		# check if the parameter is defined
		if ( ! defined $gff_file ) {
			MyX::Generic::Undef::Param->throw(
				error => "Undefined parameter value (gff_file)"
			);
		}
        
        # check if the file exists
		if ( ! -f $gff_file ) {
			MyX::Generic::DoesNotExist::File->throw(
				error => "File ($gff_file) does not exist"
			)
		}
		
		# check that the file is non empty
		if ( ! -s $gff_file ) {
			MyX::Generic::File::Empty->throw(
				error => "File ($gff_file) is empty"
			);
		}
		
		$gff_file_of{ident $self} = $gff_file;
		
		return 1;
	}
    
    sub set_genome_file {
        my ($self, $genome_file) = @_;
        
        # check if the parameter is defined
        if ( ! defined $genome_file ) {
			MyX::Generic::Undef::Param->throw(
				error => "Undefined parameter value (genome_file)"
			);
		}
		
		# check if the file exists
		if ( ! -f $genome_file ) {
			MyX::Generic::DoesNotExist::File->throw(
				error => "File ($genome_file) does not exist"
			)
		}
		
		# check that the file is non empty
		if ( ! -s $genome_file ) {
			MyX::Generic::File::Empty->throw(
				error => "File ($genome_file) is empty"
			);
		}
		
		$genome_file_of{ident $self} = $genome_file;
		
		return 1;
    }
    
    sub set_fas_file {
		my ($self, $fas_file) = @_;
		
		# check if the parameter is defined
		if ( ! defined $fas_file ) {
			MyX::Generic::Undef::Param->throw(
				error => "Undefined parameter value (fas_file)"
			);
		}
        
        # this is an output file so I don't care if it doesn't exist
		# or if it is empty.  It will be overwritten if it does exist.
		
		$fas_file_of{ident $self} = $fas_file;
		
		return 1;
	}
    
    sub set_save_feature {
        my ($self, $save_feature) = @_;
        
        if ( ! defined $save_feature ) {
            # set the default
            $save_feature = $SAVE_FEATURE;
        }
        
        $save_feature_of{ident $self} = $save_feature;
		
		return 1;
    }
	
	sub set_print_exons {
		my ($self, $print_exons_bool) = @_;
		
		if ( ! defined $print_exons_bool ) {
			$print_exons_bool = $PRINT_EXON_BOOL;
		}
		
		$print_exons_of{ident $self} = _to_bool($print_exons_bool);
		
		return 1;
	}
	
	sub set_fasta_type {
		my ($self, $fasta_type) = @_;
		
		if ( ! defined $fasta_type ) {
			$fasta_type = $FASTA_TYPE;
		}
		
		$fasta_type = lc $fasta_type;
		
		if ( ! defined $FASTA_TYPES{$fasta_type} ) {
			MyX::Generic::BadValue->throw(
				error => "Bad fasta_type: $fasta_type",
				value => $fasta_type
			);
		}
		
		$fasta_type_of{ident $self} = $fasta_type;
		
		return 1;
	}
	
	# Other Subroutines #
	sub _init {
		my ($self, $arg_href) = @_;
		
		$self->set_gff_file($arg_href->{gff_file});
        $self->set_genome_file($arg_href->{genome_file});
        $self->set_fas_file($arg_href->{fas_file});
        $self->set_save_feature($arg_href->{save_feature});
		$self->set_print_exons($arg_href->{print_exons});
		$self->set_fasta_type($arg_href->{fasta_type});
		
		return 1;
	}
	
	sub reformat {
		my ($self) = @_;

        ### save the genome
        my $genome_in = BioUtils::FastaIO->new({
            stream_type => '<',
            file => $self->get_genome_file()
        });
        
        my %contigs = ();
        while ( my $contig = $genome_in->get_next_seq() ) {
            $contigs{$contig->get_id()} = $contig;
        }
        
        
        ### go through each line of the gff file
        open my $GFF, "<", $self->get_gff_file() or
            MyX::Generic::File::CannotOpen->throw(
				error => "Cannot open file",
				file_name => $self->get_gff_file(),
			);
        
		# I store the sequences in a hash with KEY=> sequence ID,
		# VALUE => array ref of seq objects.  This allows me to paste
		# together sequences having the same ID into a single sequence.
		# If there are seqeunces that share an ID this is because they
		# are seperate exons that come from the same gene.
		my %seqs_by_id = ();  # each ID has an array ref of sequences
		
		# go through each line in the GFF and save the sequences
        foreach my $line ( <$GFF> ) {
            chomp $line;
            
            if ( $line =~ m/^#/ ) { next; } # skip comment lines
            my @vals = split(/\t/, $line);
            my $contig_name = $vals[0];
            my $feature = $vals[2];
            my $start = $vals[3];
            my $end = $vals[4];
            my $strand = $vals[6];
            my $info = $vals[8];
            
            # skip the CRISPR line
            if ( $feature eq "CRISPR" ) {
                $logger->warn("Skipping CRISPR line in GFF");
                next;
            }
            
            # check if this is a feature we want to save and print
			my $save_feature = $self->get_save_feature();
            if ( $save_feature !~ m/all/i and $save_feature ne $feature ) {
                $logger->debug("Skipping line because $feature is not a $save_feature line");
                next;
            }
            
            # NOTE: It looks like the numbers in the GFF file are on a 1-based
            #       start, so I need to subtract 1 from the start and end to
            #       transform them to a 0-based start that is required for the
            #       BioUtils::FastaSeq substr function.
            my $gene_seq_obj = _get_seq(\%contigs, $contig_name, $start-1,
                                    $end-1, $strand, $logger);
            my ($header, $id) = _get_header($self, $info, $start, $end, $strand,
											$contig_name, $logger);
            
            if ( $header ne "0" ) {
                $gene_seq_obj->set_header($header);
            }
            else {
                $logger->warn("Skipping GFF line because can't make header: $line");
                next;
			}
			
			$logger->debug("Saving sequence to output later");
			if ( defined $seqs_by_id{$id} ) {
				push @{$seqs_by_id{$id}}, $gene_seq_obj;
			}
			else {
				my @arr = ($gene_seq_obj);
				$seqs_by_id{$id} = \@arr;
			}
        }
		
		# output the sequences
		_print_seqs($self, \%seqs_by_id);
        
        $logger->info("Finished!");
		
		return 1;
	}
	
	sub _print_seqs {
		my ($self, $seqs_by_id_href) = @_;
		
		my $fasta_out = BioUtils::FastaIO->new({stream_type => '>',
                                                file => $self->get_fas_file()});
		
		my $seq_obj;
		foreach my $id ( keys %{$seqs_by_id_href} ) {
			# remember the value assigned to each ID is an array ref of
			# BioUtils::FastaSeq objects
			
			if ( $self->get_print_exons() ) {
				foreach $seq_obj ( @{$seqs_by_id_href->{$id}} ) {
					
					if ( $self->get_fasta_type() eq "prot" ) {
						$seq_obj->translate();
					}
					
					$logger->debug("Writting seq");
					$fasta_out->write_seq($seq_obj);
				}
			}
			else {
				# print the genes.  This means combining sequences with
				# the same ID (ie multiple exons in the same gene)
				$seq_obj = _get_gene_seq($self, $seqs_by_id_href->{$id});
				
				if ( $self->get_fasta_type() eq "prot" ) {
					$seq_obj->translate();
				}
				
				$logger->debug("Writting seq");
				$fasta_out->write_seq($seq_obj);
			}
		}
		
		return 1;
	}
    
    sub _get_seq {
        my ($contig_seqs_href, $contig_name, $start, $end, $strand, $logger) = @_;
        
        if ( ! defined $contig_seqs_href->{$contig_name} ) {
            $logger->warn("Cannot find contig: $contig_name");
            return(undef);
        }
        
        # there are some cases where the start is greater than the end!!  WTH!
        # at this point I just switch them
        if ( $start > $end ) {
            my $tmp = $start;
            $start = $end;
            $end = $tmp;
        }
        
        my $new_seq_obj = $contig_seqs_href->{$contig_name}->substr($start, $end);
        
        if ( $strand =~ m/-/ ) {
            $logger->debug("On minus strand so reverse complementing");
            $new_seq_obj->rev_comp();
        }
        
        return($new_seq_obj);
    }
    
    sub _get_header {
        my ($self, $info, $start, $end, $strand, $contig, $logger) = @_;
        
        my ($id, $locus_tag, $product);
        my $other = "";
        
        if ( defined $info ) {
            my @info_vals = split(/;/, $info);
            foreach my $val ( @info_vals ) {
                $logger->debug("Val: $val");
                if ( $val =~ m/ID=(\S+)/ ) {
                    $id = $1;
                }
                elsif ( $val =~ m/proteinId=(\S+)/ ) {
                    # this is a second possible way to define the ID value
                    $id = $1;
                }
                elsif ( $val =~ m/locus_tag=(\S+)/ ) {
                    $locus_tag = $1;
                }
                elsif ( $val =~ m/product=(\S+)/ ) {
                    $product = $1;
                }
                elsif ( $val =~ m/(\S+=\S+)/ ) {
                    $other .= $1 . " ";
                }
            }
        }
        
        if ( ! defined $id ) {
            my $message = "ID not defined ";
            $message .= "for info ($info) " if (defined $info);
            $message .= "on contig ($contig) ";
            $message .= "in file (" . $self->get_gff_file() . ")";
            $logger->warn($message);
            return(0);
        }
        
        # the locus tag is not required, but if it is not defined set it to NA
        if ( ! defined $locus_tag ) {
            $locus_tag = "NA";
            $logger->debug("locus_tag not defined: $info on contig: $contig");
        }
        
        my $header = "$id $locus_tag ";
        $header .= "product=$product " if (defined $product);
        $header .= "$other " if (defined $other);
        $header .= "start=$start end=$end strand=$strand contig=$contig";
        
        $logger->debug("Header: $header");
        
        return($header, $id);
    }
	
	sub _get_gene_seq {
		# remove the class (ie $self) string if it is there
		# I do this so I can test the functions by calling:
		# BioUtils::FileReformat::GffToFna->_compare_plus
		#if ( scalar @_ > 1 ) {shift @_;}
		
		my ($self, $seqs_aref) = @_;
		
		# if there is only one sequence return it
		if ( scalar @$seqs_aref == 1 ) {
			# make a new refernece for a BioUtils::FastaSeq object to return
			# because that is consistant with what I do below.  However,
			# note that this wastes some memory.
			my $new_seq_obj = BioUtils::FastaSeq->new({
				header=> $seqs_aref->[0]->get_header(),
				seq => $seqs_aref->[0]->get_seq()
			});
			
			if ( $self->get_fasta_type() eq "prot" ) {
				$new_seq_obj->translate();
			}
			return $new_seq_obj;
		}
		
		# if there are more than one sequence combine them
		# at this point all the sequences are 5` -> 3`
		# but I still need to order them correctly
		
		# first, get the strand info from the first seq
		my $strand = _get_strand_info($seqs_aref->[0]->get_header());
		
		my @ordered = ();
		if ( $strand eq "+" ) {
			# order from smallest to largest
			@ordered = sort { _compare_plus($a, $b) } @$seqs_aref;
		}
		elsif ( $strand eq "-" ) {
			# order from largest to smallest
			@ordered = sort { _compare_minus($a, $b) } @$seqs_aref;
		}
		
		# use the header of the first sequence as a template to build
		# a new header
		my $new_header = _get_gene_header(\@ordered);
		
		# get the new seq by concatenating all the ordered seqs
		my $new_seq_str = "";
		foreach my $seq_obj ( @ordered ) {
			$new_seq_str .= $seq_obj->get_seq();
		}
		
		my $new_seq_obj = BioUtils::FastaSeq->new({header=> $new_header,
												   seq => $new_seq_str});
		
		if ( $self->get_fasta_type() eq "prot" ) {
				$new_seq_obj->translate();
		}
		
		return($new_seq_obj);
	}
	
	sub _get_gene_header {
		# remove the class (ie $self) string if it is there
		# I do this so I can test the functions by calling:
		# BioUtils::FileReformat::GffToFna->_compare_plus
		if ( scalar @_ == 2 ) {shift @_;}
		
		my ($seqs_aref) = @_;
		
		# get the template from the last sequence
		my $header = $seqs_aref->[$#{ $seqs_aref }]->get_header();
		
		# get the start coordinate
		my $start = _get_start($seqs_aref->[0]);
		
		# substitute the start value
		$header =~ s/start=\d+/start=$start/;
		
		return($header);
	}
	
	sub _get_strand_info {
		my ($header) = @_;
		
		if ( $header =~ m/strand=(\S+)/ ) {
			return _format_strand_info($1);
		}
		else {
			$logger->warn("Cannot find strand info on header: $header");
		}
		
		return 0;
	}
	
	sub _format_strand_info {
		my ($strand) = @_;
		
		# converts the strand info into either "+" or "-"
		
		if ( $strand eq "+" ) {
			return $strand;
		}
		elsif ( $strand eq "-" ) {
			return $strand;
		}
		elsif ( $strand =~ m/plus/i ) {
			return "+";
		}
		elsif ( $strand =~ m/minus/i ) {
			return "-";
		}
		else {
			$logger->warn("Unknown strand value: $strand");
		}
		
		return 0;
	}
	
	sub _compare_plus {
		# remove the class (ie $self) string if it is there
		# I do this so I can test the functions by calling:
		# BioUtils::FileReformat::GffToFna->_compare_plus
		if ( scalar @_ == 3 ) {shift @_;}
		
		my ($a, $b) = @_;
		
		return( _get_start($a) - _get_start($b) );
	}
	
	sub _compare_minus {
		# remove the class (ie $self) string if it is there
		# I do this so I can test the functions by calling:
		# BioUtils::FileReformat::GffToFas->_compare_plus
		if ( scalar @_ == 3 ) {shift @_;}
		
		my ($a, $b) = @_;
		
		return( _get_start($b) - _get_start($a) );
	}
	
	sub _get_start {
		my ($seq_obj) = @_;
		
		if ( $seq_obj->get_header() =~ m/start=(\d+)/ ) {
			return $1;
		}
		
		$logger->warn("BAD WARNING: Cannot find start!");
		
		return 0;
	}
	
	sub _to_bool {
		my ($bool) = @_;
		
		if ( $bool eq 1 or $bool eq 0 ) {
			return $bool;
		}
		
		my %good_yes_values = map { $_ => 1 } qw(Y YES Yes y yes);
		if ( defined $good_yes_values{$bool} ) {
			return 1;
		}
		
		# else return false
		return 0;
	}
}

1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::FileReformat::GffToFas - Converts GFF to fas file


=head1 VERSION

This document describes BioUtils::FileReformat::GffToFas version 1.0.11


=head1 SYNOPSIS

    use BioUtils::FileReformat::GffToFas;

	my $gbk_to_gff = BioUtils::FileReformat::GffToFas->new({
		gff_file => "my_genome.gff",
        genome_file => "my_genome.fas",
        fas_file => "my_genes.fas",
        save_feature => "CDS",
		print_exons => 0,
		fasta_type => "nucl"
	});
	
	# run the reformating
	$gff_to_fas->reformat();
	
	# reset the gkb input files and gff ouptut file and other data
	$gff_to_fas->set_gff_file("my_new_genome.gff");
    $gff_to_fas->set_genonme_file("my_genome.fas");
    $gff_to_fas->set_fas_file("my_genes.fas");
    $gff_to_fas->set_save_feature("CDS");
	$gff_to_fas->set_print_exons(0);
	$gff_to_fas->set_fasta_type("nucl");
	
	# get the current setting for the input gff file
	$gff_to_fas->get_gff_file();
	
    # get the current setting for the input genome file
	$gff_to_fas->get_genome_file();
    
	# get the current setting for the output fas file
	$gff_to_fas->get_fas_file();
    
    # get the current setting for the save feature
    $gff_to_fas->get_save_feature();
	
	# A shorthand version
	BioUtils::FileReformat::GffToFas->new({
		gff_file => "my_genome.gff",
        genome_file => "my_genome.fas",
        fas_file => "my_genes.fas"
	})->reformat();
	
	# Print individual exons instead of gene coding sequence
	$gff_to_fas->set_print_exons(1);
	
	# Save the sequences from all the GFF features
	$gff_to_fas->set_save_feature("all");
	
	# output protein sequences rather than nucleotides
	$gff_to_fas->set_fasta_type("prot");
  
  
=head1 DESCRIPTION

BioUtils::FileReformat::GffToFas reformats a standard GFF file into a fas file.
The fas file has headers with the format: 
">[ID] [locus_tag] product=[product] [other_info] start=[start] end=[end] contig=[contig]

The other_info are other tags that are defined in the gff tag line.

Note that the contig names in the genome file must match those in the GFF file.

Sequence from any GFF feature can be saved.  The default is to only save the CDS
feature.  To save the sequences from all the features use set_save_feature("all").
However, there will be repeated sequences that are associated with multiple
features.

For eukaryotes the GFF often describes individual exons.  This script can either
transform that GFF information into sequences for each exon or a full gene coding
sequence.  By default exons from the same gene are concatentated into a coding
sequence.  To keep individual exons set print_exons to true (ie 1) either in
the constructor or by calling the method set_print_exons(1).

When the exons are combined the header is created by using mostly information
from the header of the last exon sequence.  The only change is the start value
which is taken from the header of the first exon.


=head1 DIAGNOSTICS

All errors are pass out of methods as MyX::Generic objects.  See MyX::Generic
documentaiton for infomation about catching and handling those errors.


=head1 CONFIGURATION AND ENVIRONMENT
  
BioUtils::FileReformat::GffToFas requires no configuration files or environment variables.


=head1 DEPENDENCIES

Carp
Log::Log4perl qw(:easy)
Readonly
Class::Std::Utils
List::MoreUtils qw(any)
version
MyX::Generic
Bio::SeqIO
BioUtils::FastaIO


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS
	
=over

	new
	reformat
	get_fas_file
    get_genome_file
	get_fas_file
    get_save_feature
	get_print_exons
	get_fasta_type
	set_gff_file
    set_genome_file
	set_fas_file
    set_save_feature
	set_print_exons
	set_fasta_type
	_init
	_print_seqs
	_get_seq
	_get_header
	_get_gene_seq
	_get_gene_header
	_get_strand_info
	_format_strand_info
	_compare_plus
	_compare_minus
	_get_start
	_to_bool

=back

=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: BioUtils::FileReformat::GffToFas->new({
                gff_file => "input.gff",
                genome_file => "genome.fasta",
				fas_file => "output_genes.fas",
				save_feature => "CDS"
			});
	Function: Creates a new BioUtils::FileReformat::GffToFas object
	Returns: BioUtils::FileReformat::GffToFas
	Args: -fas_file => a file path to the output gene fasta file
	      -gff_file => a file path to the input gff file
          -genome_file => a file path to the input genome fasta file
          -save_feature => Feature in the GFF that you want to output in the fas
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA
	
=head2 reformat

	Title: reformat
	Usage: $GffToFas_obj->reformat()
	Function: Reformats the GFF3 into a fasta file of gene nucleotide sequences
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_gff_file

	Title: get_gff_file
	Usage: $GffToFas_obj->get_gff_file()
	Function: Returns input GFF3 file path
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_genome_file

	Title: get_genome_file
	Usage: $GffToFas_obj->get_genome_file()
	Function: Returns input genome fasta file path
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_fas_file

	Title: get_fas_file
	Usage: $GffToFas_obj->get_fas_file()
	Function: Returns output gene fasta file path
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
    
=head2 get_save_feature

	Title: get_save_feature
	Usage: $GffToFas_obj->get_save_feature()
	Function: Returns the GFF feature that will be output as a gene
	Returns: str
	Args: NA
	Throws: NA
	Comments: Default is set to "CDS"
	See Also: NA
	
=head2 get_fasta_type

	Title: get_fasta_type
	Usage: $GffToFas_obj->get_fasta_type()
	Function: Returns the output file fasta type
	Returns: "nucl" | "prot"
	Args: NA
	Throws: NA
	Comments: Default is set to "nucl"
	See Also: NA
	
=head2 get_print_exon

	Title: get_print_exon
	Usage: $GffToFas_obj->get_print_exon()
	Function: Returns the boolean value for if exons should be printed
	Returns: bool
	Args: NA
	Throws: NA
	Comments: Default is set to 0 (ie false).  By default it will print the
			  gene sequences.  This means the program will concatenate the
			  exon DNA sequences into a single coding sequence.  If you want
			  to print each exon seperately set this to 1 (ie true);
	See Also: NA
	
=head2 set_gff_file

	Title: set_gff_file
	Usage: $GffToFas_obj->set_gff_file("file.gff")
	Function: Sets the input GFF3 file path
	Returns: 1 on success
	Args: Path string to input GFF3 file
	Throws: MyX::Generic::Undef::Param
            MyX::Generic::DoesNotExist::File
			MyX::Generic::File::Empty
	Comments: NA
	See Also: NA

=head2 set_genome_file

	Title: set_genome_file
	Usage: $GffToFas_obj->set_genome_file("genome.fasta")
	Function: Sets the input genome fasta file path
	Returns: 1 on success
	Args: Path string to input genome fasta file
	Throws: MyX::Generic::Undef::Param
            MyX::Generic::DoesNotExist::File
			MyX::Generic::File::Empty
	Comments: NA
	See Also: NA
    
=head2 set_fas_file

	Title: set_fas_file
	Usage: $GffToFas_obj->set_fas_file("file.gbk")
	Function: Sets the output gene fasta file path
	Returns: 1 on success
	Args: Path string to output gene fasta file
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA
    
=head2 set_save_feature

	Title: set_save_feature
	Usage: $GffToFas_obj->set_save_feature("CDS")
	Function: Sets the GFF feature to save as a gene
	Returns: 1 on success
	Args: Path string to output gene fasta file
	Throws: MyX::Generic::Undef::Param
	Comments: Features that can be saved are any string found in the 3rd column
              of the GFF file.  If you want to save the sequences for all the
              features use "all".
	See Also: NA
	
=head2 set_print_exons

	Title: set_print_exons
	Usage: $GffToFas_obj->set_print_exons(0)
	Function: Sets the boolean value to indicate if exons should be printed
			  seperately
	Returns: 1 on success
	Args: 0 | 1
	Throws: NA
	Comments: Default is set to 0 (ie false).  By default it will print the
			  gene sequences.  This means the program will concatenate the
			  exon DNA sequences into a single coding sequence.  If you want
			  to print each exon seperately set this to 1 (ie true);
	See Also: NA
	
=head2 set_fasta_type

	Title: set_fasta_type
	Usage: $GffToFas_obj->set_fasta_type("nucl")
	Function: Sets the type of fasta file to output
	Returns: 1 on success
	Args: "nucl" | "prot"
	Throws: MyX::Generic::BadValue
	Comments: Default is set to "nucl".
	See Also: NA
	
=head2 _init

	Title: _init
	Usage: $GffToFas_obj->_init({
	            fas_file => "output_genes.fasta",
                genome_file => "genome.fasta",
				gff_file => "input.gff",
                save_feature => "CDS",
				print_exons => 0,
				fasta_type => "nucl"
			});
	Function: Initializes the object by setting the parameter values
	Returns: 1 on success
	Args: -gff_file => a file path to the input gff file
          -genome_file => a file path to the input genome fasta file
          -fas_file => a file path to the output gene fasta file
          -save_feature => feature from the GFF file to save as a gene
		  -print_exons => boolean indicating if individual exons should be
						  printed
		  -fasta_type => "nucl" | "prot"
	Throws: NA
	Comments: Do NOT call this method.  It is private
	See Also: NA
	
=head2 _print_seqs

    Title: _print_seqs
    Usage: _print_seqs($self, $seqs_by_id_href);
    Function: prints the seqs in $seqs_by_id_href
    Returns: 1 on success
    Args: -self => the self object
		  -seqs_by_id_href => seqs_by_id hash reference
    Throws: NA
    Commants: Do NOT call this method.  It is private
    See Also: NA
    
=head2 _get_seq

    Title: _get_seq
    Usage: _get_seq($contig_seqs_href, $contig_name, $start, $end, $strand,
			$logger);
    Function: Get a nucleotide sequence from a GFF entry (ie line)
    Returns: BioUtils::FastaSeq object
    Args: --
    Throws: NA
    Commants: Do NOT call this method.  It is private
    See Also: NA
	
=head2 _get_header

    Title: _get_header
    Usage: _get_header($self, $info, $start, $end, $strand, $contig, $logger);
    Function: Get sequence header from info in a GFF entry (ie line)
    Returns: BioUtils::FastaSeq object
    Args: --
    Throws: NA
    Commants: Do NOT call this method.  It is private
    See Also: NA
	
=head2 _get_gene_seq

    Title: _get_gene_seq
    Usage: _get_gene_seq($seqs_aref);
    Function: Gets the concatenated exon sequences as a single coding sequence
    Returns: BioUtils::FastaSeq object
    Args: -seqs_aref => Array reference of BioUtils::FastaSeq objects that
						are all exons of a single gene
    Throws: NA
    Commants: Do NOT call this method.  It is private
    See Also: NA
	
=head2 _get_gene_header

    Title: _get_gene_header
    Usage: _get_gene_header($seqs_aref);
    Function: Gets the header of the concatenated exon sequences 
			  of a single gene
    Returns: Str
    Args: -seqs_aref => Array reference of BioUtils::FastaSeq objects that
						are all exons of a single gene
    Throws: NA
    Commants: Do NOT call this method.  It is private
    See Also: NA
	
=head2 _get_strand_info

    Title: _get_strand_info
    Usage: _get_strand_info($header);
    Function: Gets the strand info from a FastaSeq header
    Returns: + | -
    Args: -header => header string from a BioUtils::FastaSeq obj
    Throws: warning: "Cannot find strand info on header: $header"
    Commants: Do NOT call this method.  It is private
    See Also: NA
	
=head2 _format_strand_info

    Title: _format_strand_info
    Usage: _format_strand_info($strand);
    Function: Gets the formated strand info from a strand value
    Returns: + | -
    Args: -strand => the strand info
    Throws: warning: "Unknown strand value: $strand"
    Commants: Do NOT call this method.  It is private.
    See Also: NA
	
=head2 _compare_plus

    Title: _compare_plus
    Usage: _compare_plus($a, $b);
    Function: Orders two BioUtils::FastaSeq objs based on their start value
    Returns: negative int, 0, or positive int
    Args: -a => a BioUtils::FastaSeq obj with start=\d+ in the header
	      -b => a BioUtils::FastaSeq obj with start=\d+ in the header
    Throws: NA
    Commants: Do NOT call this method.  It is private.  Used to sort a list of
			  exons in order of how they fall in the genome.  This function
			  orders exons that are found on the plus strand.  
    See Also: _compare_minus for a function that orders exons on the minus
			  strand.

=head2 _compare_minus

    Title: _compare_minus
    Usage: _compare_minus($a, $b);
    Function: Orders two BioUtils::FastaSeq objs based on their start value
    Returns: negative int, 0, or positive int
    Args: -a => a BioUtils::FastaSeq obj with start=\d+ in the header
	      -b => a BioUtils::FastaSeq obj with start=\d+ in the header
    Throws: NA
    Commants: Do NOT call this method.  It is private.  Used to sort a list of
			  exons in order of how they fall in the genome.  This function
			  orders exons that are found on the minus strand.  
    See Also: _compare_plus for a function that orders exons on the plus
			  strand.

=head2 _get_start

    Title: _get_start
    Usage: _get_start($seq_obj);
    Function: Gets the start=(\d+) value from a BioUtils::FastaSeq obj
    Returns: int
    Args: -seq_obj => a BioUtils::FastaSeq obj with start=\d+ in the header
    Throws: "BAD WARNING: Cannot find start!"
    Commants: Do NOT call this method.  It is private.  
    See Also: NA
	
=head2 _to_bool

    Title: _to_bool
    Usage: _to_bool($bool);
    Function: Makes sure a boolean is formated as either 1 or 0
    Returns: bool
    Args: -bool => a boolean value
    Throws: NA
    Commants: Do NOT call this method.  It is private.  1 == TRUE, 2 == FALSE
    See Also: NA


=head1 BUGS AND LIMITATIONS

Please report any bugs or feature requests to
C<bug-<RT NAME>@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

=head2 TO DO

Add functionality to concatenate all the Exons from a eukaryotic organism.

=head1 AUTHOR

Scott Yourstone  C<< scott.yourstone81@gmail.com >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

