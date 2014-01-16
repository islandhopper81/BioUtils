package BioUtils::ConsensusBuilder::ConsensusBuilder;

use strict;
use warnings;

use BioUtils::ConsensusBuilder::FastqColumn 1.0.8;
use BioUtils::ConsensusBuilder::FastqConsensus 1.0.8;
use BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int);
use BioUtils::Codec::IUPAC qw( nuc_str_to_iupac iupac_to_nuc_str);
use BioUtils::MyX::ConsensusBuilder;
use MyX::Generic;
use Carp qw(carp croak);
use version; our $VERSION = qv('1.0.8');
use Exporter qw( import );
our @EXPORT_OK = qw(buildFromClustalwFile buildFromSimpleAlign build_consensus);

###############
# Subroutines #
###############
sub buildFromClustalwFile;
sub _parseAlignedSeqsFile;
sub buildFromSimpleAlign;
sub build_consensus;
sub _build_consensus_from_href;
sub _build_consensus_from_aref;
sub _check_seq_count;
sub _check_seq_ref;
sub _check_aln_len;






###############
# Subroutines #
###############
sub build_consensus {
	my ($seqs_ref) = @_;
	
	if ( ref $seqs_ref eq "HASH" ) {
		return _build_consensus_from_href($seqs_ref);
	}
	elsif ( ref $seqs_ref eq "ARRAY" ) {
		return _build_consensus_from_aref($seqs_ref);
	}
	else {
		MyX::Generic::Ref::UnsupportedType->throw(
			error => "build_consensus expects HASH or ARRAY reference",
			this_type => ref $seqs_ref,
			supported_types => "HASH or ARRAY reference",
		);
	}
}

sub _build_consensus_from_href {
	my ($seqs_href) = @_;
	
	# The input here is a hash where KEY => seq ID and VALUE => FastqSeq
	
	# convert to an array of sequences and use build_consensus_from_arr
	my @seqs_arr;
	foreach my $key ( keys %$seqs_href ) {
		push @seqs_arr, $seqs_href->{$key};
	}
	
	return _build_consensus_from_aref(\@seqs_arr);
}

sub _build_consensus_from_aref {
	my ($seqs_arr) = @_;
	
	# the number of sequence from which to build a consensus
	my $seq_count = _check_seq_count(scalar @$seqs_arr);
	
	# NOTE: right now I assume that the seqs_arr is square for both sequences
	# and quality values
	my @seqs_matrix;
	my @quals_matrix;
	my $alignment_len;
	my $i;
	for ( $i = 0; $i < $seq_count; $i++ ) {
		_check_seq_ref($seqs_arr->[$i]);
		
		my @base_arr = split //, $seqs_arr->[$i]->get_seq();
		push @seqs_matrix, \@base_arr;
		my @qual_arr = split //, $seqs_arr->[$i]->get_quals_str();
		push @quals_matrix, \@qual_arr;
		
		$alignment_len = _check_aln_len($alignment_len,
										scalar @base_arr,
										scalar @qual_arr);
	}
	
	# create some variables to use when building a consensus
	my $col = BioUtils::ConsensusBuilder::FastqColumn->new();
	my ($base, $qual, $iupac, $j);
	my $c_score = 0;
	my $conStr = q{};
    my $con_quals_str = q{};
    
    for ( $i = 0; $i < $alignment_len; $i++ ) {
		# now for each sequence add the bases to the column
        for ( $j = 0; $j < $seq_count; $j++ ) {
            $base = $seqs_matrix[$j]->[$i];
            if ( $base eq "-" ) {
                $col->addBase($base, int_to_illumina_1_8(0));  # adding a dash
            }
            else {
                $qual = $quals_matrix[$j]->[$i];
                $col->addBase($base, $qual);
            }
        }
		
        my ($conBase, $conQual) = $col->getConBaseAndQual();
        if ( defined $conStr and length $conBase > 0 ) {
			$conStr .= $conBase;
			$con_quals_str .= $conQual;
			
			# calcuate the c_score for this column and add to running c_score
			foreach my $b ( split(//, iupac_to_nuc_str($conBase)) ) {
				$c_score += illumina_1_8_to_int($conQual) *
							( $col->getBaseCount($b) / $seq_count );
			}
		}
		
		# clear the variables
		$col->clear_col();
		$base = '';
		$qual = '';
    }
    
    my $fastqConsensus = BioUtils::ConsensusBuilder::FastqConsensus->new({
        seq => $conStr,
        quals_str => $con_quals_str,
		c_score => ($c_score / $alignment_len),
	});
    
    return ($fastqConsensus);
}

sub _check_seq_count {
	my ($seq_count) = @_;
	
	if ( $seq_count < 1 ) {
		BioUtils::MyX::ConsensusBuilder::NoSeqs->throw(
			error => "No seqs to build a consensus from",
		);
	}
	
	return $seq_count;
}

sub _check_seq_ref {
	my ($seq_obj) = @_;
	
	if ( ref $seq_obj eq "BioUtils::FastaSeq" or
		 ref $seq_obj eq "BioUtils::FastqSeq" ) {
		return 1;
	}
	else {
		MyX::Generic::Ref::UnsupportedType->throw(
			error => "Reference must be BioUtils::FastaSeq or BioUtils::FastqSeq",
			this_type => ref $seq_obj,
			supported_types => "BioUtils::FastqSeq, BioUtils::FastaSeq",
		);
	}
}

sub _check_aln_len {
	my ($aln_len, $seq_len, $qual_len) = @_;
	
	# assign and check alignment_len
	if ( ! defined $aln_len ) {
		$aln_len = $seq_len;
	}
	else {
		if ( $aln_len != $seq_len ) {
			BioUtils::MyX::ConsensusBuilder::SeqsNotSqr->throw(
				error => "Alignment seqs are not square",
				alignment_len => $aln_len,
				seq_len => $seq_len,
			);
		}
	}
	
	# check the quality values are same length as seqs
	if ( $aln_len != $qual_len ) {
		BioUtils::MyX::ConsensusBuilder::QualsNotSqr->throw(
				error => "Alignment quals are not square",
				alignment_len => $aln_len,
				quals_len => $qual_len,
			);
	}
	
	return $aln_len;
}

sub buildFromClustalwFile {
    my ($alignedSeqsFile, $quals_href) = @_;
    
    # Read in the alignments created by clustalw
    my ($alignedSeqs_href, $alignmentLen) = _parseAlignedSeqsFile($alignedSeqsFile);

    # Initialize dashCount to count the number of dashes encounted at each
	# position in each sequence.  This is used to index into the original reads
	# to retrieve quality scores when adding bases to a FastqColumn
    my %dashCount = ();
    foreach my $id ( keys %{$alignedSeqs_href} ) {
        $dashCount{$id} = 0;
    }
    
    my $conStr = q{};
    my $con_quals_str = q{};
	
	# create some variables to use when building a consensus
	my $col= BioUtils::ConsensusBuilder::FastqColumn->new();
	my $base;
	my $qual;
	my $c_score = 0;
	my $iupac;
    
    for (my $i = 0; $i < $alignmentLen; $i++ ) {
        foreach my $id ( keys %{$alignedSeqs_href} ) {
            $base = $alignedSeqs_href->{$id}->[$i];
            if ( $base eq "-" ) {
                $dashCount{$id}++;
                $col->addBase($base, int_to_illumina_1_8(0));  # adding a dash
            }
            else {
                $qual = $quals_href->{$id}->[$i - $dashCount{$id}];
                $col->addBase($base, $qual);
            }
        }
		
        my ($conBase, $conQual) = $col->getConBaseAndQual();
        if ( defined $conStr and length $conBase > 0 ) {
			$conStr .= $conBase;
			$con_quals_str .= $conQual;
			
			# calcuate the c_score for this column and add to running c_score
			foreach my $b ( split(//, iupac_to_nuc_str($conBase)) ) {
				$c_score += illumina_1_8_to_int($conQual) *
							( $col->getBaseCount($b) /
							 scalar keys %{$alignedSeqs_href} );
			}
		}
		
		# clear the variables
		$col->clear_col();
		$base = '';
		$qual = '';
    }
    
    my $fastqConsensus = BioUtils::ConsensusBuilder::FastqConsensus->new({
        seq => $conStr,
        quals_str => $con_quals_str,
		c_score => ($c_score / $alignmentLen),
	});
    
    return ($fastqConsensus);
}

# DEPRECIATED -- see POD
sub buildFromSimpleAlign {
    my ($simpleAlignObj, $quals_href) = @_;
    
    # Because this method is depreciated I moved the unnecessary
    # use Bio::SimpleAlign statement into this method and made it a
    # require/import call so it is evaluated at runtime and not compile time.
    require Bio::SimpleAlign;
    import Bio::SimpleAlign;
    
    
    if ( not defined $simpleAlignObj ) {
        die "Undefined simpleAlignObj\n";
    }
    
    # Initialize the dash count for each sequence in the simpleAlignObj to 0
    my %dashCount = ();
    foreach my $seq ( $simpleAlignObj->each_seq() ) {
        $dashCount{$seq->id} = 0;
    }
    
    # Initialize the consensus string and consensus quality values to empty
    my $conStr = q{};
    #my @conQuals = ();
    my $con_quals_str = q{};
    
    # foreach position in the alignment
    #   foreach sequence in the alignment
    for ( my $i = 0; $i < $simpleAlignObj->length(); $i++ ) {
        my $col = BioUtils::ConsensusBuilder::FastqColumn->new();
        
        foreach my $seq ( $simpleAlignObj->each_seq ) { 
            my $base = $seq->subseq($i+1, $i+1);
			
			# NOTE: In the bioperl SimpleAlign object dots represent gaps
            if ( $base eq "-" or $base eq "." ) {  
                $dashCount{$seq->id}++;
                $col->addBase($base, int_to_illumina_1_8(0));  # adding a dash
            }
            else {
                my $qual = $quals_href->{$seq->id}->[$i-$dashCount{$seq->id}];
                $col->addBase($base, $qual);
            }
        }
        my ($conBase, $conQual) = $col->getConBaseAndQual();
        $conStr .= $conBase;
        $con_quals_str .= $conQual;
    }
    
    my $fastqConsensus = BioUtils::ConsensusBuilder::FastqConsensus->new({
		seq => $conStr,
		quals_str => $con_quals_str
	});
    
    return ($fastqConsensus);
}

sub _parseAlignedSeqsFile {
    my ($alignedSeqsFile) = @_;
    
    if ( ! defined $alignedSeqsFile ) {
        croak "Undefined alignedSeqsFile: $alignedSeqsFile\n";
    }
    
    open (ALN, "$alignedSeqsFile") or
		die "Cannot open " . $alignedSeqsFile . "\nERROR: $!\n";
    
    my %alignedSeqsHash = ();
    my $alignmentLen = 0;  # This will get set by the first sequence
    my $header = "";
    my $seq = "";
    
    # store the lines and remove trialing whitespace
    my @lines = <ALN>;
    chomp @lines;
    
    # parse the first header
    my $first = shift @lines;
    if ( $first =~m/^#(.*)/ ) { $header = $1; }
    else { warn "$0 -- parseAlignedSeqsFile -- bad first line\n"; }
    
    # parse the remaining lines
    my $bool = 0;
    foreach my $line ( @lines ) {
        if ( $line =~ m/^#(.*)/ and $bool ) {
            # Add the previous sequence and header to the alignedSeqsHashRef
            my @seqArr = split( //, $seq );
            $alignedSeqsHash{$header} = \@seqArr;
            
            # Set the alignment length with the first sequence
            if ( $alignmentLen == 0 ) {
                $alignmentLen = @seqArr;
            }
            
            # Check the sequence length to make sure it is the same as the alignment length
            if ( length $seq != $alignmentLen ) {
                warn "$0 -- parseAlignedSeqsFile() -- Alignment with UNEQUAL sequence lengths\n";
            }
            
            # reset seq and assing the new header that we just encounted in the regex above.
            $header = $1;
            $seq = "";
        }
        else {
            $seq .= $line;
        }
        $bool = 1;
    }
    
    # For the last sequence, add the sequence and header to the alignedSeqsHashRef
    my @seqArr = split( //, $seq );
    $alignedSeqsHash{$header} = \@seqArr;
    
    close (ALN);
    
    return (\%alignedSeqsHash, $alignmentLen);
}

1;
__END__

#######
# POD #
#######
=head1 BioUtils::ConsensusBuilder::ConsensusBuilder

ConsensusBuilder - A collection of methods for building a consensus sequence
from a multiple sequence alignment (MSA)

=head1 VERSION

This documentation refers to ConsensusBuilder version 1.0.8.

=head1 Included Modules

	BioUtils::ConsensusBuilder::FastqColumn
	BioUtils::ConsensusBuilder::FastqConsensus
	Bio::SimpleAlign  # DEPRECIATED - no longer required
	BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int);
	Carp qw(carp croak);
	version
	Exporter qw( import );

=head1 Inherit

	NA

=head1 SYNOPSIS
    
    use BioUtils::ConsensusBuilder::ConsensusBuilder;
	use BioUtils::ConsensusBuilder::ConsensusBuilder qw(
		buildFromClustalwFile
		buildFromSimpleAlign
	);
    
    # NOTE: This module is just a collection of methods, NOT a class.
    buildFromClustalwFile( $aligned_seqs_file, $quals_href );
    buildFromSimpleAlign($simple_align_obj, $quals_href);
    

=head1 DESCRIPTION

ConsensusBuilder is a collection of methods for building a consensus sequence
from a multiple sequence alignemnt (MSA).  The current approach is to build a
consensus sequence from the clustalw system command output gde file.  A previous
approach was to build an MSA using the BioPerl implemenation of clustalw.  This
approach was much slower because it simply wraps the clustalw system command.
It is also slower because to do this there are several BioPerl objects that are
created which take large amounts of memory.

New methods can be added with different ways of building a consensus sequence
or using different inputs.  However, all methods should return a FastqConsensus
object.

=head1 METHODS

=over

	buildFromClustalwFile
	_parseAlignedSeqsFile
	buildFromSimpleAlign
	build_consensus
	_build_consensus_from_href
	_build_consensus_from_aref
	_check_seq_count
	_check_seq_ref
	_check_aln_len
    
=back

=head1 METHODS DESCRIPTION

=head2 buildFromClustalwFile
    
    Title: buildFromClustalwFile
    Usage: ConsensusBuilder::buildFromClustalwFile($gde_clustalw_file, $quals_href);
    Function: Builds a FastqConsensus using the sequences from a clustalw
                generated gde file.
    Returns: FastqConsensus
    Args: -gde_clustalw_file => a file generated from clustalw in the gde format
          -quals_href => quality values stored in a hash reference with
                         KEY => clustalwId, VALUE => array reference of quals
    Throws: NA
    Comments: NA
    See Also: FastqConsensus
	
=head2 _parseAlignedSeqsFile

    Title: _parseAlignedSeqsFile
    Usage: _parseAlignedSeqsFile($file);
    Function: Parses the output gde output file from clustalw
    Returns: ($aligned_seqs_href, $alignment_length)
    Args: -file => a file with clustalw formated aligned sequences
    Throws: NA
    Comments: The aligned_seqs_href is a hash reference with KEY => clustalwId,
              VALUE => sequence string
    See Also: NA

=head2 buildFromSimpleAlign
    
    Title: buildFromSimpleAlign
    Usage: ConsensusBuilder::buildFromSimpleAlign($simple_align_obj, $quals_href);
    Function: Builds a FastqConsensus using the sequences in a Bio::SimpleAlign
    Returns: FastqConsensus
    Args: -simple_align_obj => a Bio::SimpleAlign object that contains a set of
                               aligned sequences
          -quals_href => quality values stored in a hash reference with
                         KEY => clustalwId, VALUE => array reference of quals
    Throws: NA
    Comments: DEPRECIATED
    See Also: FastqConsensus
	
=head2 build_consensus

    Title: build_consensus
    Usage: build_consensus($seqs_ref);
    Function: Builds a consensus sequences from 2 or more sequences passed in to the
			  method in either a hash ref or array ref.
    Returns: FastqConsensus object
    Args: -seqs_ref => either a hash reference or array reference with FastqSeq objects
    Throws: MyX::Generic::Ref::UnsupportedType
			BioUtils::MyX::ConsensusBuilder::NoSeqs
			BioUtils::Myx::ConsensusBuilder::SeqsNotSqr
			BioUtils::MyX::ConsensusBuilder::QualsNotSqr
    Comments: This method assumes that the sequences are in a square MSA alignment.
    See Also: NA
	
=head2 _build_consensus_from_href

    Title: _build_consensus_from_href
    Usage: _build_consensus_from_href($seqs_href);
    Function: Builds a consensus sequences from 2 or more sequences passed in to the
			  method in a hash ref.
    Returns: FastqConsensus object
    Args: -seqs_href => a hash reference with FastqSeq objects
    Throws: MyX::Generic::Ref::UnsupportedType
			BioUtils::MyX::ConsensusBuilder::NoSeqs
			BioUtils::Myx::ConsensusBuilder::SeqsNotSqr
			BioUtils::MyX::ConsensusBuilder::QualsNotSqr
    Comments: This method assumes that the sequences are in a square MSA alignment.
    See Also: NA
	
=head2 _build_consensus_from_aref

    Title: _build_consensus_from_aref
    Usage: _build_consensus_from_aref($seqs_aref);
    Function: Builds a consensus sequences from 2 or more sequences passed in to the
			  method in an array ref.
    Returns: FastqConsensus object
    Args: -seqs_aref => an array reference with FastqSeq objects
    Throws: MyX::Generic::Ref::UnsupportedType
			BioUtils::MyX::ConsensusBuilder::NoSeqs
			BioUtils::Myx::ConsensusBuilder::SeqsNotSqr
			BioUtils::MyX::ConsensusBuilder::QualsNotSqr
    Comments: This method assumes that the sequences are in a square MSA alignment.
    See Also: NA
	
=head2 _check_seq_count

    Title: _check_seq_count
    Usage: _check_seq_count($seq_count);
    Function: Ensures there are at least 2 seqs to build a consensus from.
    Returns: $seq_count
    Args: -seq_count => Int representing the number of seqs
    Throws: BioUtils::MyX::ConsensusBuilder::NoSeqs
    Comments: NA
    See Also: NA
	
=head2 _check_seq_ref

    Title: _check_seq_ref
    Usage: _check_seq_ref($seq_obj);
    Function: Ensures that the sequences are BioUtils::FastqSeq objects
    Returns: 1 on success
    Args: -seq_obj => The sequence object
    Throws: MyX::Generic::Ref::UnsupportedType
    Comments: NA
    See Also: NA
	
=head2 _check_aln_len

    Title: _check_aln_len
    Usage: _check_aln_len($aln_len, $seq_len, $qual_len);
    Function: Ensures sequences and quality values are in a square matrix
    Returns: $aln_len
    Args: -aln_len => the current alignment length based on the first seq
		  -seq_len => the current sequence length to compare to aln_len
		  -qual_len => the current quality length to compare to aln_len
    Throws: BioUtils::MyX::ConsensusBuilder::SeqsNotSqr
			BioUtils::MyX::ConsensusBuilder::QualsNotSqr
    Comments: NA
    See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.


=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


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

=cut




