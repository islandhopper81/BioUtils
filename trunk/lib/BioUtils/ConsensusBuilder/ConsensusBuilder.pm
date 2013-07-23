package BioUtils::ConsensusBuilder::ConsensusBuilder;

use strict;
use warnings;

use BioUtils::ConsensusBuilder::FastqColumn 1.0.5;
use BioUtils::ConsensusBuilder::FastqConsensus 1.0.5;
use BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int);
use Carp qw(carp croak);
use version; our $VERSION = qv('1.0.5');
use Exporter qw( import );
our @EXPORT_OK = qw( buildFromClustalwFile buildFromSimpleAlign );

###############
# Subroutines #
###############
sub buildFromClustalwFile($$);
sub buildFromSimpleAlign($$);
sub _parseAlignedSeqsFile($);





###############
# Subroutines #
###############

sub buildFromClustalwFile($$) {
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
    
    for (my $i = 0; $i < $alignmentLen; $i++ ) {
        my $col = BioUtils::ConsensusBuilder::FastqColumn->new();
        
        foreach my $id ( keys %{$alignedSeqs_href} ) {
            my $base = $alignedSeqs_href->{$id}->[$i];
            if ( $base eq "-" ) {
                $dashCount{$id}++;
                $col->addBase($base, int_to_illumina_1_8(0));  # adding a dash
            }
            else {
                my $qual = $quals_href->{$id}->[$i - $dashCount{$id}];
                $col->addBase($base, $qual);
            }
        }
        my ($conBase, $conQual) = $col->getConBaseAndQual();
        if ( defined $conStr and defined $conBase ) {
			$conStr .= $conBase;
			$con_quals_str .= $conQual;
		}
    }
    
    my $fastqConsensus = BioUtils::ConsensusBuilder::FastqConsensus->new({
        seq => $conStr,
        quals_str => $con_quals_str
	});
    
    return ($fastqConsensus);
}

# DEPRECIATED -- see POD
sub buildFromSimpleAlign($$) {
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

sub _parseAlignedSeqsFile($) {
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

This documentation refers to ConsensusBuilder version 1.0.5.

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
    buildFromSimpleAlign
    _parseAlignedSeqsFile
    
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




