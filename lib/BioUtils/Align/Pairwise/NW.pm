package BioUtils::Align::Pairwise::NW;

use warnings;
use strict;

use version; our $VERSION = qv('1.2.1');
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use Scalar::Util qw(looks_like_number);

use MyX::Generic 0.0.1;
use BioUtils::Align::FastaSeq;
use BioUtils::Align::Pairwise::Alignment;
use BioUtils::Align::Pairwise;
use base qw(BioUtils::Align::Pairwise);


{
	# Usage statement
	# see BioUtils::Align::Pairwise
	
    # Attributes #
	# see BioUtils::Align::Pairwise
    
    # Subroutines #
	sub align;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
		my $new_obj = $class->SUPER::new($arg_href);

        return $new_obj;
    }
    
    ###############
	# Subroutines #
	###############
	sub align {
		my ($self, $s1, $s2) = @_;
		
		# Pointer Code
		# diagonal => 0
		# left => 1
		# up => 2
		
		my @seq1 = split //, $s1->get_seq();
		my @seq2 = split //, $s2->get_seq();
		my $match = $self->get_match_score();
		my $mismatch = $self->get_mismatch_score();
		my $gap = $self->get_gap_score();
		
		# step 1: Inializtaion
		my @scores;
		my @ptrs;
		$scores[0][0] = 0;
		$ptrs[0][0] = "none";
		for (my $j = 1; $j <= scalar @seq1; $j++) {  # initalize column
			$scores[0][$j] = $gap * $j;
			$ptrs[0][$j] = 1;
		}
		for ( my $i = 1; $i <= scalar @seq2; $i++) {
			$scores[$i][0] = $gap * $i;
			$ptrs[$i][0] = 2;
		}
		
		# step 2: Populate the matrix
		my ($diagonal_score, $left_score, $up_score);
		
		for ( my $i = 1; $i <= scalar @seq2; $i++ ) {
			for ( my $j = 1; $j <= scalar @seq1; $j++ ) {
				
				# calculate match score
				if ( $seq1[$j-1] eq $seq2[$i-1] ) {
					$diagonal_score = $scores[$i-1][$j-1] + $match;
				}
				else {
					$diagonal_score = $scores[$i-1][$j-1] + $mismatch;
				}
				
				# calculate gap scores
				$up_score = $scores[$i-1][$j] + $gap;
				$left_score = $scores[$i][$j-1] + $gap;
				
				# Choose the best score and set the apropriate information
				if ( $diagonal_score >= $up_score) {
					if ( $diagonal_score >= $left_score ) {
						$scores[$i][$j] = $diagonal_score;
						$ptrs[$i][$j] = 0;
					}
					else {
						$scores[$i][$j] = $left_score;
						$ptrs[$i][$j] = 1;
					}
				}
				else {
					if ( $up_score >= $left_score ) {
						$scores[$i][$j] = $up_score;
						$ptrs[$i][$j] = 2;
					}
					else {
						$scores[$i][$j] = $left_score;
						$ptrs[$i][$j] = 1;
					}
				}
			}
		}
		
		
		# step 3: Trace-back
		my $align1 = '';
		my $align2 = '';
		my $mismatch_count = 0;
		my $indel_count = 0;
		
		# start at the last cell in the matrix
		my $j = scalar @seq1;
		my $i = scalar @seq2;
		my $max_j = $j;
		my $max_i = $i;
		
		while (1) {
			last if $ptrs[$i][$j] eq "none";  # this is the last cell in the matrix
			
			if ( $ptrs[$i][$j] eq 0 ) {
				$align1 .= $seq1[$j-1];
				$align2 .= $seq2[$i-1];
				if ( $seq1[$j-1] ne $seq2[$i-1] ) { $mismatch_count++; }
				$i--;
				$j--;
			}
			elsif ( $ptrs[$i][$j] eq 1 ) {
				$align1 .= $seq1[$j-1];
				$align2 .= "-";
				$indel_count++;
				$j--;
			}
			elsif ( $ptrs[$i][$j] eq 2 ) {
				$align1 .= "-";
				$align2 .= $seq2[$i-1];
				$indel_count++;
				$i--;
			}
		}
		
		$align1 = reverse $align1;
		$align2 = reverse $align2;
		
		#print "$align1\n";
		#print "$align2\n";
		
		my $aln_seq1 = BioUtils::Align::FastaSeq->new({
			header => $s1->get_header(),
			seq => $align1,
			start => 0,
			end => $max_j - 1,
		});
		my $aln_seq2 = BioUtils::Align::FastaSeq->new({
			header => $s2->get_header(),
			seq => $align2,
			start => 0,
			end => $max_i - 1,
		});
		
		my $aln = BioUtils::Align::Pairwise::Alignment->new({
			seq1 => $aln_seq1,
			seq2 => $aln_seq2,
			score => $scores[$max_i][$max_j],
			perc_iden => 1 - (($mismatch_count + $indel_count) / length $align1),
			mismatch_count => $mismatch_count,
			indel_count => $indel_count,
		});
		
		return $aln
	}
}
	
	

1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::Align::Pairwise::NW - A module for running Needleman-Wunsch


=head1 VERSION

This document describes BioUtils::Align::Pairwise::NW version 1.2.1


=head1 SYNOPSIS

    use BioUtils::Align::Pairwise::NW;
	
	my $nw = BioUtils::Align::Pairwise::NW->new();
	# OR #
	my $nw = BioUtils::Align::Pairwise::NW->new({
				match_score => 5,
				mismatch_score => -1,
				gap_score => -2
			});
	
	# Set new parameters
	$nw->set_match_score(2);
	$nw->set_mismatch_score(0);
	$nw->get_gap_score(-1);
	
	# Align two sequences
	$nw->align($seq1, $seq2);

  
=head1 DESRCIPTION

The Needleman-Wunsch algrithm is a dynammic programming algorithm for pairwise
sequence alignment.  See BioUtils::Align::Pairwise for more documentation on
inherited methods.


=head1 DEPENDENCIES

	version
	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	Scalar::Util qw(looks_like_number)

	BioUtils::Align::FastaSeq
	BioUtils::Align::Pairwise::Alignment
	BioUtils::Align::Pairwise
	MyX::Generic 1.2.1


=head1 INCOMPATIBILITIES

	None reported.

=head1 METHODS
	
	align


=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: 	my $nw = BioUtils::Align::Pairwise::NW->new({
						match_score => 5,
						mismatch_score => -1,
						gap_score => -2
					});
			# OR #
			my $nw = BioUtils::Align::Pairwise::NW->new();
	Function: Create a new NW object
	Returns: Reference to a NW object
	Args: -match_score 	  => the match score parameter
		  -mismatch_score => the mismatch score parameter
		  -gap_score	  => the gap score parameter
	Throws: MyX::Generic::Digit::MustBeDigit
	Comments: NA
	See Also: NA
	
=head2 align

	Title: align
	Usage: $nw->align($seq1, $seq2)
	Function: Aligns two sequences
	Returns: BioUtils::Align::Pairwise::Alignment
	Args: -seq1 => BioUtils::FastqSeq OR BioUtils::FastaSeq
		  -seq2 => BioUtils::FastqSeq OR BioUtils::FastaSeq
	Throws: NA
	Comments: NA
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to author


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
