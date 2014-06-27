package BioUtils::Align::Pairwise::Alignment;

use warnings;
use strict;

use version; our $VERSION = qv('1.0.11');
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use Scalar::Util qw(looks_like_number);

use MyX::Generic 1.0.11;


{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new( {
										seq1 => ,
										seq2 => ,
										score => ,
										perc_iden => ,
										mismatch_count => ,
										indel_count => ,
										} ) };
	
    # Attributes #
	my %seq1_of;
	my %seq2_of;
	my %score_of;
	my %perc_iden_of;
	my %mismatch_count_of;
	my %indel_count_of;
    
    # Subroutines #
	sub _init;
	sub get_score;
	sub get_seq1;
	sub get_seq2;
	sub get_seq;
	sub get_perc_iden;
	sub get_mismatch_count;
	sub get_indel_count;
	sub get_str;
	
	sub set_score;
	sub set_seq1;
	sub set_seq2;
	sub set_perc_iden;
	sub set_mismatch_count;
	sub set_indel_count;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
		
		# Initialize the parameters
		$new_obj->_init($arg_href);
        
        return $new_obj;
    }
    
    ###############
	# Subroutines #
	###############
	sub _init {
		my ($self, $arg_href) = @_;
		
		$self->set_seq1($arg_href->{seq1});
		$self->set_seq2($arg_href->{seq2});
		$self->set_score($arg_href->{score});
		$self->set_perc_iden($arg_href->{perc_iden});
		$self->set_mismatch_count($arg_href->{mismatch_count});
		$self->set_indel_count($arg_href->{indel_count});
		
		return 1;
	}
	
	sub set_score {
		my ($self, $score) = @_;
		
		if ( ! defined $score ) {
			MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => "set_score(NUMBER)",
            );
		}
		
		if ( ! looks_like_number($score) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Match Score must be a digit",
				value => $score,
			);
		}
		
		$score_of{ident $self} = $score;
		
		return 1;
	}
	
	sub set_seq1 {
		my ($self, $seq) = @_;
		
		if ( ! defined $seq ) {
			MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => 'set_score($seq)',
            );
		}
		
		if ( ref $seq ne "BioUtils::Align::FastaSeq" and
			 ref $seq ne "BioUtils::Align::FastqSeq" ) {
			MyX::Generic::Ref::UnsupportedType->throw(
				error => "Unsupported Type for Seq1",
				this_type => ref $seq,
				supported_types => "BioUtils::Align::FastqSeq or " .
									"BioUtils::Align::FastaSeq",
			);
		}
		
		$seq1_of{ident $self} = $seq;
		
		return 1;
	}
	
	sub set_seq2 {
		my ($self, $seq) = @_;
		
		if ( ! defined $seq ) {
			MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => 'set_score($seq)',
            );
		}
		
		if ( ref $seq ne "BioUtils::Align::FastaSeq" and
			 ref $seq ne "BioUtils::Align::FastqSeq" ) {
			MyX::Generic::Ref::UnsupportedType->throw(
				error => "Unsupported Type for Seq2",
				this_type => ref $seq,
				supported_types => "BioUtils::Align::FastqSeq or " .
								   "BioUtils::Align::FastaSeq",
			);
		}
		
		$seq2_of{ident $self} = $seq;
		
		return 1;
	}
	
	sub set_perc_iden {
		my ($self, $p) = @_;
		
		if ( ! defined $p ) {
			MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => 'set_perc_iden($perc)',
            );
		}
		
		if ( ! looks_like_number($p) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Percent Identity must be a digit",
				value => $p,
			);
		}
		
		if ( $p < 0 or $p > 1 ) {
			MyX::Generic::Digit::OOB->throw(
				error => "Percent Identity must be between 0 and 1",
				value => $p,
				MIN => 0,
				MAX => 1,
			);
		}
		
		$perc_iden_of{ident $self} = $p;
		
		return 1;
	}
	
	sub set_mismatch_count {
		my ($self, $c) = @_;
		
		if ( ! defined $c ) {
			MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => 'set_mismatch_count($count)',
            );
		}
		
		if ( ! looks_like_number($c) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Mismatch count must be a digit",
				value => $c,
			);
		}
		
		if ( $c < 0 ) {
			MyX::Generic::Digit::TooSmall->throw(
				error => "Set mismtach count to greater than 0",
				value => $c,
				MIN => 0,
			);
		}
		
		$mismatch_count_of{ident $self} = $c;
		
		return 1;
	}
	
	sub set_indel_count {
		my ($self, $c) = @_;
		
		if ( ! defined $c ) {
			MyX::Generic::Undef::Param->throw(
                error => 'Undefined parameter value',
                usage => 'set_indel_count($count)',
            );
		}
		
		if ( ! looks_like_number($c) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
				error => "Indel count must be a digit",
				value => $c,
			);
		}
		
		if ( $c < 0 ) {
			MyX::Generic::Digit::TooSmall->throw(
				error => "Set indel count to greater than 0",
				value => $c,
				MIN => 0,
			);
		}
		
		$indel_count_of{ident $self} = $c;
		
		return 1;
	}
	
	sub get_score {
		my ($self) = @_;
		
		return $score_of{ident $self};
	}
	
	sub get_seq1 {
		my ($self) = @_;
		
		return $seq1_of{ident $self};
	}
	
	sub get_seq2 {
		my ($self) = @_;
		
		return $seq2_of{ident $self};
	}
	
	sub get_seq {
		my ($self, $id) = @_;
		
		if ( $seq1_of{ident $self}->get_id() eq $id ) {
			return $seq1_of{ident $self};
		}
		elsif ( $seq2_of{ident $self}->get_id() eq $id ) {
			return $seq2_of{ident $self};
		}
		else {
			croak("Alignment::get_seq() -- cannot find id: $id");
		}
	}
	
	sub get_perc_iden {
		my ($self) = @_;
		
		return $perc_iden_of{ident $self};
	}
	
	sub get_mismatch_count {
		my ($self) = @_;
		
		return $mismatch_count_of{ident $self};
	}
	
	sub get_indel_count {
		my ($self) = @_;
		
		return $indel_count_of{ident $self};
	}
	
	sub get_str {
		my ($self) = @_;
		
		my $str;
		$str .= "######\n";
		$str .= "+" . $self->get_seq1()->get_id() . "\n";
		$str .= "*" . $self->get_seq2()->get_id() . "\n";
		$str .= "\tScore: " . $self->get_score() . "\n";
		$str .= "\tPercent Identity: " . $self->get_perc_iden() . "\n";
		$str .= "\tMismatch Count: " . $self->get_mismatch_count() . "\n";
		$str .= "\tIndel Count: " . $self->get_indel_count() . "\n";
		$str .= "\t+Start: " . $self->get_seq1()->get_start() . "\n";
		$str .= "\t+End: " . $self->get_seq1()->get_end() . "\n";
		$str .= "\t*Start: " . $self->get_seq2()->get_start() . "\n";
		$str .= "\t*End: " . $self->get_seq2()->get_end() . "\n";
		$str .= "\n";
		$str .= "+" . $self->get_seq1()->get_seq() . "\n";
		$str .= "*" . $self->get_seq2()->get_seq() . "\n";
		$str .= "______\n";
	}
}
	
	

1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::Align::Pairwise::Alignment - A module for storing a pairwise alignment


=head1 VERSION

This document describes BioUtils::Align::Pairwise::Alignment version 1.0.11


=head1 SYNOPSIS

    use BioUtils::Align::Pairwise::Alignment;
	
	my $aln = BioUtils::Align::Pairwise::Alignment->new({
				seq1 => $seq1,
				seq2 => $seq2,
				score => $score,
				perc_iden => $perc_iden,
				mismatch_count => $mismatch_count,
				indel_count => $indel_count,
			});
	
	# Set new parameters
	$aln->set_seq1($seq1);
	$aln->set_seq2($seq2);
	$aln->set_score(10);
	$aln->set_perc_iden(.9);
	$aln->set_mismatch_count(2);
	$aln->set_indel_count(1);
	
	# get a seq via it's id
	my $seq = $aln->get_seq($id);
	my $score = $aln->get_score();
	my $perc_iden = $aln->get_perc_iden();
	my $mismatch_count = $aln->get_mismatch_count();
	my $indel_count = $aln->get_indel_count();
	
	# get this object as a sting
	my $str = $aln->get_str();

  
=head1 DESRCIPTION

The Alignment object is a data structure to store aligned sequences.


=head1 DEPENDENCIES

	version
	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	Scalar::Util qw(looks_like_number)

	MyX::Generic 1.0.11


=head1 INCOMPATIBILITIES

	None reported.

=head1 METHODS
	
	_init
	get_seq1
	get_seq2
	get_score
	get_perc_iden
	get_seq
	get_mismatch_count
	get_indel_count
	get_str
	
	set_seq1
	set_seq2
	set_score
	set_prec_iden
	set_mismatch_count
	set_indel_count

=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: 	my $nw = BioUtils::Align::Pairwise::Alignment->new({
						seq1 => $seq1,
						seq2 => $seq2,
						score => $score,
						prec_iden => $perc_iden,
						mismatch_count => $mismatch_count,
						indel_count => $indel_count,
					});
	Function: Create a new NW object
	Returns: Reference to a NW object
	Args: -seq1		=> an aligned sequence
		  -seq2		=> another aligned sequence
		  -score	=> score of the alignment
		  -perc_iden => percent identity of the alignment
		  -mismatch_count => the number of mismatches in the alignment
		  -indel_count => the number of indels in the alignment
	Throws: MyX::Generic::Ref::UnsupportedType
	Comments: NA
	See Also: NA

=head2 _init

	Title: _init
	Usage: $nw->_init($arg_href)
	Function: Initializes the Alignment attributes
	Returns: 1 on completion
	Args: arg_href => arguments hash reference
	Throws: MyX::Generic::Ref::UnsupportedType
	Comments: NA
	See Also: NA

=head2 set_seq1

	Title: set_seq1
	Usage: $aln->set_seq1($seq1)
	Function: Sets an aligned sequence
	Returns: 1 on successful completion
	Args: seq1 => an aligned sequence
	Throws: MyX::Generic::Ref::UnsupportedType
	Comments: NA
	See Also: NA
	
=head2 set_seq2

	Title: set_seq2
	Usage: $aln->set_seq2($seq2)
	Function: Sets an aligned sequence
	Returns: 1 on successful completion
	Args: seq2 => another aligned sequence
	Throws: MyX::Generic::Ref::UnsupportedType
	Comments: NA
	See Also: NA
	
=head2 set_score

	Title: set_score
	Usage: $aln->set_score($score)
	Function: Sets the alignment score
	Returns: 1 on successful completion
	Args: score => the alignment score (INT)
	Throws: MyX::Generic::Digit::MustBeDigit
	Comments: NA
	See Also: NA
	
=head2 set_perc_iden

	Title: set_perc_iden
	Usage: $aln->set_perc_iden($perc)
	Function: Sets the percent identity of the alignment
	Returns: 1 on successful completion
	Args: perc => the percent identity
	Throws: MyX::Generic::Digit::OOB
	Comments: Must be between 0 and 1
	See Also: NA
	
=head2 set_mismatch_count

	Title: set_mismatch_count
	Usage: $aln->set_mismatch_count($count)
	Function: Sets the number of mismatches in the alignment
	Returns: 1 on successful completion
	Args: count => the number of mismatches
	Throws: MyX::Generic::Digit::MustBeDigit
			MyX::Generic::Undef::Param
			MyX::Generic::Digit::TooSmall
	Comments: NA
	See Also: NA
	
=head2 set_indel_count

	Title: set_indel_count
	Usage: $aln->set_indel_count($count)
	Function: Sets the number of indels in the alignment
	Returns: 1 on successful completion
	Args: count => the number of indels
	Throws: MyX::Generic::Digit::MustBeDigit
			MyX::Generic::Undef::Param
			MyX::Generic::Digit::TooSmall
	Comments: NA
	See Also: NA
	
=head2 get_seq1
	
	Title: get_seq1
	Usage: my $seq = $aln->get_seq1()
	Function: gets the first aligned sequence
	Returns: BioUtils::FastqSeq or BioUtils::FastaSeq
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_seq2
	
	Title: get_seq2
	Usage: my $seq = $aln->get_seq2()
	Function: gets the first aligned sequence
	Returns: BioUtils::FastqSeq or BioUtils::FastaSeq
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_seq
	
	Title: get_seq
	Usage: my $seq = $aln->get_seq($id)
	Function: gets the first aligned sequence
	Returns: BioUtils::FastqSeq or BioUtils::FastaSeq
	Args: id => the sequence id to look for in seq1 or seq2
	Throws: croaks when $id is not found
	Comments: NA
	See Also: NA

=head2 get_score
	
	Title: get_score
	Usage: my $score = $aln->get_score();
	Function: gets the score of the alignment
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_perc_iden
	
	Title: get_perc_iden
	Usage: my $perc = $alnget_perc_iden();
	Function: gets the percent identity of the alignment
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_mismatch_count
	
	Title: get_mismatch_count
	Usage: my $mismatch_count = $aln->get_mismatch_count();
	Function: gets the score of number of mismatches in the alignment
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_indel_count
	
	Title: get_indel_count
	Usage: my $indel_count = $aln->get_indel_count();
	Function: gets the score of number of indels in the alignment
	Returns: number
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_str
	
	Title: get_str
	Usage: my $str = $aln->get_str();
	Function: gets this object as a string
	Returns: string
	Args: NA
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
