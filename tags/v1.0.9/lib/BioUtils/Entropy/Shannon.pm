package BioUtils::Entropy::Shannon;

use warnings;
use strict;
use Carp;
use Readonly;
use Class::Std::Utils;
use List::MoreUtils qw(any);
use version; our $VERSION = qv('1.0.9');

# Other recommended modules (uncomment to use):
# use MyX::Generic;
use BioUtils::FastaIO;
use BioUtils::FastqIO;

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new()};
	
	Readonly::Hash my %FA_TYPES => ('fasta' => 1,
                             'FASTA' => 1,
                             'fas' => 1,
                             'FAS' => 1,
                             'fa' => 1,
                             'FA' => 1,);
	Readonly::Hash my %FQ_TYPES => ('fastq' => 1,
                                'FASTQ' => 1,
                                'fq' => 1,
                                'FQ' => 1);

	# Attributes #
	
	# Subroutines #
	sub compute_file;
	sub compute_seq;
	sub compute_str;
	sub _compute_file_FA;
	sub _compute_file_FQ;
	sub _compute_file_line;

	# Constructor #
	sub new {
		my ($class, $arg_href) = @_;

		# Croak if calling new on already blessed reference
		croak 'Constructor called on existing object instead of class'
			if ref $class;

		# Make sure the required parameters are defined
		#if ( any {!defined $_} $arg_href->{arg1},
		#	) {
		#	MyX::Generic::Undef::Params->throw(
		#		error => 'Undefined parameter value',
		#		usage => $NEW_USAGE,
		#	);
		#}

		# Bless a scalar to instantiate an object
		my $new_obj = bless \do{my $anon_scalar}, $class;

		# Attributes

		return $new_obj;
	}

	# Subroutines #
	sub compute_file {
		my ($self, $file, $type) = @_;
		
		if ( defined $FA_TYPES{$type} ) {
			$self->_compute_file_FA($file);
		}
		elsif ( defined $FQ_TYPES{$type} ) {
			$self->_compute_file_FQ($file);
		}
		else {
			$self->_compute_file_line($file);
		}
	}
	
	sub _compute_file_FA {
		my ($self, $file) = @_;
		
		my $io = BioUtils::FastaIO->new({stream_type => '<', file => $file});
		my $e;
		while ( my $seq = $io->get_next_seq() ) {
			$e = $self->compute_seq($seq);
			print $seq->get_header(), "\t", $e, "\n";
		}
	}
	
	sub _compute_file_FQ {
		my ($self, $file) = @_;
		
		my $io = BioUtils::FastqIO->new({stream_type => '<', file => $file});
		my $e;
		while ( my $seq = $io->get_next_seq() ) {
			$e = $self->compute_seq($seq);
			print $seq->get_header(), "\t", $e, "\n";
		}
	}
	
	sub _compute_file_line {
		my ($self, $file) = @_;
		
		open my $IN, "<", "$file" or croak "$!";
		
		my $i = 0;
		my $e;
		while ( my $line = <$IN> ) {
			chomp $line;
			$e = $self->compute_str($line);
			print $i, "\t", $e, "\n";
			$i++;
		}
		
		close($IN);
	}
	
	sub compute_seq {
		my ($self, $seq) = @_;
		
		my $e = $self->compute_str($seq->get_seq());
		
		return $e;
	}
	
	sub compute_str {
		my ($self, $str) = @_;
		
		my $e;  # the returned shannon entropy
		my $char_tot = 0;
		my %char_count;
		
		foreach my $char ( split //, $str) {
			if ( defined $char_count{$char} ) {
				$char_count{$char}++;
			}
			else {
				$char_count{$char} = 1;
			}
			$char_tot++;
		}
		
		foreach my $char ( keys %char_count ) {
			# p is the probability of seeing that character
			my $p = $char_count{$char} / $char_tot;
			$e += $p * log($p);
		}
		
		$e = -$e/log(2);
		
		return $e;
	}
}

1; # Magic true value required at end of module
__END__

=head1 NAME

Entropy::Shannon - Computes the Shannon Entropy for a sequence or set of seqs


=head1 VERSION

This document describes Entropy::Shannon version 1.0.9


=head1 SYNOPSIS

    use Entropy::Shannon;
	
	my $shannon = BioUtils::Entropy::Shannon->new();
	
	# Compute the entropy of all the seqs in a file
	$shannon->compute_file("my_file.fastq", "fastq");
	
	# Compute the entropy of a BioUtils::Fast*Seq object
	$shannon->compute_seq($seq);
	
	# Compute the entropy of a string
	$shannon->compute_str($str);
  
  
=head1 DESCRIPTION

	Shannon entropy is a measure of uncertainty.  It can be used to estimate the
	"randomness" of a sequence (e.g. a DNA sequence).
	
	
=head1 DEPENDENCIES

	Carp
	Readonly
	Class::Std::Utils
	List::MoreUtils qw(any)
	version our $VERSION = qv('1.0.9')
	BioUtils::FastaIO
	BioUtils::FastqIO


=head1 INTERFACE 
	
	compute_file
	compute_seq
	compute_str
	_compute_file_FA
	_compute_file_FQ
	_compute_file_line


=head1 DIAGNOSTICS

	None


=head1 CONFIGURATION AND ENVIRONMENT
  
Entropy::Shannon requires no configuration files or environment variables.


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: my $shannon = BioUtils::Entropy::Shannon->new();
	Function: Create a new BioUtils::Entropy::Shannon object
	Returns: Reference to BioUtils::Entropy::Shannon object
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 compute_file

	Title: compute_file
	Usage: $shannon->compute_file($file, $type);
	Function: Computes the shannon entropy for each sequence in a file
	Returns: NA
	Args: -file => path to fasta or fastq file
		  -type => fastq or fasta
	Throws: NA
	Comments: Prints output to screen
	See Also: NA
	
=head2 compute_seq

	Title: compute_seq
	Usage: $shannon->compute_seq($seq);
	Function: Computes the shannon entropy for a sequence
	Returns: Number
	Args: -seq => BioUtils::Fasta or BioUtils::Fastq object
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 compute_str

	Title: compute_str
	Usage: $shannon->compute_str($seq);
	Function: Computes the shannon entropy for a string
	Returns: Number
	Args: -str => any perl string
	Throws: NA
	Comments: NA
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to the author


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

