package BioUtils::FileReformat::GbkToGff;

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
use Bio::Tools::GFF;

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new( {
		gbk_file => ,
		gff_file => ,
	})};

	# Attributes #
	my %gbk_file_of;
	my %gff_file_of;
	
	# Getters #
	sub get_gbk_file;
	sub get_gff_file;
	
	# Setters #
	sub set_gbk_file;
	sub set_gff_file;
	
	# Others #
	sub reformat;
	sub _init;

	# Constructor #
	sub new {
		my ($class, $arg_href) = @_;

		# Croak if calling new on already blessed reference
		croak 'Constructor called on existing object instead of class'
			if ref $class;

		# Make sure the required parameters are defined
		if ( any {!defined $_}
			$arg_href->{gbk_file},
			$arg_href->{gff_file},
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
	sub get_gbk_file {
		my ($self) = @_;
		
		return $gbk_file_of{ident $self};
	}
	
	sub get_gff_file {
		my ($self) = @_;
		
		return $gff_file_of{ident $self};
	}
	
	# Setters #
	sub set_gbk_file {
		my ($self, $gbk_file) = @_;
		
		# check if the parameter is defined
		if ( ! defined $gbk_file ) {
			MyX::Generic::Undef::Param->throw(
				error => "Undefined parameter value (gbk_file)"
			);
		}
		
		# check if the file exists
		if ( ! -f $gbk_file ) {
			MyX::Generic::DoesNotExist::File->throw(
				error => "File ($gbk_file) does not exist"
			)
		}
		
		# check that the file is non empty
		if ( ! -s $gbk_file ) {
			MyX::Generic::File::Empty->throw(
				error => "File ($gbk_file) is empty"
			);
		}
		
		$gbk_file_of{ident $self} = $gbk_file;
		
		return 1;
	}
	
	sub set_gff_file {
		my ($self, $gff_file) = @_;
		
		# check if the parameter is defined
		if ( ! defined $gff_file ) {
			MyX::Generic::Undef::Param->throw(
				error => "Undefined parameter value (gff_file)"
			);
		}
		
		# this is an output file so I don't care if it doesn't exist
		# or if it is empty.  It will be overwritten if it does exist.
		
		$gff_file_of{ident $self} = $gff_file;
		
		return 1;
	}
	
	# Other Subroutines #
	sub _init {
		my ($self, $arg_href) = @_;
		
		$self->set_gbk_file($arg_href->{gbk_file});
		$self->set_gff_file($arg_href->{gff_file});
		
		return 1;
	}
	
	sub reformat {
		my ($self) = @_;
		
		# open the input gbk file
		my $IN = Bio::SeqIO->new(-file => $self->get_gbk_file(),
								  '-format' => 'GenBank');
		
		# open the output gff file
		open my $OUT, ">", $self->get_gff_file() or 
			MyX::Generic::File::CannotOpen->throw(
				error => "Cannot open file",
				file_name => $self->get_gff_file(),
			);
		
		# get each sequence.  here a sequence is likely either a contig or scaffold
		my $gffIO;
		while ( my $seq = $IN->next_seq() ) {
			for my $feat_obj ( $seq->get_SeqFeatures() ) {
				
				# skip things that are not CDS or tRNA
				if ( $feat_obj->primary_tag() !~ m/CDS|tRNA/i ) {
					next;
				}
				
				# remove the translation feature to reduce the size of the output
				if ( $feat_obj->has_tag("translation") ) {
					$feat_obj->remove_tag("translation");
				}
				
				$gffIO = Bio::Tools::GFF->new(-gff_version => 3);
				print $OUT $feat_obj->gff_string($gffIO), "\n";
			}
		}
		
		return 1;
	}
}

1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::FileReformat::GbkToGff - Converts Genbank to GFF


=head1 VERSION

This document describes BioUtils::FileReformat::GbkToGff version 1.0.11


=head1 SYNOPSIS

    use BioUtils::FileReformat::GbkToGff;

	my $gbk_to_gff = BioUtils::FileReformat::GbkToGff->new({
		gbk_file => "my_genome.gbk",
		gff_file => "my_genome.gff"
	});
	
	# run the reformating
	$gbk_to_gff->reformat();
	
	# reset the gkb input file and gff ouptut file
	$gbk_to_gff->set_gbk_file("my_new_genome.gbk");
	$gbk_to_gff->set_gff_file("my_new_genome.gff");
	
	# get the current setting for the input gbk file
	$gbk_to_gff->get_gbk_file();
	
	# get the current setting for the output gff file
	$gff_to_gff->get_gff_file();
	
	# A shorthand version
	BioUtils::FileReformat::GbkToGff->new({
		gbk_file => "my_genome.gbk",
		gff_file => "my_genome.gff"
	})->reformat();
  
  
=head1 DESCRIPTION

BioUtils::FileReformat::GbkToGff reformats a standard genbank file into a
GFF3 file.  Currently only CDS and tRNA tags in the genbank file are being
output in the GFF3 output file.  Also the /translation tag is not being
output to reduce the size of the output files.


=head1 DIAGNOSTICS

All errors are pass out of methods as MyX::Generic objects.  See MyX::Generic
documentaiton for infomation about catching and handling those errors.


=head1 CONFIGURATION AND ENVIRONMENT
  
BioUtils::FileReformat::GbkToGff requires no configuration files or environment variables.


=head1 DEPENDENCIES

Carp
Log::Log4perl qw(:easy)
Readonly
Class::Std::Utils
List::MoreUtils qw(any)
version
MyX::Generic
Bio::SeqIO
Bio::Tools::GFF


=head1 INCOMPATIBILITIES

None reported.


=head1 METHODS
	
=over

	new
	reformat
	get_gbk_file
	get_gff_file
	set_gbk_file
	set_gff_file
	_init

=back

=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: BioUtils::FileReformat::GbkToGff->new({
				gbk_file => "input_genbank.gbk",
				gff_file => "output.gff"
			});
	Function: Creates a new BioUtils::FileReformat::GbkToGff object
	Returns: BioUtils::FileReformat::GbkToGff
	Args: -gbk_file => a file path to the input genbank file
	      -gff_file => a file path to the output gff file
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA
	
=head2 reformat

	Title: reformat
	Usage: $GbkToGff_obj->reformat()
	Function: Reformats the Genbank file into GFF3
	Returns: 1 on success
	Args: NA
	Throws: NA
	Comments: Right now this function only prints CDS and tRNA tags.  And it
	          removes the translation tag from the GFF3 output to reduce the
			  output file size.
	See Also: NA
	
=head2 get_gbk_file

	Title: get_gbk_file
	Usage: $GbkToGff_obj->get_gbk_file()
	Function: Returns input genbank file path
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 get_gff_file

	Title: get_gff_file
	Usage: $GbkToGff_obj->get_gff_file()
	Function: Returns ouptut GFF3 file path
	Returns: str
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 set_gbk_file

	Title: set_gbk_file
	Usage: $GbkToGff_obj->set_gbk_file("file.gbk")
	Function: Sets the input genbank file path
	Returns: 1 on success
	Args: Path string to input genbank file
	Throws: MyX::Generic::Undef::Param
	        MyX::Generic::DoesNotExist::File
			MyX::Generic::File::Empty
	Comments: NA
	See Also: NA
	
=head2 set_gff_file

	Title: set_gff_file
	Usage: $GbkToGff_obj->set_gff_file("file.gff")
	Function: Sets the ouptut GFF3 file path
	Returns: 1 on success
	Args: Path string to output GFF3 file
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA
	
=head2 _init

	Title: _init
	Usage: $GbkToGff_obj->_init({
	            gbk_file => "input_genbank.gbk",
				gff_file => "output.gff"
			});
	Function: Initializes the object by setting the parameter values
	Returns: 1 on success
	Args: -gbk_file => a file path to the input genbank file
	      -gff_file => a file path to the output gff file
	Throws: NA
	Comments: Do NOT call this method.  It is private
	See Also: NA


=head1 BUGS AND LIMITATIONS

Please report any bugs or feature requests to
C<bug-<RT NAME>@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

=head2 TO DO

- make options to include other tags besides CDS and tRNA
- make an option to include the /translation tag

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

