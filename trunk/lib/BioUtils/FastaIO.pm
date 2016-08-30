package BioUtils::FastaIO;

use warnings;
use strict;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use Scalar::Util qw(openhandle);
use MyX::Generic 0.0.1;
use BioUtils::FastaSeq 1.2.0;

use version; our $VERSION = qv('1.2.0');


{
    # Global variables
    Readonly my $NEW_USAGE => q{ new( {file => , stream_type => < | >, } ) };
    
    # Class Attributes
    my %stream_type_of;
    my %file_of;
    my %fh_of;
    
    # Public Class Methods
    sub get_next_seq;
    sub write_seq;
    
    # Private Class Methods
    sub _open_fh;
	sub _close_fh;
    sub DESTROY;
    
    # Constructor
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Check that the required parameters are defined (only seq and quals_str)
        if ( any { !defined $_ } $arg_href->{file}, $arg_href->{stream_type} ) {
            MyX::Generic::Undef::Param->throw(
                                                   error => "Undefined parameter",
                                                   usage => $NEW_USAGE,
                                                   );
        }
        
        # Initialize the objects attributes
        $stream_type_of{ident $new_obj} = $arg_href->{stream_type};
        $file_of{ident $new_obj} = $arg_href->{file};
        
        # Open and store the file handle
        $fh_of{ident $new_obj} = $new_obj->_open_fh($arg_href->{file});
    
        return $new_obj;
    }

    sub get_next_seq {
        my ($self) = @_;
        
        my $fh = $fh_of{ident $self};  # stored to avoid repeated calls to ident
        
        # Check for some error in the file handle
        # ... do I really want to do these each time I call get_next_seq??
        if ( eof $fh == 1 ) {
            return 0;  # Indicates the end of the file
        }
        
        # Store the header
        my $header;
        {
            $header = readline $fh;
            chomp $header;
            
            # remove '>' symbol from header if needed
            $header =~ tr/>//d;
        }
        
        # Store the sequence up until the next header line.
        my $seq;
        my $pos = tell $fh;
        LINE: while ( my $line = readline $fh ) {
            chomp $line;
            
            # check for the next header
            if ( $line =~ m/^>/ ) {
                last;
            }
            
            $seq .= $line;  # save the growing sequence
            $pos = tell $fh;  # save the fh position
        }
        seek $fh, $pos, 0;  # move the fh pointer back before the header
        
        # Create a new FastqSeq
        my $fasta_seq = undef;
        $fasta_seq = BioUtils::FastaSeq->new({header => $header, seq => $seq});
       
        
        return $fasta_seq;
    }
    
    sub write_seq {
        my ($self, $fastq_seq) = @_;
        
        # check if the object is a fastaSeq
        
        # check if the values are all defined
        
        # Get the header
        my $header = undef;
        eval {
            $header = $fastq_seq->get_header();
        };
        if ( my $e = MyX::Generic::Undef::Attribute->caught() ) {
            MyX::Generic::Undef::Attribute->throw(
                error => 'Printing a FastqSeq without a header'
            );
        }
        
        # Get the base sequence string
        my $seq = $fastq_seq->get_seq();
        
        # Print the gathered information
        my $fh = $fh_of{ident $self};
        print $fh   qw{>}, $header, "\n",
                                    $seq, "\n",
                                    ;
        
        return 1;
    }
    
    sub _open_fh {
        my ($self, $file) = @_;
        
        open my $fh, $stream_type_of{ident $self}, $file or
            croak("Cannot open file: $file\nERROR: $!\n");

        return $fh;
    }
	
	sub _close_fh {
		my ($self) = @_;
		
		my $fh = $fh_of{ident $self};
		if ( defined $fh and openhandle($fh) ) {
			close($fh_of{ident $self});
		}
		
		return 1;
	}
    
    sub DESTROY {
        my ($self) = @_;
        
		$self->_close_fh();
        delete $stream_type_of{ident $self};
        delete $file_of{ident $self};
        delete $fh_of{ident $self};
    }
}


1; # Magic true value required at end of module
__END__

=head1 NAME

BioUtils::FastaIO - An object for reading and writing Fasta files and
BioUtils::FastaSeq objects


=head1 VERSION

This document describes BioUtils::FastaIO version 1.2.0


=head1 SYNOPSIS

    use BioUtils::FastaIO;

    my $fasta_in = BioUtils::FastaIO->new({stream_type => '<', file => $filename});
    
    # Read in all the sequences
    while ( my $seq = $fasta_in->get_next_seq() ) {
        # Do something with $seq
    }
    
    
    # Build a parser for writing FASTQ files
    my $fasta_out = BioUtils::FastaIO->new( {stream_type => '>', file => $file} );
    
    # Write a FastqSeq object to $out
    my $fasta_seq = BioUtils::FastaSeq->new( {
		header => 'Seq1',
        seq => 'ATCG'
	});
    $fasta_out->write_seq($fastq_seq);
  
  
=head1 DESCRIPTION

BioUtils::FastaIO is an object for reading and writing Fasta files and
BioUtils::FastaSeq objects.  Fasta sequences have a header which is preceeded by
the '>' character.  The next line (and in some cases multiple lines) is the
sequence string.  BioUtils::FastaIO has can handle Fasta files that have
sequence strings that span multiple lines.


=head1 METHODS 

	get_next_seq
	write_seq
	_open_fh
	_close_fh
	DESTROY


=head1 DIAGNOSTICS

=over

=item C<< Printing a FastqSeq without a header >>

In BioUtils::FastaSeq object the header attribute is not required.  However it
is requiredfor printing a BioUtils::FastaSeq object using a BioUtils::FastaIO
object.

=item C<< Undefined parameter >>

When a method is called that is missing a required parameter this error is
thown.

=back


=head1 CONFIGURATION AND ENVIRONMENT

BioUtils::FastaIO requires no configuration files or environment variables.


=head1 DEPENDENCIES

	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(croak)
	Scalar::Util
	MyX::Generic
	BioUtils::FastaSeq


=head1 INCOMPATIBILITIES

	None reported.


=head1 METHODS DESCRIPTION

=head2 new

	Title: new
	Usage: my $fasta_in = BioUtils::FastaIO->new({
				stream_type => '>', file => $filename
			});
	Function: Returns a BioUtils::FastaIO object
	Returns: BioUtils::FastaIO
	Args: -stream_type => '>' for output '<' for input
          -file => a string representing the file path
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA

=head2 get_next_seq

	Title: get_next_seq
	Usage: my $fasta_in->get_next_seq();
	Function: Returns a BioUtils::FastaSeq object or 0 when finished
	Returns: BioUtils::FastaIO or 0
	Args: NA
	Throws: MyX::Generic::Undef::Param
	Comments: This function should most commonly be encased in a while loop
	See Also: NA

=head2 write_seq

	Title: write_seq
	Usage: my $fasta_in->write_seq($seq);
	Function: Write the BioUtils::FastaSeq to a file
	Returns: 1 on successful completion
	Args: -seq => a BioUtils::FastaSeq object
	Throws: MyX::Generic::Undef::Attribute
	Comments: NA
	See Also: NA
	
=head2 _open_fh

	Title: _open_fh
	Usage: my $fasta_in->_open_fh($file);
	Function: Opens a file handle using the given file
	Returns: the file handle variable
	Args: -file => path to a file
	Throws: NA
	Comments: When a FastaIO object is create this method is automatically
			  called and the file handle is stored as an attribute to this
			  object.
	See Also: NA
	
=head2 _close_fh

	Title: _close_fh
	Usage: my $fasta_in->_close_fh();
	Function: Closes the file handle stored in this object
	Returns: 1 on successful completion
	Args: NA
	Throws: NA
	Comments: This method is called in DESTROY.  In rare cases when you open a
			  file for writing and later use a differrent BioUtils::FastaIO
			  object to read that file, it will be necessary to close the file
			  handle of the original object if it is still in scope. 
	See Also: NA

=head2 DESTROY
    
    Title: DESTROY
    Usage: $my_fasta_seq->DESTROY();
    Function: Deletes attribute hashes
    Returns: NA
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to
C<bug-fastaio@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.


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
