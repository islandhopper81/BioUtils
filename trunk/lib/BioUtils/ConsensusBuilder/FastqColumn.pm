package BioUtils::ConsensusBuilder::FastqColumn;

# This obect is a data structure that holds an aligned column and their quality
# values.  It is also repsonsible for calling the consensus of that column.

use strict;
use warnings;

use Class::Std::Utils;
use Carp;
use Readonly;
use version; our $VERSION = qv('1.0.3');
use BioUtils::Codec::QualityScores qw( int_to_illumina_1_8 illumina_1_8_to_int);

{
    Readonly my $NEW_USAGE => q{ new() };
    
    # Attributes of #
    my %bases_of;

    # Setters #
    sub addBase;
    
    # Getters #
    sub getConBaseAndQual;
    sub _by_probability;
    sub _base_probability;
    sub _qual_to_prob;
    sub getIupacCode;
    sub _getMeanQual;
    sub _getSortedTopBases;
    sub getBaseCount;
    sub getBaseQuals;
    
    # Others #
    sub _resolveTie;
    

    ###############
    # Constructor #
    ###############
    sub new {
        my ($class) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        my %bases = ();
        $bases{A}{count} = 0;
        $bases{T}{count} = 0;
        $bases{C}{count} = 0;
        $bases{G}{count} = 0;
        $bases{N}{count} = 0;
        $bases{dash}{count} = 0;
        
        my @A_quals = ();
        $bases{A}{quals} = \@A_quals;
        my @T_quals = ();
        $bases{T}{quals} = \@T_quals;
        my @C_quals = ();
        $bases{C}{quals} = \@C_quals;
        my @G_quals = ();
        $bases{G}{quals} = \@G_quals;
        my @N_quals = ();
        $bases{N}{quals} = \@N_quals;
        my @dash_quals = ();
        $bases{dash}{quals} = \@dash_quals;
        
        $bases_of{ident $new_obj} = \%bases;
    
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub addBase($$) {
        my ($self, $base, $qual) = @_;
        
        my $_bases_of = $bases_of{ident $self};
        
        $base =~ tr/[a-z]/[A-Z]/;
        $base =~ tr/\\./-/; # the Bioperl clustalw SimpleAlign object stores gaps as dots.
        
        # The base should be either A, T, C, G, N, or -
        if ( $base !~ m/[ATCGN-]/i ) {
            warn "$0 -- base not added to column: $base\n";
        }
        
        # The quality value is encoded.  Any ASCII enc-8 value is acceptable.
    
        # Substitute '-' characters to the word dash.  It makes a more clear key
        $base =~ s/-/dash/;
        
        $_bases_of->{$base}{count}++;
        push @{$_bases_of->{$base}{quals}}, $qual;
        
        return 1;
    }
    
    
    
    ###########
    # Getters #
    ###########
    sub getConBaseAndQual() {
        my ($self) = @_;
        
        my $_bases_of = $bases_of{ident $self};
        
        my $topBases;  # First return value
        my $meanQual;  # Second return value
        my $topCount = 0;
        
        foreach my $base ( keys %{$_bases_of} ) {
            if ( $_bases_of->{$base}{count} > $topCount ) {
                $topBases = $base;
                $topCount = $_bases_of->{$base}{count};
            }
            elsif ( $_bases_of->{$base}{count} == $topCount ) {
                $topBases .= $base;
            }
        }
        
        # substitute dash for - so length test in next if statement will work
        $topBases =~ s/dash/-/;
        
        if ( length $topBases == 1 ) {
            # if the topBases is only a dash then return empty
            if ( $topBases =~ m/-|dash/ ) {return ("", "");}
            
            $meanQual = $self->_getMeanQual($topBases);
        }
        else {
            $topBases = $self->_resolveTie($topBases);
            if ( length $topBases == 1 ) {
                $meanQual = $self->_getMeanQual($topBases);
            }
            else {
                $meanQual = $self->_getMeanQual($topBases);
                $topBases = $self->getIupacCode($topBases);
            }
        }
        
        # here topBases is a single letter, but it may be an IUPAC coded base
        return ($topBases, $meanQual);
    }
    
    sub _by_probability {
        my ($self) = @_;
        
        my $_bases_of = $bases_of{ident $self};
        
        my $max_base;
        my $max_prob = 0;
        
        foreach my $true_base ( qw(A T C G) ) {
            my $temp_prob = 0;
            
            foreach my $observed ( keys %{$_bases_of} ) {
                foreach my $qual ( @{$_bases_of->{$observed}{quals}} ) {
                    my $temp_qual = illumina_1_8_to_int($qual);
                    $temp_prob += _base_probability($true_base,
                                                    $observed,
                                                    $temp_qual);
                }
                
            }
            
            if ( $temp_prob > $max_prob ) {
                $max_prob = $temp_prob;
                $max_base = $true_base;
            }
        }
        
        return ( $max_base, $max_prob );
    }
    
    sub _base_probability {
        my ($realized, $observed, $qual) = @_;
        
        if ( $realized eq $observed ) {
            return 1 - _qual_to_prob($qual);
        }
        
        return _qual_to_prob($qual) / 3;
    }
    
    sub _qual_to_prob {
        my ($qual) = @_;
        
        return 10**( -$qual / 10 );
    }
    
    sub getIupacCode($) {
        my ($self, $topBases) = @_;
        
        my $sorted = $self->_getSortedTopBases($topBases);
        
        if ( $sorted eq "A" ) { return "A"; }
        elsif ( $sorted eq "C" ) { return "C"; }
        elsif ( $sorted eq "G" ) { return "G"; }
        elsif ( $sorted eq "T" ) { return "T"; }
        elsif ( $sorted eq "AG" ) { return "R"; }
        elsif ( $sorted eq "CT" ) { return "Y"; }
        elsif ( $sorted eq "CG" ) { return "S"; }
        elsif ( $sorted eq "AT" ) { return "W"; }
        elsif ( $sorted eq "GT" ) { return "K"; }
        elsif ( $sorted eq "AC" ) { return "M"; }
        elsif ( $sorted eq "CGT" ) { return "B"; }
        elsif ( $sorted eq "AGT" ) { return "D"; }
        elsif ( $sorted eq "ACT" ) { return "H"; }
        elsif ( $sorted eq "ACG" ) { return "V"; }
        else { return "N"; }
    }
    
    sub _getSortedTopBases($) {
        my ($self, $topBases) = @_;
    
        $topBases =~ tr/[a-z]/[A-Z]/;
        
        my @bases = split(//, $topBases);
        my @sorted = sort {uc($a) cmp uc($b)} @bases;
        
        return (join('', @sorted));
    }
    
    sub _getMeanQual($) {
        my ($self, $topBases) = @_;
        
        my $_bases_of = $bases_of{ident $self};
        
        $topBases =~ tr/[a-z]/[A-Z]/;
        
        if ( $topBases !~ m/[ATCGN-]/ ) {
            warn "$0 -- _getMeanQual() -- unrecignized base in $topBases\n";
            return -1;  # this signifies an error
        }
        
        my @topBasesArr = split(//, $topBases);
        my $count = 0;
        my $total = 0;
        foreach my $base ( @topBasesArr ) {
            $base =~ s/-/dash/;
            foreach ( @{$_bases_of->{$base}{quals}} ) {
                $count++;
                $total += illumina_1_8_to_int($_);
            }
        }
        
        # NOTE: this is NOT rounding it.  It takes the int and not the remainder.
        # ie 16.6 == 16
        return int_to_illumina_1_8(int($total / $count));
    }
    
    sub getBaseCount($) {
        my ($self, $base) = @_;
        # This is mostly used for testing purposes
        
        my $_bases_of = $bases_of{ident $self};
        
        $base =~ tr/[a-z]/[A-Z]/;
        
        if ( length $base != 1 ) {
            return 0;
        }
        if ( $base !~ m/[ATCGN-]/ ) {
            return 0;
        }
        $base =~ s/-/dash/;
            
        return $_bases_of->{$base}{count};
    }
    
    sub getBaseQuals($) {
        my ($self, $base) = @_;
        # This is mostly used for testing purposes
        
        my $_bases_of = $bases_of{ident $self};
        
        $base =~ tr/[a-z]/[A-Z]/;
        
        if ( length $base != 1 ) {
            return 0;
        }
        if ( $base !~ m/[ATCGN-]/ ) {
            return 0;
        }
        $base =~ s/-/dash/;
            
        return $_bases_of->{$base}{quals};
    
    }
    
    
    ##########
    # Others #
    ##########
    sub _resolveTie($) {
        my ($self, $topBases) = @_;
        
        # change dash back to a symbol -
        $topBases =~ s/dash/-/;
        
        my @topBasesArr = split(//, $topBases);
        my $topQual = '!';
        my $newTopBases;
        
        
        # chose the base that has the highest avg quality scores
        foreach my $base (@topBasesArr) {
            my $nextQual = $self->_getMeanQual($base);
            if ( $nextQual gt $topQual ) {
                $topQual = $nextQual;
                $newTopBases = $base;
            }
            elsif ( $nextQual eq $topQual ) {
                $newTopBases .= $base;
            }
        }
        
        $newTopBases = $self->_getSortedTopBases($newTopBases);
        
        return $newTopBases;
    }
}

1;
__END__

#######
# POD #
#######
=head1 BioUtils::ConsensusBuilder::FastqColumn

FastqColumn - A class that stores an aligned column and quality values.  It is
also repsonsible for calling the consensus of that column.

=head1 VERSION

This documentation refers to FastqColumn version 1.0.3.

=head1 Included Modules

    NA

=head1 Inherit

    NA

=head1 SYNOPSIS
    
    use BioUtils::ConsensusBuilder::FastqColumn;
    my $my_fastq_col = BioUtils::ConsensusBuilder::FastqColumn->new();
    
    $my_fastq_col->addBase($base, $qual);
    my ($topBases, $meanQual) = $my_fastq_col->getConBaseAndQual();
    $my_fastq_col->getIupacCode($topBases);

    

=head1 DESCRIPTION

FastqColumn holds the bases and corresponding quality values for an aligned
column in a multiple sequence alignment (MSA).  Once all the bases and their
quality values are added to the object, getConBaseAndQual() can be called to
get the column consensus bases and its quality value.

To determine which base should be the consensus base the algorithm first counts
the number of different bases (ie. A, T, C, G, N, etc).  If there is a clear
majority the majority bases is the consensus base and the quality value for that
base is the average of all its quality scores.

For example, if a column has bases:
A, A, A, A, T

with corresponding quality values:
40, 40, 40, 38, 10

The consensus base would be A.  The consensus base quality score would be
calculated as follows:
(40 + 40 + 40 + 38) / 4 = int(39.5) => 39.

Note that the int() operation in perl does not round.

In some cases there may be a tie between the number bases.  For example, if a
column has bases:
A, A, A, T, T, T

with corresonding quality values:
40, 40, 40, 38, 38, 38

In this case the algorithm will choose the base with the highest average quality
score (A => 40, T => 38).  So A is chosen as the consensus base as the consensus
quality value is calculated as explained above.

In very rare but possible cases (especially when the column has very few bases),
there could be a tie between bases and their corresponding average quality scores.
For example, if column has bases:
A, A, T, T

with corresponding quality values:
40, 30, 40, 30

In this case the algorithm will return the string AT with an average quality
score of all the top bases (so in this case 35).  The program using this
FastqColumn will have to decide what to do at that point.  Should it randomly
pick one of the top bases?  Should it simply pick the first top base?  Should it
use an IUPAC coding?  In the event that the calling program wants to use an
IUPAC code a method (getIupacCode) has been provided with FastqColumn.

=head1 METHODS

=over

    new
    addBase
    getConBaseAndQual
    getIupacCode
    _getMeanQual
    _getSortedTopBases
    _resolveTie
    getBaseCount
    getBaseQuals
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: BioUtils::ConsensusBuilder::FastqColumn->new();
    Function: Creates a new FastqColumn object
    Returns: FastqColumn
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 addBase

    Title: addBase
    Usage: $my_fastq_col->addBase($base, $qual);
    Function: Adds a base and its quality value
    Returns: 1 on successful completion
    Args: -base => A, T, C, G, N, -, or .
          -qual => digit representing the quality value
    Throws: NA
    Comments: NA
    See Also: NA

=head2 getConBaseAndQual

    Title: getConBaseAndQual
    Usage: $my_fastq_col->getConBaseAndQual();
    Function: Gets the consensus base and it quality score
    Returns: ($base, qual)
    Args: None
    Throws: NA
    Comments: This returns an array with two items (1. base, 2. quality score).
              See the DESCRIPTION section for detail on how a consensus bases
              and its quality score are calculated.
    See Also: DESCRIPTION section

=head2 getIupacCode

    Title: getIupacCode
    Usage: $my_fastq_con->getIupacCode($topBases);
    Function: Gets the IUPAC code for a string of bases
    Returns: String
    Args: -topBases => a string of bases
    Throws: NA
    Comments: NA
    See Also: NA

=head2 _getMeanQual

    Title: _getMeanQual
    Usage: $my_fastq_con->_getMeanQual($topBases);
    Function: Gets the mean quality score for a string of bases
    Returns: int
    Args: -topBases => a string of bases
    Throws: NA
    Comments: NA
    See Also: DESCRIPTION section

=head2 _getSortedTopBases

    Title: _getSortedTopBases
    Usage: $my_fastq_con->_getSortedTopBases($topBases);
    Function: Gets the sorted string of bases
    Returns: String
    Args: -topBases => a string of bases
    Throws: NA
    Comments: Used in getIupacCode so that not all combinations of strings had
              to be checked.
    See Also: DESCRIPTION section

=head2 _resolveTie

    Title: _resolveTie
    Usage: $my_fastq_con->_resolveTie($topBases);
    Function: Attempts to resolve a tie between bases using their quality values
    Returns: String
    Args: -topBases => a string of bases
    Throws: NA
    Comments: NA
    See Also: DESCRIPTION section

=head2 getBaseCount

    Title: getBaseCount
    Usage: $my_fastq_con->getBaseCount($base);
    Function: Gets the number of times $base is in this FastqColumn
    Returns: int
    Args: -base => one base
    Throws: NA
    Comments: Mostly used for testing purposes
    See Also: NA

=head2 getBaseQuals

    Title: getBaseQuals
    Usage: $my_fastq_con->getBaseQuals($base);
    Function: Gets the array ref of quality values for $base is in this FastqColumn
    Returns: array reference
    Args: -base => one base
    Throws: NA
    Comments: Mostly used for testing purposes
    See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.


=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2012, Scott Yourstone C<< <scott.yourstone81@gmail.com> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


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
