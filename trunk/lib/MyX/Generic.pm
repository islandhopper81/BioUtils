package MyX::Generic;

use version; our $VERSION = qv('1.0.0');

use Exception::Class (
    'MyX::Generic' => {
    },
    
    'MyX::Generic::OutOfBounds' => {
        isa => 'MyX::Generic',
        fields => [ 'index' ],
    },
    
    'MyX::Generic::Undef' => {
        isa => 'MyX::Generic',
    },
    
    'MyX::Generic::Undef::Param' => {
        isa => 'MyX::Generic::Undef',
        fields => [ 'usage' ],
    },
    
    'MyX::Generic::Undef::Attribute' => {
        isa => 'MyX::Generic::Undef',
        fields => [ 'att_name' ],
    },
    
    'MyX::Generic::UnmatchedRegex' => {
        isa => 'MyX::Generic',
        fields => [ 'string_to_match', 'regex' ],
    },
    
    'MyX::Generic::DoesNotExist' => {
        isa => 'MyX::Generic',
    },
    
    'MyX::Generic::DoesNotExist::File' => {
        isa => 'MyX::Generic::DoesNotExist',
        fields => [ 'file_name', ],
    },
    
    'MyX::Generic::DoesNotExist::Dir' => {
        isa => 'MyX::Generic::DoesNotExist',
        fields => [ 'dir_name', ],
    },
    
    'MyX::Generic::Digit' => {
        isa => 'MyX::Generic',
    },
    
    'MyX::Generic::Digit::MustBeDigit' => {
        isa => 'MyX::Generic::Digit',
        fields => ['value'],
    },
    
    'MyX::Generic::Digit::TooSmall' => {
        isa => 'MyX::Generic::Digit',
        fields => ['value', 'MIN'],
    },
    
    'MyX::Generic::Digit::TooBig' => {
        isa => 'MyX::Generic::Digit',
        fields => ['value', 'MAX'],
    },
    
    # OOB == Out Of Bounds
    'MyX::Generic::Digit::OOB' => {
        isa => 'MyX::Generic::Digit',
        fields => ['value', 'MAX', 'MIN'],
    },
    
    'MyX::Generic::File' => {
        isa => 'MyX::Generic',
    },
    
    'MyX::Generic::File::Empty' => {
        isa => 'MyX::Generic::File',
        fields => ['file_name'],
    },
    
    'MyX::Generic::File::BadExtension' => {
        isa => 'MyX::Generic::File',
        fields => ['file_name', 'ext'],
    },
    
    'MyX::Generic::BadValue' => {
        isa => 'MyX::Generic',
        fields => ['value'],
    },
);

# Causes trace to be printed with message.  Doesn't always work!
MyX::Generic->Trace(1);

1;
__END__


#######
# POD #
#######
=head1 NAME

MyX::Generic - A hierarchy of exceptions that can be used in many
different objects

=head1 VERSION

This documentation refers to MyX::Generic version 1.0.0.

=head1 Included Modules

    Exception::Class

=head1 Inherit

NA

=head1 SYNOPSIS

    # Throw the most basic of Generic exceptions
    use MyX::Generic;
    if ( ... ) {   # Some code looking for an error
        MyX::Generic->throw(
            error => 'A generic exception'
        );
    }
    
    # In caller catch the Generic exception
    eval { ... };
    if ( my $e = MyX::Generic->caught() ) {
        # Do something to handle the exception like print an error message
        print $e->error(), " via package ", $e->package(), " at ", $e->file,
            " line ", $e->line();
    }
    

=head1 DESCRIPTION

MyX::Generic holds a hierarchy of generic exception classes.  A generic
execption is an exception that is commonly found in many different classes
reguardless of class function.  These exceptions are included in this file to
limit copy and pasting to exception files that are specific for individual
objects.  Use caution in changing the behaviour of already writen exceptions
because they are probably being used my many different classes according to
their current implementation.  Careless implementation changes could have
drastic effects on the functionality of other classes that use these generic
exceptions.

For more information what can be done when throwing and catching an exception
see Exception::Class and Exception::Class::Base.

=head1 CLASSES

=over

    MyX::Generic
    MyX::Generic::OutOfBounds
    MyX::Generic::Undef
    MyX::Generic::Undef::Param
    MyX::Generic::Undef::Attribute
    MyX::Generic::UnmatchedRegex
    MyX::Generic::DoesNotExist
    MyX::Generic::DoesNotExist::File
    MyX::Generic::DoesNotExist::Dir
    MyX::Generic::Digit
    MyX::Generic::Digit::MustBeDigit
    MyX::Generic::Digit::TooSmall
    MyX::Generic::Digit::TooBig
    MyX::Generic::Digit::OOB
    MyX::Generic::File
    MyX::Generic::File::Empty
    MyX::Generic::File::BadExtension
    MyX::Generic::BadValue
    
=back

=head1 CLASSES DESCRIPTION

=head2 MyX::Generic
    
    Title: MyX::Generic
    Throw Usage: MyX::Generic->throw(
                    error => 'Any generic error message'
                 );
    Catch Usage: if ( my $e = MyX::Generic->caught() ) { ... }
    Function: Throw/Catch a MyX::Generic exception
    Fields: error => an error message
    Inherits: NA
    Comments: NA
    See Also: NA

=head2 MyX::Generic::OutOfBounds

    Title: MyX::Generic::OutOfBounds
    Throw Usage: MyX::Generic::OutOfBounds->throw(
                    index => $i
                 );
    Catch Usage: if ( my $e = MyX::Generic::OutOfBounds->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::OutOfBounds exception when
              trying to access something beyond the bounds of a data structure
              (typically and array).
    Fields: index => The index that is out of bounds
    Inherits: MyX::Generic
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Undef

    Title: MyX::Generic::Undef
    Throw Usage: MyX::Generic::Undef->throw();
    Catch Usage: if ( my $e = MyX::Generic::Undef->caught() ) { ... }
    Function: Throw/Catch a MyX::Generic::Undef exception when trying
              to operate on something that is undefined.
    Fields: See MyX::Generic
    Inherits: MyX::Generic
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Undef::Param

    Title: MyX::Generic::Undef::Param
    Throw Usage: MyX::Generic::Undef::Param->throw(
                    usage => 'my_function($file)'
                 );
    Catch Usage: if ( my $e = MyX::Generic::Undef::Param->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::Undef::Param exception when
              a required parameter in a function is undefined.
    Fields: usage => a usage statement on how a function should correctly be
                    called
    Inherits: MyX::Generic::Undef
    Comments: NA
    See Also: NA

=head2 MyX::Generic::Undef::Attribute

    Title: MyX::Generic::Undef::Attribute
    Throw Usage: MyX::Generic::Undef::Attribute->throw();
    Catch Usage: if ( my $e = MyX::Generic::Undef::Attribute->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::Undef::Attribute exception
              when a required object attribute is undefined.
    Fields: att_name => the attribute that is undefined
    Inherits: MyX::Generic::Undef
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::UnmatchedRegex

    Title: MyX::Generic::UnmatchedRegex
    Throw Usage: MyX::Generic::UnmatchedRegex->throw();
    Catch Usage: if ( my $e = MyX::Generic::UnmatchedRegex->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::UnmatchedRegex exception
              when a regular expression fails to match in an anticipated manner.
              Usually this means there is an unexpected line format.  It could
              also mean there is bad user input.
    Fields: string => The string we are trying to match
            regex => A string showing the regular expression we are using
    Inherits: MyX::Generic
    Comments: Using the string and regex one can perform indepth testing using
              an online regex tester.
    See Also: NA
    
=head2 MyX::Generic::DoesNotExist

    Title: MyX::Generic::DoesNotExist
    Throw Usage: MyX::Generic::DoesNotExist->throw();
    Catch Usage: if ( my $e = MyX::Generic::DoesNotExist->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::DoesNotExist exception when
              trying to operate on something that does not exist.
    Fields: See MyX::Generic
    Inherits: MyX::Generic
    Comments: NA
    See Also: NA

=head2 MyX::Generic::DoesNotExist::File

    Title: MyX::Generic::DoesNotExist::File
    Throw Usage: MyX::Generic::DoesNotExist::File->throw();
    Catch Usage: if ( my $e = MyX::Generic::DoesNotExist::File->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::DoesNotExist::File exception
              when trying to operate on a file that does not exist.
    Fields: file_name => the file path that cannot be found
    Inherits: MyX::Generic::File
    Comments: NA
    See Also: NA

=head2 MyX::Generic::DoesNotExist::Dir

    Title: MyX::Generic::DoesNotExist::Dir
    Throw Usage: MyX::Generic::DoesNotExist::Dir->throw();
    Catch Usage: if ( my $e = MyX::Generic::DoesNotExist::Dir->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::DoesNotExist::Dir exception
              when trying to operate on a directory that does not exist.
    Fields: dir_name => the directory path that does not exist
    Inherits: MyX::Generic::Dir
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Digit

    Title: MyX::Generic::Digit
    Throw Usage: MyX::Generic::Digit->throw();
    Catch Usage: if ( my $e = MyX::Generic::Digit->caught() ) { ... }
    Function: Throw/Catch a MyX::Generic::Digit exception when
              opperating on digits
    Fields: See MyX::Generic
    Inherits: MyX::Generic
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Digit::MustBeDigit

    Title: MyX::Generic::DigitMustBeDigit
    Throw Usage: MyX::Generic::Digit::MustBeDigit->throw();
    Catch Usage: if ( my $e = MyX::Generic::Digit::MustBeDigit->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::Digit::MustBeDigit exception
              when opperating on digits that should be digits but are not.
    Fields: value => the value that is not a digit
    Inherits: MyX::Generic::Digit
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Digit::TooSmall

    Title: MyX::Generic::Digit::TooSmall
    Throw Usage: MyX::Generic::TooSmall->throw();
    Catch Usage: if ( my $e = MyX::Generic::Digit::TooSmall->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::Digit::TooSmall exception
              when opperating on digits are smaller than they should be.
    Fields: value => the value that is too small
            MIN => the minimum accepted value
    Inherits: MyX::Generic::Digit
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Digit::TooBig

    Title: MyX::Generic::Digit::TooBig
    Throw Usage: MyX::Generic::TooBig->throw();
    Catch Usage: if ( my $e = MyX::Generic::Digit::TooBig->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::Digit::TooBig exception when
              opperating on digits are bigger than they should be.
    Fields: value => the value that is too big
            MAX => the maximum accepted value
    Inherits: MyX::Generic::Digit
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::Digit::OOB

    Title: MyX::Generic::Digit::OOB
    Throw Usage: MyX::Generic::Digit::OOB->throw();
    Catch Usage: if ( my $e = MyX::Generic::Digit::OOB->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::Digit::OOB exception when
              opperating on digits beyond their bounds.
    Fields: value => the value that is too big
            MIN => the minimum accepted value
            MAX => the minimum accepted value
    Inherits: MyX::Generic::Digit
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::File

    Title: MyX::Generic::File
    Throw Usage: MyX::Generic::File->throw();
    Catch Usage: if ( my $e = MyX::Generic::File->caught() ) { ... }
    Function: Throw/Catch a MyX::Generic::File exceptions when
              opperating on files.
    Fields: See MyX::Generic
    Inherits: MyX::Generic
    Comments: NA
    See Also: NA

=head2 MyX::Generic::File::Empty

    Title: MyX::Generic::File::Empty
    Throw Usage: MyX::Generic::File::Empty->throw();
    Catch Usage: if ( my $e = MyX::Generic::File::Empty->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::File::Empty exceptions when
              opperating on empty files.
    Fields: file_name => the file path wich has an empty file
    Inherits: MyX::Generic::File
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::File::BadExtension

    Title: MyX::Generic::File::BadExtension
    Throw Usage: MyX::Generic::File::BadExtension->throw();
    Catch Usage: if ( my $e = MyX::Generic::File::BadExtension->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::File::BadExtension
              exceptions when opperating on files with unrecignized exceptions.
    Fields: file_name => the file with the bad extention
            ext => the extention that is bad
    Inherits: MyX::Generic::File
    Comments: NA
    See Also: NA
    
=head2 MyX::Generic::BadValue

    Title: MyX::Generic::BadValue
    Throw Usage: MyX::Generic::BadValue->throw();
    Catch Usage: if ( my $e = MyX::Generic::BadValue->caught() )
                    { ... }
    Function: Throw/Catch a MyX::Generic::BadValue exceptions when
              opperating on values that are bad for one reason or another.
    Fields: value => the value that is bad
    Inherits: MyX::Generic::File
    Comments: NA
    See Also: NA

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




