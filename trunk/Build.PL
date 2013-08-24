use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'BioUtils',
    license             => 'perl',
    dist_author         => 'Scott Yourstone <scott.yourstone81@gmail.com>',
    dist_version_from   => 'lib/BioUtils/FastaSeq.pm',
    requires => {
        'Test::More' => 0,
        'version'    => 0,
    },
	build_requires => {
		'Carp' => 0,
		'Chart::Graph::Gnuplot' => 0,
		'Class::Std::Utils' => 0,
		'Config::Std' => 0,
		'Data::Dumper' => 0,
		'Exception::Class' => 0,
		'Getopt::Long' => 0,
		'Exporter' => 0,
		'File::Temp' => 0,
		'File::Basename' => 0,
		'List::MoreUtils' => 0,
		'Pod::Usage' => 0,
		'Readonly' => 0,
		'Scalar::Util' => 0,
		'Statistics::Descriptive' => 0,
		'Test::Exception' => 0,
		'Test::More' => 0,
		'Test::Pod' => 0,
		'Test::Pod::Coverage' => 0,
		'Test::Warn' => 0,
		'version' => 0,
	},
    add_to_cleanup      => [ 'BioUtils-*' ],
    recursive_test_files => 1,
);

$builder->create_build_script();
