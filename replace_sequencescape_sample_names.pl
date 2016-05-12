#!/usr/bin/env perl

# PODNAME: replace_sequencescape_sample_names.pl
# ABSTRACT: Replace Sequencescape samples name with readable sample names

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-05-11

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;
use DBI;

# Constants
Readonly our $HOST => 'seqw-db';
Readonly our $PORT => '3379';
Readonly our $USER => 'warehouse_ro';
Readonly our $PASS => q{};
Readonly our $NAME => 'sequencescape_warehouse';

# Default options
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connect to database
my $dsn = "dbi:mysql:database=$NAME;host=$HOST;port=$PORT";
my $dbh = DBI->connect( $dsn, $USER, $PASS );

# Cache sample names
my %name_for;
my $ary_ref = $dbh->selectall_arrayref(
    <<"SQL"
    SELECT name, supplier_name
    FROM   current_samples
    WHERE  (name LIKE '%STDY%' OR supplier_name LIKE '%STDY%')
    AND    supplier_name IS NOT NULL
SQL
);
foreach ( @{$ary_ref} ) {
    my ( $name, $supplier_name ) = @{$_};
    $name_for{$name} = $supplier_name;
}

# Iterate over STDIN
while ( my $line = <> ) {
    $line =~ s/\b(\d+STDY\d+)\b/$name_for{$1} || $1/xmsge;
    printf '%s', $line;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'debug' => \$debug,
        'help'  => \$help,
        'man'   => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

replace_sequencescape_sample_names.pl

Replace Sequencescape samples name with readable sample names

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a file and replaces Sequencescape sample names (e.g.
3407STDY6401339) with more readable sample names (e.g. zmp_ph250_E4_p2).

=head1 EXAMPLES

    perl \
        replace_sequencescape_sample_names.pl \
        < input.txt > output.txt

=head1 USAGE

    replace_sequencescape_sample_names.pl
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

None

=head1 AUTHOR

=over 4

=item *

Ian Sealy <ian.sealy@sanger.ac.uk>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2016 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
