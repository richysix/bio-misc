#!/usr/bin/env perl

# PODNAME: add_sum.pl
# ABSTRACT: Add sum to table of data

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-11-12

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use List::Util qw( sum );

# Default options
my @sum_fields;
my $has_header = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Print header line if present
if ($has_header) {
    my $line = <>;
    chomp $line;
    printf "%s\tSum\n", $line;

    if ($debug) {
        my @fields = split /\t/xms, $line;
        warn "Sum field(s):\n";
        foreach my $sum_field (@sum_fields) {
            warn q{  }, $fields[ $sum_field - 1 ], "\n";
        }
    }
}

# Iterate over STDIN and add sum
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $sum = sum( map { $fields[ $_ - 1 ] } @sum_fields );
    printf "%s\t%s\n", $line, $sum;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'sum=i@{1,}' => \@sum_fields,
        'header'     => \$has_header,
        'debug'      => \$debug,
        'help'       => \$help,
        'man'        => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !@sum_fields ) {
        pod2usage("--sum must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

add_sum.pl

Add sum to table of data

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file on STDIN and adds a new final column
containing the sum of the specified fields.

=head1 EXAMPLES

    perl \
        add_sum.pl \
        --sum 1 2 3 4 5 6  \
        < input.tsv > output.tsv

=head1 USAGE

    add_sum.pl
        [--sum field...]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--sum FIELD...>

The field (or fields) to be summed.

=item B<--header>

Input has a header line.

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

This software is Copyright (c) 2014 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
