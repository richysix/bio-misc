#!/usr/bin/env perl

# PODNAME: add_log2_fold_change.pl
# ABSTRACT: Add log2 fold change to table of data

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-11-04

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use List::Util qw( sum );

# Default options
my @numerator_fields;
my @denominator_fields;
my $has_header = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Print header line if present
if ($has_header) {
    my $line = <>;
    chomp $line;
    printf "%s\tLog2 fold change\n", $line;

    if ($debug) {
        my @fields = split /\t/xms, $line;
        warn "Numerator field(s):\n";
        foreach my $numerator_field (@numerator_fields) {
            warn q{  }, $fields[ $numerator_field - 1 ], "\n";
        }
        warn "Denominator field(s):\n";
        foreach my $denominator_field (@denominator_fields) {
            warn q{  }, $fields[ $denominator_field - 1 ], "\n";
        }
    }
}

# Iterate over STDIN and add log2 fold change
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $mean_numerator =
      sum( map { $fields[ $_ - 1 ] } @numerator_fields ) /
      scalar @numerator_fields;
    my $mean_denominator =
      sum( map { $fields[ $_ - 1 ] } @denominator_fields ) /
      scalar @denominator_fields;
    my $log2_fold_change = q{-};
    if ( $mean_numerator && $mean_denominator ) {
        $log2_fold_change = log( $mean_numerator / $mean_denominator ) / log 2;
    }
    printf "%s\t%s\n", $line, $log2_fold_change;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'numerator=i@{1,}'   => \@numerator_fields,
        'denominator=i@{1,}' => \@denominator_fields,
        'header'             => \$has_header,
        'debug'              => \$debug,
        'help'               => \$help,
        'man'                => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !@numerator_fields ) {
        pod2usage("--numerator must be specified\n");
    }
    if ( !@denominator_fields ) {
        pod2usage("--denominator must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

add_log2_fold_change.pl

Add log2 fold change to table of data

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file on STDIN and adds a new final column
containing the log2 fold change of the specified numerator (e.g. mutant) and
denominator (e.g. sibling) fields.

=head1 EXAMPLES

    perl \
        add_log2_fold_change.pl \
        --numerator 1 3 5 --denominator 2 4 6 \
        < input.tsv > output.tsv

=head1 USAGE

    add_log2_fold_change.pl
        [--numerator field...]
        [--denominator field...]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--numerator FIELD...>

The field (or fields) that specify the numerator (e.g. mutant).

=item B<--denominator FIELD...>

The field (or fields) that specify the denominator (e.g. sibling).

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
