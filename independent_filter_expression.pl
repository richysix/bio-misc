#!/usr/bin/env perl

# PODNAME: independent_filter_expression.pl
# ABSTRACT: Independently filter gene expression data based on mean expression

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-05-12

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use List::Util qw(sum);
use List::MoreUtils qw(indexes);

# Default options
my $threshold;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Extract header and identify normalised count columns
my $header = <>;
chomp $header;
my @headings = split /\t/xms, $header;
my @count_columns = indexes { m/\s normalised \s count \z/xms } @headings;
if ( !@count_columns ) {

    # Find columns after gene description
    my $i   = 0;
    my $add = 0;
    foreach my $heading (@headings) {
        if ($add) {
            push @count_columns, $i;
        }
        if ( $heading eq 'Description' || $heading eq q{"Gene description"} ) {
            $add = 1;
        }
        $i++;
    }
}
printf "%s\n", $header;
printf {*STDERR} "Filtering based on mean of %d columns:\n%s\n",
  scalar @count_columns, join "\n", @headings[@count_columns];

# Iterate over STDIN
my $filtered_count = 0;
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;

    # If second field blank then assume additional header
    if ( $fields[1] eq q{} ) {
        printf "%s\n", $line;
        next;
    }

    my $mean = sum( @fields[@count_columns] ) / scalar @count_columns;
    if ( $mean >= $threshold ) {
        printf "%s\n", $line;
    }
    else {
        $filtered_count++;
    }
}
printf {*STDERR} "Filtered %d lines\n", $filtered_count;

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'threshold=i' => \$threshold,
        'debug'       => \$debug,
        'help'        => \$help,
        'man'         => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$threshold ) {
        pod2usage("--threshold must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

independent_filter_expression.pl

Independently filter gene expression data based on mean expression

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a file of tab-delimited expression data and removes
genes/regions below a threshold mean expression.

=head1 EXAMPLES

    perl \
        independent_filter_expression.pl \
        --threshold 100 \
        < input.txt > output.txt

=head1 USAGE

    replace_sequencescape_sample_names.pl
        [--threshold INT]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--threshold INT>

Mean expression level to filter out below.

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
