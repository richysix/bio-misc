#!/usr/bin/env perl

# PODNAME: keep_best_duplicate.pl
# ABSTRACT: Identify duplicates by specified fields and keep best duplicate

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-07-13

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Scalar::Util qw( looks_like_number );

# Default options
my @key_fields;
my $best_field;
my $best_type         = 'max';
my $non_numeric_value = 0;
my $has_header        = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Print header line if present
if ($has_header) {
    my $line = <>;
    chomp $line;
    printf "%s\n", $line;

    if ($debug) {
        my @fields = split /\t/xms, $line;
        warn "Key field(s):\n";
        foreach my $key_field (@key_fields) {
            warn q{  }, $fields[ $key_field - 1 ], "\n";
        }
        warn "Best field:\n";
        warn q{  }, $fields[ $best_field - 1 ], "\n";
        warn "Best type:\n";
        warn q{  }, $best_type, "\n";
        warn "Value of non-numeric best:\n";
        warn q{  }, $non_numeric_value, "\n";
    }
}

# Iterate over STDIN and eliminate duplicates
my %data_for;
my %best_for;
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $key = join "\n", map { $fields[ $_ - 1 ] } @key_fields;
    my $best = $fields[ $best_field - 1 ];
    if ( !looks_like_number($best) ) {
        $best = $non_numeric_value;
    }
    if ( !exists $data_for{$key} ) {
        $data_for{$key} = $line;
        $best_for{$key} = $best;
    }
    elsif (( $best_type eq 'min' && $best < $best_for{$key} )
        || ( $best_type eq 'max' && $best > $best_for{$key} ) )
    {
        $data_for{$key} = $line;
        $best_for{$key} = $best;
    }
}

# Output data (now ordered by key)
foreach my $key ( sort keys %data_for ) {
    printf "%s\n", $data_for{$key};
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'key=i@{1,}'          => \@key_fields,
        'best=i'              => \$best_field,
        'type=s'              => \$best_type,
        'non_numeric_value=i' => \$non_numeric_value,
        'header'              => \$has_header,
        'debug'               => \$debug,
        'help'                => \$help,
        'man'                 => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !@key_fields ) {
        pod2usage("--key must be specified\n");
    }
    if ( !$best_field ) {
        pod2usage("--best must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

keep_best_duplicate.pl

Identify duplicates by specified fields and keep best duplicate

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file on STDIN and looks for duplicates
specified by a key of one or more fields. Duplicates are removed and the best
duplicate is kept according to the value of another field.

=head1 EXAMPLES

    perl \
        keep_best_duplicate.pl \
        < input.tsv > output.tsv

    perl \
        keep_best_duplicate.pl \
        --header --key 10 --best 7 --type min \
        < input.tsv > output.tsv

=head1 USAGE

    keep_best_duplicate.pl
        [--key field...]
        [--best field]
        [--type min|max]
        [--non_numeric_value int]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--key FIELD...>

The field (or fields) that specify the key to identify duplicates.

=item B<--best FIELD>

The field that identifies the best duplicate.

=item B<--type min|max>

Whether the best duplicate is a minimum or a maximum (defaults to max).

=item B<--non_numeric_value INT>

Value assigned for non-numeric values when identifying best duplicate.

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
