#!/usr/bin/env perl

# PODNAME: merge_duplicates_by_type.pl
# ABSTRACT: Identify duplicates by specified fields and merge all other fields

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-04-13

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Sort::Naturally;

# Default options
my @key_fields;
my @sum_fields;
my @mean_fields;
my $has_header = 0;
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
        warn "Sum field(s):\n";
        foreach my $sum_field (@sum_fields) {
            warn q{  }, $fields[ $sum_field - 1 ], "\n";
        }
        warn "Mean field(s):\n";
        foreach my $mean_field (@mean_fields) {
            warn q{  }, $fields[ $mean_field - 1 ], "\n";
        }
    }
}

# Iterate over STDIN and eliminate duplicates
my %data_for;
my $num_fields;
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $key = join "\t", map { $fields[ $_ - 1 ] } @key_fields;
    push @{ $data_for{$key} }, \@fields;
    $num_fields = scalar @fields;
}

output( \%data_for, $num_fields );

sub output {
    my ( $data_for, $num_fields ) = @_;    ## no critic (ProhibitReusedNames)

    my %is_key  = map { $_ => 1 } @key_fields;
    my %is_sum  = map { $_ => 1 } @sum_fields;
    my %is_mean = map { $_ => 1 } @mean_fields;

    # Output data (now ordered by key)
    foreach my $key ( sort { ncmp( $a, $b ) } keys %{$data_for} ) {
        if ( scalar @{ $data_for->{$key} } == 1 ) {

            # If key is unique then just print line as is
            printf "%s\n", join "\t", @{ $data_for->{$key}->[0] };
            next;
        }
        my @output;
        foreach my $i ( 1 .. $num_fields ) {
            if ( $is_key{$i} ) {
                push @output, $data_for->{$key}->[0][ $i - 1 ];
                next;
            }
            my @values;
            foreach my $line ( @{ $data_for->{$key} } ) {
                if ( defined $line->[ $i - 1 ] && $line->[ $i - 1 ] ne q{} ) {
                    push @values, $line->[ $i - 1 ];
                }
            }
            my $output_value = 0;
            if ( $is_sum{$i} || $is_mean{$i} ) {
                foreach my $value (@values) {
                    $output_value += $value;
                }
                if ( $is_mean{$i} && scalar @values ) {
                    $output_value /= scalar @values;
                }
                elsif ( !scalar @values ) {
                    $output_value = q{};
                }
            }
            else {
                my $all_same = 1;
                foreach my $j ( 1 .. scalar @values - 1 ) {
                    if ( $values[0] ne $values[$j] ) {
                        $all_same = 0;
                        last;
                    }
                }
                if ($all_same) {
                    $output_value = $values[0];
                }
                else {
                    $output_value = join q{,}, @values;
                }
            }
            push @output, $output_value;
        }
        printf "%s\n", join "\t", @output;
    }

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'key=i@{1,}'  => \@key_fields,
        'sum=i@{1,}'  => \@sum_fields,
        'mean=i@{1,}' => \@mean_fields,
        'header'      => \$has_header,
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

    if ( !@key_fields ) {
        pod2usage("--key must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

merge_duplicates_by_type.pl

Identify duplicates by specified fields and merge all other fields

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file on STDIN and looks for duplicates
specified by a key of one or more fields. Duplicates are removed and all other
fields are merged. By default fields are merged by comma separation (collapsed
to a single value if all the same), but fields to be merged by sum or mean can
be specified.

=head1 EXAMPLES

    perl \
        merge_duplicates_by_type.pl \
        < input.tsv > output.tsv

    perl \
        merge_duplicates_by_type.pl \
        --header --key 1 --sum 2 3 4 --mean 6 \
        < input.tsv > output.tsv

=head1 USAGE

    merge_duplicates_by_type.pl
        [--key field...]
        [--sum field...]
        [--mean field...]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--key FIELD...>

The field (or fields) that specify the key to identify duplicates.

=item B<--sum FIELD...>

The field (or fields) to be merged by summing the values.

=item B<--mean FIELD...>

The field (or fields) to be merged by taking the mean of the values.

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
