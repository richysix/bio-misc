#!/usr/bin/env perl

# PODNAME: merge_duplicates.pl
# ABSTRACT: Identify duplicates by specified fields and merge other fields

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-07-29

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my @key_fields;
my @merge_fields;
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
        warn "Merge field(s):\n";
        foreach my $merge_field (@merge_fields) {
            warn q{  }, $fields[ $merge_field - 1 ], "\n";
        }
    }
}

# Iterate over STDIN and eliminate duplicates
my %data_for;
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $key = join "\t", map { $fields[ $_ - 1 ] } @key_fields;
    foreach my $merge_field (@merge_fields) {
        push @{ $data_for{$key}{$merge_field} }, $fields[ $merge_field - 1 ];
    }
}

# Output data (now ordered by key)
foreach my $key ( sort keys %data_for ) {
    ## no critic (RequireCheckedSyscalls)
    print $key;
    foreach my $merge_field (@merge_fields) {
        print "\t", ( join q{,}, @{ $data_for{$key}{$merge_field} } );
    }
    print "\n";
    ## use critic
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'key=i@{1,}'   => \@key_fields,
        'merge=i@{1,}' => \@merge_fields,
        'header'       => \$has_header,
        'debug'        => \$debug,
        'help'         => \$help,
        'man'          => \$man,
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
    if ( !@merge_fields ) {
        pod2usage("--merge must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

merge_duplicates.pl

Identify duplicates by specified fields and merge other fields

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file on STDIN and looks for duplicates
specified by a key of one or more fields. Duplicates are removed and other
specified fields are merged. Data not specified by either the key or merge
options will not be output. The order of the fields will not be preserved.

=head1 EXAMPLES

    perl \
        merge_duplicates.pl \
        < input.tsv > output.tsv

    perl \
        merge_duplicates.pl \
        --header --key 1 2 3 --merge 4 \
        < input.tsv > output.tsv

=head1 USAGE

    merge_duplicates.pl
        [--key field...]
        [--merge field...]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--key FIELD...>

The field (or fields) that specify the key to identify duplicates.

=item B<--merge FIELD...>

The field (or fields) to be merged.

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
