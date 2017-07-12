#!/usr/bin/env perl

# PODNAME: upset_rnaseq.pl
# ABSTRACT: Convert RNA-Seq output into UpSet format

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-04-24

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my @input_files;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

my @sets = remove_common_prefix_suffix(@input_files);

my %upset;

my $i = 0;
foreach my $file (@input_files) {
    my $set = $sets[$i];    ## no critic (ProhibitAmbiguousNames)
    open my $fh, '<', $file;
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        chomp $line;
        my @fields = split /\t/xms, $line;
        my $id = $fields[0];
        $upset{$id}{$set} = 1;
    }
    close $fh;
    $i++;
}

printf "%s\r\n", join q{,}, 'Gene', @sets;
foreach my $id ( sort keys %upset ) {
    printf "%s\r\n", join q{,}, $id, map { $upset{$id}{$_} || 0 } @sets;
}

sub remove_common_prefix_suffix {
    my @strings = @_;

    foreach ( 0 .. 1 ) {
        @strings = map { scalar reverse } @strings;
        my $common_prefix_idx = 0;
      CHAR: while (1) {
            my $char = substr $strings[0], $common_prefix_idx, 1;
            foreach my $string (@strings) {
                last CHAR if substr( $string, $common_prefix_idx, 1 ) ne $char;
            }
            $common_prefix_idx++;
        }
        if ($common_prefix_idx) {
            @strings = map { substr $_, $common_prefix_idx } @strings;
        }
    }

    return @strings;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_files=s@{2,}' => \@input_files,
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

    # Check options
    if ( !@input_files ) {
        pod2usage("--input_files must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

upset_rnaseq.pl

Convert RNA-Seq output into UpSet format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes tab-delimited files of RNA-Seq output produced by
https://gist.github.com/iansealy/2dca28d07c0764e014df or
https://gist.github.com/iansealy/b9cbc56bd1affe10d37a and converts to the format
required by UpSet.

=head1 EXAMPLES

    perl \
        upset_rnaseq.pl \
        --input_file sig1.tsv sig2.tsv sig3.tsv > upset.csv

=head1 USAGE

    upset_rnaseq.pl
        [--input_files file...]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_files FILE...>

RNA-Seq output files (e.g. sig.tsv).

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

This software is Copyright (c) 2017 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
