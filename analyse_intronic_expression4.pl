#!/usr/bin/env perl

# PODNAME: analyse_intronic_expression4.pl
# ABSTRACT: Analyse intronic expression (part 4) - normalised counts

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2018-07-03

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my $size_factor;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over STDIN
my $header = <>;
chomp $header;
printf "%s\n", join "\t", $header,
  'Intron Enclosed Normalised Count',
  'Repeat Overlap Normalised Total Count',
  'Intron Segment Enclosed Normalised Total Count';
while ( my $line = <> ) {
    chomp $line;
    if ( $line =~ m/\A [#]/xms ) {
        printf "%s\n", $line;
        next;
    }
    my @fields = split /\t/xms, $line;
    ## no critic (ProhibitMagicNumbers)
    my $intron_normalised_count =
      $fields[10] eq q{-} ? q{-} : $fields[10] / $size_factor;
    my $repeat_normalised_count =
      $fields[26] eq q{-} ? q{-} : $fields[26] / $size_factor;
    my $intron_segment_normalised_count =
      $fields[31] eq q{-} ? q{-} : $fields[31] / $size_factor;
    ## use critic
    printf "%s\n", join "\t", $line, $intron_normalised_count,
      $repeat_normalised_count, $intron_segment_normalised_count;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'size_factor=f' => \$size_factor,
        'debug'         => \$debug,
        'help'          => \$help,
        'man'           => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$size_factor ) {
        pod2usage("--size_factor must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

analyse_intronic_expression4.pl

Analyse intronic expression (part 4) - normalised counts

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script analyses the output from either of the previous intronic expression
analysis scripts to incorporate normalised read counts.

=head1 EXAMPLES

    perl analyse_intronic_expression4.pl \
        --size_factor 0.9 \
        < input.txt > output.txt

=head1 USAGE

    analyse_intronic_expression4.pl
        [--size_factor float]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--size_factor FLOAT>

Size factor (from DESeq2).

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

This software is Copyright (c) 2018 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
