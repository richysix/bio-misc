#!/usr/bin/env perl

# PODNAME: analyse_intronic_expression2.pl
# ABSTRACT: Analyse intronic expression (part 2) - decision tree

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2018-03-22

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;

# Constants
Readonly our $COUNT_THRESHOLD => 10;
Readonly our $FPKM_THRESHOLD  => 1;

# Default options
my $fpkm;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

my $threshold = $fpkm ? $FPKM_THRESHOLD : $COUNT_THRESHOLD;
my $name      = $fpkm ? 'FPKM'          : 'Count';

# Iterate over STDIN
my $header = <>;
chomp $header;
printf "%s\n", join "\t", $header,
  "Intron Segment Enclosed Total $name > $threshold",
  "Prev Exon Overlap $name & Next Exon Overlap $name > $threshold",
  "Repeat Overlap Total $name > $threshold",
  'Category';
while ( my $line = <> ) {
    chomp $line;
    if ( $line =~ m/\A [#]/xms ) {
        printf "%s\n", $line;
        next;
    }
    my @fields = split /\t/xms, $line;
    ## no critic (ProhibitMagicNumbers)
    my $counts_or_fpkm_in_intron =
        $fields[ 31 + ( $fpkm ? 1 : 0 ) ] eq q{-}      ? q{-}
      : $fields[ 31 + ( $fpkm ? 1 : 0 ) ] > $threshold ? q{y}
      :                                                  q{n};
    my $counts_or_fpkm_in_exons = $fields[ 15 + ( $fpkm ? 1 : 0 ) ] > $threshold
      && $fields[ 20 + ( $fpkm ? 1 : 0 ) ] > $COUNT_THRESHOLD ? q{y} : q{n};
    my $counts_or_fpkm_in_repeats =
        $fields[ 26 + ( $fpkm ? 1 : 0 ) ] eq q{-}      ? q{-}
      : $fields[ 26 + ( $fpkm ? 1 : 0 ) ] > $threshold ? q{y}
      :                                                  q{n};
    ## use critic
    my $category;
    if ( $counts_or_fpkm_in_intron eq q{-} ) {
        $category = 'intron-entirely-repeat';
    }
    elsif ($counts_or_fpkm_in_intron eq q{y}
        && $counts_or_fpkm_in_exons eq q{y} )
    {
        $category = 'novel-exon-or-intron-retention';
    }
    elsif ($counts_or_fpkm_in_intron eq q{y}
        && $counts_or_fpkm_in_exons eq q{n} )
    {
        $category =
          $counts_or_fpkm_in_repeats eq q{y}
          ? 'repeats-expressed-independently'
          : 'novel-exon';
    }
    elsif ($counts_or_fpkm_in_intron eq q{n}
        && $counts_or_fpkm_in_exons eq q{y} )
    {
        $category =
          $counts_or_fpkm_in_repeats eq q{y}
          ? 'repeats-expressed-independently'
          : 'only-exons-expressed';
    }
    elsif ($counts_or_fpkm_in_intron eq q{n}
        && $counts_or_fpkm_in_exons eq q{n} )
    {
        $category =
          $counts_or_fpkm_in_repeats eq q{y}
          ? 'repeats-expressed-independently'
          : 'transcriptionally-silent';
    }
    printf "%s\n", join "\t", $line, $counts_or_fpkm_in_intron,
      $counts_or_fpkm_in_exons,
      $counts_or_fpkm_in_repeats, $category;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'fpkm'  => \$fpkm,
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

analyse_intronic_expression2.pl

Analyse intronic expression (part 2) - decision tree

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script analyses the output from the previous intronic expression analysis
script to incorporate decision tree information.

=head1 EXAMPLES

    perl analyse_intronic_expression2.pl \
        < input.txt > output.txt

=head1 USAGE

    analyse_intronic_expression2.pl
        [--fpkm]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--fpkm>

Use FPKM threshold, rather than count threshold.

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
