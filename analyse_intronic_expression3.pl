#!/usr/bin/env perl

# PODNAME: analyse_intronic_expression3.pl
# ABSTRACT: Analyse intronic expression (part 3) - repeat distribution

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2018-05-21

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;

# Default options
my $file;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over file and count introns and repeats
my %gene_intron_count;
my %transcript_intron_count;
my %gene_repeat_count;
my %transcript_repeat_count;
open my $fh, q{<}, $file;    ## no critic (RequireBriefOpen)
my $header = <$fh>;
while ( my $line = <$fh> ) {
    next if $line =~ m/\A [#]/xms;
    my @fields = split /\t/xms, $line;
    my $gene = $fields[0];
    ## no critic (ProhibitMagicNumbers)
    my $transcript = $fields[3];
    my $repeats    = $fields[23];
    ## use critic
    $gene_intron_count{$gene}++;
    $transcript_intron_count{$transcript}++;
    next if $repeats eq q{-};
    my $count = scalar split /,/xms, $repeats;
    $gene_repeat_count{$gene}             += $count;
    $transcript_repeat_count{$transcript} += $count;
}
close $fh;

# Iterate over file and add counts
open $fh, q{<}, $file;    ## no critic (RequireBriefOpen)
$header = <$fh>;
chomp $header;
printf "%s\n", join "\t", $header,
  'Introns In Transcript',
  'Introns In Gene',
  'Repeats In Intron',
  'Repeats In Transcript',
  'Repeats In Gene',
  'Repeats In Intron / Transcript',
  'Repeats In Intron / Gene';
while ( my $line = <$fh> ) {
    chomp $line;
    if ( $line =~ m/\A [#]/xms ) {
        printf "%s\n", $line;
        next;
    }
    my @fields = split /\t/xms, $line;
    my $gene = $fields[0];
    ## no critic (ProhibitMagicNumbers)
    my $transcript = $fields[3];
    my $repeats    = $fields[23];
    ## use critic
    my $count = scalar split /,/xms, $repeats;
    if ( $repeats eq q{-} ) {
        $count = 0;
    }
    my $transcript_repeat_count = $transcript_repeat_count{$transcript} || 0;
    my $gene_repeat_count       = $gene_repeat_count{$gene}             || 0;
    my $transcript_ratio        = q{-};
    if ($transcript_repeat_count) {
        $transcript_ratio = $count / $transcript_repeat_count;
    }
    my $gene_ratio = q{-};
    if ($gene_repeat_count) {
        $gene_ratio = $count / $gene_repeat_count;
    }
    printf "%s\n", join "\t", $line, $transcript_intron_count{$transcript},
      $gene_intron_count{$gene}, $count, $transcript_repeat_count,
      $gene_repeat_count, $transcript_ratio, $gene_ratio;
}
close $fh;

# Get and check command line options
get_and_check_options();

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'file=s' => \$file,
        'debug'  => \$debug,
        'help'   => \$help,
        'man'    => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$file ) {
        pod2usage("--file must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

analyse_intronic_expression3.pl

Analyse intronic expression (part 3) - repeat distribution

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script analyses the output from either of the previous intronic expression
analysis scripts to incorporate repeat distribution information.

=head1 EXAMPLES

    perl analyse_intronic_expression3.pl \
        --file input.txt > output.txt

=head1 USAGE

    analyse_intronic_expression3.pl
        [--file file]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--file FILE>

File to analyse.

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
