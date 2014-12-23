#!/usr/bin/env perl

# PODNAME: make_gap_map_baits.pl
# ABSTRACT: Make GAPmap baits from reference genome

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-12-22

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;
use List::Util qw( sum );
use Bio::SeqIO;

# Constants
Readonly my $BAIT_LENGTH => 120;
Readonly my $MAX_BAITS   => 57_750;    # Maximum number of unique 120mer baits

# Zv9 GAPmap used 52,834 120mer baits

# Default options
my $ref_fasta;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Gather statistics about reference genome
my $scaffold_count          = 0;
my $total_scaffold_length   = 0;
my $total_unmasked_length   = 0;
my $total_masked_length     = 0;
my $total_unmasked_gc       = 0;
my $total_masked_gc         = 0;
my $contig_count            = 0;
my $total_contig_length     = 0;
my $designable_contig_count = 0;
my @contig_lengths;

my $in = Bio::SeqIO->new(
    -file     => $ref_fasta,
    -format   => 'fasta',
    -alphabet => 'dna',
);
while ( my $scaffold = $in->next_seq ) {

    # Scaffold-based stats

    $scaffold_count++;
    $total_scaffold_length += $scaffold->length;
    $total_unmasked_length += ( $scaffold->seq =~ tr/AGCT/AGCT/ );
    $total_masked_length   += ( $scaffold->seq =~ tr/agct/agct/ );
    $total_unmasked_gc     += ( $scaffold->seq =~ tr/GC/GC/ );
    $total_masked_gc       += ( $scaffold->seq =~ tr/gc/gc/ );

    # Contig-based stats (assuming contigs are separated by 2+ Ns)

    my @contigs = split /NN+/xmsi, $scaffold->seq;
    $total_contig_length += sum( map { length $_ } @contigs );
    $contig_count += scalar @contigs;

    # Ignore contigs where a 120mer bait can't be designed
    @contigs = grep { m/[AGCT]{$BAIT_LENGTH}/xms } @contigs;

    $designable_contig_count += scalar @contigs;
    push @contig_lengths, map { length $_ } @contigs;
}

my $ideal_bait_spacing = int $total_contig_length / $MAX_BAITS;

printf {*STDERR} "Scaffold count:\t%d\n",        $scaffold_count;
printf {*STDERR} "Total scaffold length:\t%d\n", $total_scaffold_length;
printf {*STDERR} "Total contig length:\t%d\n",   $total_contig_length;
printf {*STDERR} "Total unmasked length:\t%d\n", $total_unmasked_length;
## no critic (ProhibitMagicNumbers)
printf {*STDERR} "GC:\t%.1f%%\n",
  ( $total_unmasked_gc + $total_masked_gc ) / $total_contig_length * 100;
printf {*STDERR} "Unmasked GC:\t%.1f%%\n",
  $total_unmasked_gc / $total_unmasked_length * 100;
## use critic
printf {*STDERR} "Contig count:\t%d\n", $contig_count;
printf {*STDERR}
"Designable contig count:\t%d (excluding where %dmer definitely can't be designed)\n",
  $designable_contig_count, $BAIT_LENGTH;
printf {*STDERR} "=> Ideal bait spacing: %d bp (%d baits)\n",
  $ideal_bait_spacing, $MAX_BAITS;

my $gc_threshold = $total_unmasked_gc / $total_unmasked_length;

my $bait_spacing = $ideal_bait_spacing;

# Iterate to find optimal bait spacing
# (i.e. where bait spacing in contigs with 2+ baits matches bait spacing in shorter contigs)
while (1) {
    my $contig_count_one_bait =
      grep { $_ <= $bait_spacing * 2 } @contig_lengths;
    my $total_contig_length_multiple_baits =
      sum( grep { $_ > $bait_spacing * 2 } @contig_lengths );

    # Calculate bait spacing in longer contigs
    my $actual_bait_spacing = $total_contig_length_multiple_baits /
      ( $MAX_BAITS - $contig_count_one_bait );

    if ( $actual_bait_spacing <= $bait_spacing ) {
        last;
    }
    else {
        $bait_spacing += 100;    ## no critic (ProhibitMagicNumbers)
    }
}

my $contig_count_one_bait =
  grep { $_ <= $bait_spacing * 2 } @contig_lengths;
my $contig_count_multiple_baits =
  grep { $_ > $bait_spacing * 2 } @contig_lengths;
my $total_contig_length_one_bait =
  sum( grep { $_ <= $bait_spacing * 2 } @contig_lengths );
my $total_contig_length_multiple_baits =
  sum( grep { $_ > $bait_spacing * 2 } @contig_lengths );

printf {*STDERR} "Contig count with 1 bait (i.e. <= %d):\t%d\n",
  $bait_spacing * 2, $contig_count_one_bait;
printf {*STDERR} "Total contig length with 1 bait:\t%d\n",
  $total_contig_length_one_bait;
printf {*STDERR} "Contig count with 2+ baits (i.e. > %d):\t%d\n",
  $bait_spacing * 2, $contig_count_multiple_baits;
printf {*STDERR} "Total contig length with 2+ baits:\t%d\n",
  $total_contig_length_multiple_baits;
printf {*STDERR} "=> Actual bait spacing: %d bp (%d baits)\n", $bait_spacing,
  $MAX_BAITS;

# Design baits
my $undesignable_contig_count;
my $undesignable_scaffold_count;
my $designed_bait_count;
$in = Bio::SeqIO->new(
    -file     => $ref_fasta,
    -format   => 'fasta',
    -alphabet => 'dna',
);
while ( my $scaffold = $in->next_seq ) {

    my @scaffold_baits;

    # Get position of each contig in genomic coordinates
    my $scaffold_seq = $scaffold->seq;
    my @coords       = (1);              # Start of first contig
    while ( $scaffold_seq =~ m/(NN+)/xmsgi ) {
        push @coords, pos($scaffold_seq) - length $1;    # End of contig
        push @coords, pos($scaffold_seq) + 1;            # Start of next contig
    }
    push @coords, $scaffold->length;                     # End of last contig

    # Get each contig
    while (@coords) {
        my $start = shift @coords;
        my $end   = shift @coords;
        next if $start == $end;    # i.e. Ns at start or end of scaffold

        my @contig_baits;

        my $baits_needed = int( ( $end - $start ) / $bait_spacing );
        if ( $baits_needed <= 1 ) {

            # Just one bait needed
            my $midpoint = int( ( $end - $start ) / 2 - $BAIT_LENGTH / 2 );
            my ( $bait_start, $bait_end, $bait_seq ) =
              design_bait( $scaffold_seq, $start, $end, $midpoint,
                $gc_threshold );
            if ( defined $bait_seq ) {
                push @contig_baits,
                  [ $scaffold->id, $bait_start, $bait_end, $bait_seq ];
            }
        }
        else {
            # 2+ baits needed

            # Design bait at far left
            my ( $bait_left_start, $bait_left_end, $bait_left_seq ) =
              design_bait( $scaffold_seq, $start, $end, $start, $gc_threshold );
            if ( defined $bait_left_seq ) {
                push @contig_baits,
                  [
                    $scaffold->id,  $bait_left_start,
                    $bait_left_end, $bait_left_seq
                  ];
                $start =
                  $bait_left_end + 1;    # Adjust start so baits don't overlap
                $baits_needed--;
            }

            # Design bait at far right
            my ( $bait_right_start, $bait_right_end, $bait_right_seq ) =
              design_bait( $scaffold_seq, $start, $end, $end - $BAIT_LENGTH + 1,
                $gc_threshold );
            if ( defined $bait_right_seq ) {
                push @contig_baits,
                  [
                    $scaffold->id,   $bait_right_start,
                    $bait_right_end, $bait_right_seq
                  ];
                $end =
                  $bait_right_start - 1;    # Adjust end so baits don't overlap
                $baits_needed--;
            }

            # Evenly space remaining baits
            while ($baits_needed) {
                my $pos =
                  int( ( $end - $start ) / ( $baits_needed + 1 ) +
                      $start -
                      $BAIT_LENGTH / 2 );
                my ( $bait_start, $bait_end, $bait_seq ) =
                  design_bait( $scaffold_seq, $start, $end, $pos,
                    $gc_threshold );
                if ( defined $bait_seq ) {
                    push @contig_baits,
                      [ $scaffold->id, $bait_start, $bait_end, $bait_seq ];
                    $start =
                      $bait_end + 1;    # Adjust start so baits don't overlap
                    $baits_needed--;
                }
                else {
                    last;
                }
            }
        }

        push @scaffold_baits, @contig_baits;
        if ( !@contig_baits ) {
            $undesignable_contig_count++;
        }
    }

    foreach my $bait ( sort { $a->[1] <=> $b->[1] } @scaffold_baits ) {
        printf "%s\n", join "\t", @{$bait};
    }
    $designed_bait_count += scalar @scaffold_baits;
    if ( !@scaffold_baits ) {
        $undesignable_scaffold_count++;
    }
}

# Final stats
printf {*STDERR}
  "Undesignable scaffold count:\t%d (where %dmer couldn't be designed)\n",
  $undesignable_scaffold_count, $BAIT_LENGTH;
printf {*STDERR}
  "Undesignable contig count:\t%d (where %dmer couldn't be designed)\n",
  $undesignable_contig_count, $BAIT_LENGTH;
printf {*STDERR} "Designed baits count:\t%d\n", $designed_bait_count;

# Design bait
sub design_bait {
    my ( $full_seq, $from, $to, $begin, $max_gc ) = @_;

    if ( $begin < $from ) {
        $begin = $from;
    }

    # No point searching if impossible to design 120mer bait
    my $seq = substr $full_seq, $from - 1, $to - $from + 1;
    if ( $seq !~ m/[AGCT]{$BAIT_LENGTH}/xms ) {
        return;
    }

    my $current_left_start  = $begin;
    my $current_right_start = $begin + 1;

    while ($current_left_start >= $from
        || $current_right_start + $BAIT_LENGTH - 1 <= $to )
    {
        foreach my $bait_start ( $current_left_start, $current_right_start ) {
            next if $bait_start < $from;
            next if $bait_start + $BAIT_LENGTH - 1 > $to;
            $seq = substr $full_seq, $bait_start - 1, $BAIT_LENGTH;
            my $gc = ( $seq =~ tr/GC/GC/ ) / $BAIT_LENGTH;
            if ( $gc < $max_gc && $seq =~ m/\A [AGCT]+ \z/xms ) {
                return $bait_start, $bait_start + $BAIT_LENGTH - 1, $seq;
            }
        }
        $current_left_start--;
        $current_right_start++;
    }

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'ref=s' => \$ref_fasta,
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

make_gap_map_baits.pl

Make GAPmap baits from reference genome in FASTA format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script makes GAPmap baits for pulldown enrichment based on a reference
genome in FASTA format.

=head1 EXAMPLES

    perl make_gap_map_baits.pl \
        --ref ref.fa
        > baits.tsv

=head1 USAGE

    make_gap_map_baits.pl
        [--ref file]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--ref FILE>

FASTA file containing reference genome.

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
