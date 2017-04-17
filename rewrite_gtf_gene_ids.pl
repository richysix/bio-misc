#!/usr/bin/env perl

# PODNAME: rewrite_gtf_gene_ids.pl
# ABSTRACT: Rewrite gene IDs in GTF to account for overlapping transcripts

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-04-17

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my $gtf_file;
my $gene_prefix = 'GENE_';
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Organise transcripts by chromosome and strand
my %transcript;
open my $fh, '<', $gtf_file;    ## no critic (RequireBriefOpen)
while ( my $line = <$fh> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    ## no critic (ProhibitMagicNumbers)
    my ( $chr, $feature, $start, $end, $strand, $attributes ) =
      @fields[ 0, 2, 3, 4, 6, 8 ];
    ## use critic
    next if $feature ne 'transcript';
    my ($transcript_id) = $attributes =~ m/transcript_id \s "(\w+)"/xms;
    confess "No transcript ID found in: $line" if !$transcript_id;
    push @{ $transcript{$chr}{$strand} }, [ $start, $end, $transcript_id ];
}
close $fh;

# Assign gene IDs
my $next_id = 1;
my %transcript_to_gene;
foreach my $chr ( sort keys %transcript ) {
    foreach my $strand ( sort keys %{ $transcript{$chr} } ) {
        my @transcripts = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
          @{ $transcript{$chr}{$strand} };
        my $current_gene_id = $next_id++;
        my $current_start   = $transcripts[0][0];
        my $current_end     = $transcripts[0][1];
        $transcript_to_gene{ $transcripts[0][2] } = $current_gene_id;
        shift @transcripts;
        while ( my $transcript = shift @transcripts ) {
            my ( $start, $end, $transcript_id ) = @{$transcript};
            if ( $start > $current_end ) {

                # New gene
                $current_gene_id = $next_id++;
                $current_start   = $start;
                $current_end     = $end;
            }
            else {
                # Same gene
                $current_end = $end;
            }
            $transcript_to_gene{$transcript_id} = $current_gene_id;
        }
    }
}

# Rewrite GTF
my $zero_pad = length --$next_id;
open $fh, '<', $gtf_file;    ## no critic (RequireBriefOpen)
while ( my $line = <$fh> ) {
    chomp $line;
    my ($transcript_id) = $line =~ m/transcript_id \s "(\w+)"/xms;
    confess "No transcript ID found in: $line" if !$transcript_id;
    my $gene_id = sprintf "$gene_prefix%0${zero_pad}d",
      $transcript_to_gene{$transcript_id};
    confess "No gene ID for transcript ID ($transcript_id)" if !$gene_id;
    $line =~ s/gene_id \s "$transcript_id"/gene_id "$gene_id"/xms;
    printf "%s\n", $line;
}
close $fh;

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'gtf_file=s'    => \$gtf_file,
        'gene_prefix=s' => \$gene_prefix,
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

    if ( !$gtf_file ) {
        pod2usage("--gtf_file must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

rewrite_gtf_gene_ids.pl

Rewrite gene IDs in GTF to account for overlapping transcripts

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script rewrites gene IDs in a GTF file and merges overlapping transcripts
into a single gene.

=head1 EXAMPLES

    perl \
        rewrite_gtf_gene_ids.pl --gtf_file in.gtf > out.gtf

=head1 USAGE

    rewrite_gtf_gene_ids.pl
        [--gtf_file file]
        [--gene_prefix prefix]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--gtf_file FILE>

GTF file to be rewritten.

=item B<--gene_prefix PREFIX>

Prefix to be prepended to ordinals to make gene IDs.

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
