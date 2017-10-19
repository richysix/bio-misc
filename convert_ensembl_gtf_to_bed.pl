#!/usr/bin/env perl

# PODNAME: convert_ensembl_gtf_to_bed.pl
# ABSTRACT: Convert an Ensembl GTF file to BED format

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2014-07-26

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over GTF on STDIN
my @exons;
while ( my $line = <> ) {

    # Skip header
    next if $line =~ m/\A [#]/xms;

    chomp $line;
    my (
        $seq,  undef,   $feature, $start, $end,
        undef, $strand, undef,    $attributes
    ) = split /\t/xms, $line;

    # Ignore everything other than exons
    next if $feature ne 'exon';

    my ( $gene_id, $transcript_id, $exon_id, $exon_ordinal ) =
      parse_attributes($attributes);

    # Print if starting a new transcript
    ## no critic (ProhibitMagicNumbers)
    if ( @exons && $exons[0]->[5] ne $transcript_id ) {
        ## use critic
        print_exons( \@exons );
        @exons = ();
    }

    push @exons,
      [
        $seq,     $start,         $end,     $strand,
        $gene_id, $transcript_id, $exon_id, $exon_ordinal
      ];
}
print_exons( \@exons );

# Parse attributes field and extract specific types
sub parse_attributes {
    my ($attributes) = @_;

    # Remove trailing ;
    $attributes =~ s/;\z//xms;

    my @attributes = split /;\s/xms, $attributes;
    my %attribute;
    foreach my $attribute (@attributes) {
        my ( $type, $value ) = $attribute =~ m/\A (\S+) \s "? (.+?) "? \z/xms;
        $attribute{$type} = $value;
    }

    return $attribute{gene_id}, $attribute{transcript_id}, $attribute{exon_id},
      $attribute{exon_number};
}

# Print exons making up a transcript
sub print_exons {
    my ($exons) = @_;

    my $max_exon_ordinal = $exons->[-1]->[-1];

    foreach my $exon ( @{$exons} ) {
        my ( $seq, $start, $end, $strand, $gene_id, $transcript_id, $exon_id,
            $exon_ordinal )
          = @{$exon};
        my $name = sprintf '%s-%s-%s-%d-%d', $gene_id, $transcript_id, $exon_id,
          $exon_ordinal, $max_exon_ordinal;
        printf "%s\t%d\t%d\t%s\t.\t%s\n", $seq, $start - 1, $end, $name,
          $strand;
    }

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
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

convert_ensembl_gtf_to_bed.pl

Convert an Ensembl GTF file to BED format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a GTF file on STDIN and outputs a BED file with each line
corresponding to an exon.

=head1 EXAMPLES

    perl \
        convert_ensembl_gtf_to_bed.pl \
        < Danio_rerio.Zv9.75.gtf > Danio_rerio.Zv9.75.bed

=head1 USAGE

    convert_ensembl_gtf_to_bed.pl
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

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
