#!/usr/bin/env perl

# PODNAME: convert_gff3_to_gtf.pl
# ABSTRACT: Convert GFF3 to GTF

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-11-17

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Sort::Naturally;

# Default options
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over STDIN
my %fields;
my %parent_of;
my %is_feature;
while ( my $line = <> ) {
    chomp $line;
    next if $line =~ m/\A \#/xms;
    my ( $chr, $source, $type, $start, $end, undef, $strand, undef,
        $attributes ) = split /\t/xms, $line;
    next if $type ne 'gene' && $type ne 'transcript' && $type ne 'exon';
    my %attr;
    foreach my $attr ( split /;/xms, $attributes ) {
        my ( $tag, $value ) = split /=/xms, $attr;
        $attr{$tag} = $value;
    }
    my $id = $attr{ID} || join q{:}, $chr, $source, $type, $start, $end,
      $strand;
    $fields{$id} = [ $chr, $source, $type, $start, $end, $strand ];
    $is_feature{$type}{$id} = 1;
    my $parent = $attr{Parent};
    if ( !$parent && $type eq 'transcript' ) {
        $parent = $attr{geneID};
        if ($parent) {
            $is_feature{gene}{$parent} = 1;
        }
    }
    if ( !$parent && $type eq 'exon' ) {
        $parent = $attr{transcriptID};
        if ($parent) {
            $is_feature{transcript}{$parent} = 1;
        }
    }
    if ($parent) {
        $parent_of{$id} = $parent;
    }
}

# Index by parent and sort along chromosome
my %child_of;
foreach my $id ( keys %parent_of ) {
    push @{ $child_of{ $parent_of{$id} } }, $id;
}
foreach my $id ( keys %child_of ) {
    ## no critic (ProhibitMagicNumbers)
    @{ $child_of{$id} } = sort {
             ncmp( $fields{$a}->[0], $fields{$b}->[0] )
          || $fields{$a}->[3] <=> $fields{$b}->[3]
          || $fields{$a}->[4] <=> $fields{$b}->[4]
    } @{ $child_of{$id} };
    ## use critic
}

# Add any missing features
foreach my $type (qw(transcript gene)) {
    my @features = keys %{ $is_feature{$type} };
    foreach my $id (@features) {
        next if exists $fields{$id};
        my ( $start, $end );
        my ( $chr, $source, undef, $min_start, $max_end, $strand );
        foreach my $child ( @{ $child_of{$id} } ) {
            ( $chr, $source, undef, $start, $end, $strand ) =
              @{ $fields{$child} };
            $min_start = !defined $min_start ? $start : $min_start;
            $min_start = $start < $min_start ? $start : $min_start;
            $max_end   = !defined $max_end   ? $end   : $max_end;
            $max_end   = $end > $max_end     ? $end   : $max_end;
        }
        $fields{$id} = [ $chr, $source, $type, $min_start, $max_end, $strand ];
    }
}

# Output GTF
foreach my $gene (
    ## no critic (ProhibitMagicNumbers)
    sort {
             ncmp( $fields{$a}->[0], $fields{$b}->[0] )
          || $fields{$a}->[3] <=> $fields{$b}->[3]
          || $fields{$a}->[4] <=> $fields{$b}->[4]
    } grep { $fields{$_}->[2] eq 'gene' } keys %fields
    ## use critic
  )
{
    my %gene_attr = ( gene_id => $gene );
    output_line( @{ $fields{$gene} }, %gene_attr );
    foreach my $transcript ( @{ $child_of{$gene} } ) {
        my %transcript_attr = %gene_attr;
        $transcript_attr{transcript_id} = $transcript;
        output_line( @{ $fields{$transcript} }, %transcript_attr );
        foreach my $exon ( @{ $child_of{$transcript} } ) {
            output_line( @{ $fields{$exon} }, %transcript_attr );
        }
    }
}

sub output_line {    ## no critic (ProhibitManyArgs)
    my ( $chr, $source, $type, $start, $end, $strand, %attr ) = @_;
    my @attributes;
    foreach my $key ( sort keys %attr ) {
        push @attributes, sprintf '%s "%s";', $key, $attr{$key};
    }
    my $attributes = join q{ }, @attributes;
    printf "%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n", $chr, $source, $type, $start,
      $end, $strand, $attributes;
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

convert_gff3_to_gtf.pl

Convert GFF3 to GTF

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script converts GFF3 files to GTF files, including adding missing data.

=head1 EXAMPLES

    perl \
        convert_gff3_to_gtf.pl < file.gff3 > file.gtf

=head1 USAGE

    convert_gff3_to_gtf.pl
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

This software is Copyright (c) 2017 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
