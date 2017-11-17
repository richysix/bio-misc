#!/usr/bin/env perl

# PODNAME: convert_deseq2_to_gsea.pl
# ABSTRACT: Convert DESeq2 output to GSEA input

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-10-26

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use File::Spec;
use File::Path qw( make_path );

# Default options
my $dir;
my $all_file;
my $samples_file;
my $exp_condition;
my $con_condition;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get samples
my %condition_for;
my %is_condition;
my @all_samples;
open my $samples_fh, '<', $samples_file;
my $header = <$samples_fh>;
while ( my $line = <$samples_fh> ) {
    chomp $line;
    my ( $sample, $condition ) = split /\t/xms, $line;
    $condition_for{$sample}   = $condition;
    $is_condition{$condition} = 1;
    push @all_samples, $sample;
}
close $samples_fh;

# Remove common prefix from conditions
my $prefix =
  scalar keys %is_condition > 1 ? get_common_prefix( keys %is_condition ) : q{};
%is_condition = ();
foreach my $sample ( keys %condition_for ) {
    $condition_for{$sample} =~ s/\A $prefix //xms;
    $is_condition{ $condition_for{$sample} } = 1;
}

# Guess conditions if necessary
if ( ( !$exp_condition || !$con_condition ) && $all_file =~ m/_vs_/xms ) {
    my ( undef, undef, $filename ) = File::Spec->splitpath($all_file);
    ( $exp_condition, $con_condition ) =
      $filename =~ m/\A (\w+)_vs_(\w+)[.]/xms;
}

# Check conditions
if ( !$exp_condition ) {
    confess 'Experimental condition not guessable and must be specified';
}
if ( !$con_condition ) {
    confess 'Control condition not guessable and must be specified';
}
if ( !exists $is_condition{$exp_condition} ) {
    confess "Specified experimental condition ($exp_condition) unknown";
}
if ( !exists $is_condition{$con_condition} ) {
    confess "Specified control condition ($con_condition) unknown";
}

# Remove samples with other conditions
my @samples = grep {
         $condition_for{$_} eq $exp_condition
      || $condition_for{$_} eq $con_condition
} @all_samples;

# Write CLS file
my $cls_file = File::Spec->catfile( $dir, 'samples.cls' );
open $samples_fh, '>', $cls_file;
printf {$samples_fh} "%d 2 1\n", scalar @samples;
printf {$samples_fh} "# %s %s\n", $exp_condition, $con_condition;
my @classes = map { $condition_for{$_} eq $exp_condition ? 0 : 1 } @samples;
printf {$samples_fh} "%s\n", ( join q{ }, @classes );
close $samples_fh;

# Get headings, normalised counts and log2 fold change
open my $all_fh, '<', $all_file;    ## no critic (RequireBriefOpen)
$header = <$all_fh>;
chomp $header;
my @headings = split /\t/xms, $header;
my %sample_to_col;
my $col = -1;                       ## no critic (ProhibitMagicNumbers)
foreach my $heading (@headings) {
    $col++;
    if ( $heading =~ m/\s normalised \s count \z/xms ) {
        $heading =~ s/\s normalised \s count \z//xms;
        $sample_to_col{$heading} = $col;
    }
}
my @genes;
my %counts_for;
my %l2fc_for;
while ( my $line = <$all_fh> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    next if $fields[2] eq 'NA';    # No adjusted p-value
    push @genes, $fields[0];
    foreach my $sample (@samples) {
        push @{ $counts_for{$sample} }, $fields[ $sample_to_col{$sample} ];
    }
    $l2fc_for{ $fields[0] } = $fields[3];    ## no critic (ProhibitMagicNumbers)
}
close $all_fh;

# Write GCT file
my $gct_file = File::Spec->catfile( $dir, 'counts.gct' );
open my $gct_fh, '>', $gct_file;             ## no critic (RequireBriefOpen)
printf {$gct_fh} "%s\n", '#1.2';
printf {$gct_fh} "%d\t%d\n", scalar @genes, scalar @samples;
printf {$gct_fh} "NAME\tDescription\t%s\n", ( join "\t", @samples );
foreach my $i ( 0 .. ( scalar @genes ) - 1 ) {
    my @counts;
    foreach my $sample (@samples) {
        push @counts, $counts_for{$sample}->[$i];
    }
    printf {$gct_fh} "%s\tNA\t%s\n", $genes[$i], ( join "\t", @counts );
}
close $gct_fh;

# Write RNK file
my $rnk_file = File::Spec->catfile( $dir, 'genes.rnk' );
open my $rnk_fh, '>', $rnk_file;
foreach my $gene (@genes) {
    printf {$rnk_fh} "%s\t%s\n", $gene, $l2fc_for{$gene};
}
close $rnk_fh;

sub get_common_prefix {
    my @strings = @_;

    my $common_prefix_idx = 0;
  CHAR: while (1) {
        my $char = substr $strings[0], $common_prefix_idx, 1;
        foreach my $string (@strings) {
            last CHAR if substr( $string, $common_prefix_idx, 1 ) ne $char;
        }
        $common_prefix_idx++;
    }

    my $common_prefix = substr $strings[0], 0, $common_prefix_idx;

    # Hack to stop xx_het and xx_hom and het and hom reducing to et and om
    $common_prefix =~ s/_h \z/_/xms;
    if ( length $common_prefix == 1 ) {
        $common_prefix = q{};
    }

    return $common_prefix;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'           => \$dir,
        'all_file=s'      => \$all_file,
        'samples_file=s'  => \$samples_file,
        'exp_condition=s' => \$exp_condition,
        'con_condition=s' => \$con_condition,
        'debug'           => \$debug,
        'help'            => \$help,
        'man'             => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    $dir =~ s/\/ \z//xms;

    if ( !$all_file ) {
        $all_file = $dir . '.tsv';
    }
    if ( !$samples_file ) {
        $samples_file = File::Spec->catfile( $dir, 'samples.txt' );
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

convert_deseq2_to_gsea.pl

Convert DESeq2 output to GSEA input

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes an RNA-Seq output file and samples file and converts them for
use with GSEA.

=head1 EXAMPLES

    perl \
        convert_deseq2_to_gsea.pl \
        --dir hom_vs_wt

    perl \
        convert_deseq2_to_gsea.pl \
        --dir hom_vs_wt \
        --all_file hom_vs_wt.tsv \
        --samples_file hom_vs_wt/samples.txt \
        --exp_condition hom \
        --con_condition wt

=head1 USAGE

    convert_deseq2_to_gsea.pl
        [--dir dir]
        [--all_file file]
        [--samples_file file]
        [--exp_condition condition]
        [--con_condition condition]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--dir DIR>

Directory in which to create output files (and which should become the RNA-seq
output file if .tsv is appended and should contain the samples.txt file if
neither is specified explicitly).

=item B<--all_file FILE>

RNA-Seq output file (e.g. all.tsv).

=item B<--samples_file FILE>

DESeq2 samples file (e.g. samples.txt). The order of samples in the samples file
determines the order of the columns in the output.

=item B<--exp_condition CONDITION>

Experimental condition from sample files. Can be guessed from all file name if
in form x_vs_y.tsv.

=item B<--con_condition CONDITION>

Control condition from sample files. Can be guessed from all file name if in
form x_vs_y.tsv.

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
