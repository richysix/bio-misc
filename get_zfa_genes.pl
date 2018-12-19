#!/usr/bin/env perl

# PODNAME: get_zfa_genes.pl
# ABSTRACT: Get genes driving ZFA enrichment

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2018-10-22

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;

# Constants
Readonly our $SHORT_OUTPUT_FILE  => 'zfa.sig.tsv';
Readonly our $MEDIUM_OUTPUT_FILE => 'zfa.sig.genes.sig.tsv';
Readonly our $LONG_OUTPUT_FILE   => 'zfa.sig.genes.tsv';
Readonly our $SIG_LEVEL          => 0.05;

# Default options
my $all_file;
my $sig_file;
my $ensembl_zfin_file;
my $annotations_file;
my $ontologizer_file = 'table-sig-zfin-genes-Parent-Child-Union-Bonferroni.txt';

my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get all Ensembl genes and some associated info
my %is_ens_gene;
my %gene_pval_for;
my %gene_adj_pval_for;
my %gene_log2fc_for;
my %gene_name_for;
my %gene_description_for;
open my $all_fh, '<', $all_file;    ## no critic (RequireBriefOpen)
while ( my $line = <$all_fh> ) {
    chomp $line;
    my (
        $gene, $pval, $adj_pval, $log2fc, undef, undef,
        undef, undef, undef,     $name,   $description
    ) = split /\t/xms, $line;
    next if $gene !~ m/\A ENS[[:upper:]]*G\d{11} \z/xms;
    $is_ens_gene{$gene}          = 1;
    $gene_pval_for{$gene}        = $pval;
    $gene_adj_pval_for{$gene}    = $adj_pval;
    $gene_log2fc_for{$gene}      = $log2fc;
    $gene_name_for{$gene}        = $name;
    $gene_description_for{$gene} = $description;
}
close $all_fh;
confess "No Ensembl IDs in first column of $all_file\n"
  if !scalar keys %is_ens_gene;

# Get significant Ensembl genes
my %is_sig_ens_gene;
open my $sig_fh, '<', $sig_file;
while ( my $gene = <$sig_fh> ) {
    $gene =~ s/\s+.*//xms;
    next if $gene !~ m/\A ENS[[:upper:]]*G\d{11} \z/xms;
    $is_sig_ens_gene{$gene} = 1;
}
close $sig_fh;
confess "No Ensembl IDs in first column of $sig_file\n"
  if !scalar keys %is_sig_ens_gene;

# Link ZFIN IDs to Ensembl IDs
my %zfin2ens;
open my $ensembl_zfin_fh, '<', $ensembl_zfin_file;
while ( my $line = <$ensembl_zfin_fh> ) {
    chomp $line;
    my ( $ensembl_id, $zfin_id ) = split /\s+/xms, $line;
    next if $ensembl_id !~ m/\A ENS[[:upper:]]*G\d{11} \z/xms;
    next if $zfin_id !~ m/\A ZDB-GENE-\d{6}-\d+ \z/xms;
    $zfin2ens{$zfin_id}{$ensembl_id} = 1;
}
close $ensembl_zfin_fh;
confess "No Ensembl IDs or ZFIN IDs in $ensembl_zfin_file\n"
  if !scalar keys %zfin2ens;

# Link ZFA terms to ZFIN IDs
my %zfa2zfin;
open my $annotations_fh, '<', $annotations_file;
while ( my $line = <$annotations_fh> ) {
    my ( undef, $zfin_id, undef, $zfa_term ) = split /\s+/xms, $line;
    next if $zfin_id !~ m/\A ZDB-GENE-\d{6}-\d+ \z/xms;
    next if $zfa_term !~ m/\A ZFA:\d{7} \z/xms;
    $zfa2zfin{$zfa_term}{$zfin_id} = 1;
}
close $annotations_fh;
confess "No ZFIN IDs or ZFA terms in $annotations_file\n"
  if !scalar keys %zfa2zfin;

# Get ZFA terms
my %name_for;
my %pval_for;
my %adj_pval_for;
open my $ontologizer_fh, '<', $ontologizer_file; ## no critic (RequireBriefOpen)
while ( my $line = <$ontologizer_fh> ) {
    chomp $line;
    my (
        $zfa_term, undef,    undef, undef,     undef, undef, undef,
        undef,     $trivial, $pval, $adj_pval, undef, $name
    ) = split /\t/xms, $line;
    next if $trivial ne 'false';
    next if $zfa_term !~ m/\A ZFA:\d{7} \z/xms;
    next if $adj_pval > $SIG_LEVEL;
    $name =~ s/\A "//xms;
    $name =~ s/\" \z//xms;
    $name_for{$zfa_term}     = $name;
    $pval_for{$zfa_term}     = $pval;
    $adj_pval_for{$zfa_term} = $adj_pval;
}
close $ontologizer_fh;

## no critic (RequireBriefOpen)
open my $short_fh,  '>', $SHORT_OUTPUT_FILE;
open my $medium_fh, '>', $MEDIUM_OUTPUT_FILE;
open my $long_fh,   '>', $LONG_OUTPUT_FILE;
## use critic
printf {$short_fh} "%s\n", join "\t", 'ID', 'p-value', 'Adjusted p-value',
  'Name', 'Sig Genes', 'Other Genes';
printf {$medium_fh} "%s\n", join "\t", 'ID', 'p-value', 'Adjusted p-value',
  'Name', 'Gene', 'Gene Adjusted p-value', 'log2fc', 'Name', 'Description';
printf {$long_fh} "%s\n", join "\t", 'ID', 'p-value', 'Adjusted p-value',
  'Name', 'Gene', 'Gene Adjusted p-value', 'log2fc', 'Name', 'Description';

foreach my $zfa_term (
    sort {
             $adj_pval_for{$a} <=> $adj_pval_for{$b}
          || $pval_for{$a} <=> $pval_for{$b}
    } keys %adj_pval_for
  )
{
    my @zfin_ids = keys %{ $zfa2zfin{$zfa_term} };
    my @ensembl_ids =
      get_valid_ensembl_ids_from_zfin_ids( \%zfin2ens, \%is_ens_gene,
        @zfin_ids );
    my $sig_ensembl_ids = join q{,}, grep { $is_sig_ens_gene{$_} } @ensembl_ids;
    my $other_ensembl_ids = join q{,},
      grep { !$is_sig_ens_gene{$_} } @ensembl_ids;

    printf {$short_fh} "%s\n", join "\t", $zfa_term, $pval_for{$zfa_term},
      $adj_pval_for{$zfa_term}, $name_for{$zfa_term}, $sig_ensembl_ids || q{-},
      $other_ensembl_ids || q{-};

    foreach my $ensembl_id (
        sort {
            sort_ensembl_ids( $a, $b, \%gene_adj_pval_for, \%gene_pval_for )
        } @ensembl_ids
      )
    {
        printf {$long_fh} "%s\n", join "\t", $zfa_term, $pval_for{$zfa_term},
          $adj_pval_for{$zfa_term}, $name_for{$zfa_term}, $ensembl_id,
          $gene_adj_pval_for{$ensembl_id}, $gene_log2fc_for{$ensembl_id},
          $gene_name_for{$ensembl_id},     $gene_description_for{$ensembl_id};
        next
          if $gene_adj_pval_for{$ensembl_id} eq 'NA'
          || $gene_adj_pval_for{$ensembl_id} > $SIG_LEVEL;
        printf {$medium_fh} "%s\n", join "\t", $zfa_term, $pval_for{$zfa_term},
          $adj_pval_for{$zfa_term}, $name_for{$zfa_term}, $ensembl_id,
          $gene_adj_pval_for{$ensembl_id}, $gene_log2fc_for{$ensembl_id},
          $gene_name_for{$ensembl_id},     $gene_description_for{$ensembl_id};
    }
}
close $short_fh;
close $medium_fh;
close $long_fh;

sub get_valid_ensembl_ids_from_zfin_ids {
    my ( $zfin2ens, $is_ens_gene, @zfin_ids ) = @_;

    my %ens;
    foreach my $zfin_id (@zfin_ids) {
        foreach my $ensembl_id ( keys %{ $zfin2ens->{$zfin_id} } ) {
            $ens{$ensembl_id} = 1;
        }
    }
    my @ensembl_ids = sort grep { $is_ens_gene->{$_} } keys %ens;

    return @ensembl_ids;
}

sub sort_ensembl_ids {
    my ( $a, $b, $gene_adj_pval_for, $gene_pval_for ) = @_;

    my $pval_a = $gene_adj_pval_for->{$a};
    my $pval_b = $gene_adj_pval_for->{$b};
    if ( $pval_a eq $pval_b ) {
        $pval_a = $gene_pval_for->{$a};
        $pval_b = $gene_pval_for->{$b};
    }

    if ( $pval_a eq 'NA' && $pval_b eq 'NA' ) {
        return $a cmp $b;
    }
    elsif ( $pval_a eq 'NA' ) {
        return 1;
    }
    elsif ( $pval_b eq 'NA' ) {
        return -1;    ## no critic (ProhibitMagicNumbers)
    }

    return $pval_a <=> $pval_b;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'all_file=s'          => \$all_file,
        'sig_file=s'          => \$sig_file,
        'ensembl_zfin_file=s' => \$ensembl_zfin_file,
        'annotations_file=s'  => \$annotations_file,
        'ontologizer_file=s'  => \$ontologizer_file,
        'debug'               => \$debug,
        'help'                => \$help,
        'man'                 => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$all_file ) {
        pod2usage("--all_file must be specified\n");
    }
    if ( !$sig_file ) {
        pod2usage("--sig_file must be specified\n");
    }
    if ( !$ensembl_zfin_file ) {
        pod2usage("--ensembl_zfin_file must be specified\n");
    }
    if ( !$annotations_file ) {
        pod2usage("--annotations_file must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

get_zfa_genes.pl

Get genes driving ZFA enrichment

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes the inputs  and outputs to a ZFA enrichment with Ontologizer
and gives a list of genes driving each enriched ZFA term.

=head1 EXAMPLES

    perl \
        get_zfa_genes.pl \
        --all_file all.tsv \
        --sig_file sig.tsv \
        --ensembl_zfin_file danio_rerio_e94_zfin.txt \
        --annotations_file annotations.txt

    perl \
        get_zfa_genes.pl \
        --all_file all.tsv \
        --sig_file sig.tsv \
        --ensembl_zfin_file danio_rerio_e94_zfin.txt \
        --annotations_file annotations.txt \
        --ontologizer_file \
            table-sig-zfin-genes-Parent-Child-Union-Bonferroni.txt

=head1 USAGE

    get_zfa_genes.pl
        [--all_file file]
        [--sig_file file]
        [--ensembl_zfin_file file]
        [--annotations_file file]
        [--ontologizer_file file]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--all_file FILE>

File containing list of all Ensembl IDs in the experiment.

=item B<--sig_file FILE>

File containing list of all significant Ensembl IDs in the experiment.

=item B<--ensembl_zfin_file FILE>

File linking Ensembl IDs to ZFIN IDs (for a specific Ensembl version).

=item B<--annotations_file FILE>

File linking ZFIN IDs to ZFA terms.

=item B<--ontologizer_file FILE>

Ontologizer output file listing all ZFA terms.

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
