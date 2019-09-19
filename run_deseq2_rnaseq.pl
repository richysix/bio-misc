#!/usr/bin/env perl

# PODNAME: run_deseq2_rnaseq.pl
# ABSTRACT: Run DESeq2 on RNA-Seq counts

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-05-16

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use English qw( -no_match_vars );
use POSIX qw( WIFEXITED);
use File::Spec;
use File::Path qw( make_path );

# Default options
my $output_dir = 'deseq2';
my $counts_file;
my $samples_file;
my @comparisons;
my $remove_other_conditions;
my $interaction;
my $lfc_threshold          = 0;
my $r_cmd                  = 'Rscript';
my $memory                 = 4000;        ## no critic (ProhibitMagicNumbers)
my $strip_condition_prefix = 1;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get samples
my %condition_for;
my %is_condition;
my %group_for;
my @all_samples;
my %is_sample;
open my $samples_fh, '<', $samples_file;    ## no critic (RequireBriefOpen)
my $header = <$samples_fh>;
chomp $header;
my ( undef, undef, @group_names ) = split /\t/xms, $header;

while ( my $line = <$samples_fh> ) {
    chomp $line;
    my ( $sample, $condition, @groups ) = split /\t/xms, $line;
    $condition_for{$sample}   = $condition;
    $is_condition{$condition} = 1;
    my $i = 0;
    foreach my $group_name (@group_names) {
        $group_for{$group_name}{$sample} = $groups[$i];
        $i++;
    }
    push @all_samples, $sample;
    $is_sample{$sample} = 1;
}
close $samples_fh;

# Remove common prefix from conditions
if ($strip_condition_prefix) {
    my $prefix =
      scalar keys %is_condition > 1
      ? get_common_prefix( keys %is_condition )
      : q{};
    %is_condition = ();
    foreach my $sample ( keys %condition_for ) {
        $condition_for{$sample} =~ s/\A $prefix //xms;
        $is_condition{ $condition_for{$sample} } = 1;
    }
}

my @all_conditions = sort keys %is_condition;

# Get counts
my @genes;
my %counts_for;
open my $counts_fh, '<', $counts_file;    ## no critic (RequireBriefOpen)
$header = <$counts_fh>;
chomp $header;
my ( undef, @count_samples ) = split /\t/xms, $header;
while ( my $line = <$counts_fh> ) {
    chomp $line;
    my ( $gene, @counts ) = split /\t/xms, $line;
    push @genes, $gene;
    foreach my $i ( 0 .. ( scalar @count_samples ) - 1 ) {
        next if !$is_sample{ $count_samples[$i] };
        push @{ $counts_for{ $count_samples[$i] } }, $counts[$i];
    }
}

# Assume some comparisons
if ( !@comparisons ) {
    if ( scalar @all_conditions == 1 ) {
        confess "Only one condition (@all_conditions) for $output_dir";
    }

    my ($wt)   = grep { m/(\b|_)wt   \z/xms } @all_conditions;
    my ($het)  = grep { m/(\b|_)het  \z/xms } @all_conditions;
    my ($hom)  = grep { m/(\b|_)hom  \z/xms } @all_conditions;
    my ($hemi) = grep { m/(\b|_)hemi \z/xms } @all_conditions;
    my ($sib)  = grep { m/(\b|_)sib  \z/xms } @all_conditions;
    my ($mut)  = grep { m/(\b|_)mut  \z/xms } @all_conditions;

    if ( $wt && $het && $hom ) {
        push @comparisons, "$het:$wt";
        push @comparisons, "$hom:$wt";
        push @comparisons, "$hom:$het";
        push @comparisons, "$hom:$het,$wt";
        push @comparisons, "$hom,$het:$wt";
    }
    if ( $wt && $het && $hemi ) {
        push @comparisons, "$het:$wt";
        push @comparisons, "$hemi:$wt";
        push @comparisons, "$hemi:$het";
        push @comparisons, "$hemi:$het,$wt";
        push @comparisons, "$hemi,$het:$wt";
    }
    if ( $wt && $het && $hom && $hemi ) {
        push @comparisons, "$het:$wt";
        push @comparisons, "$hom:$wt";
        push @comparisons, "$hemi:$wt";
        push @comparisons, "$hom:$het";
        push @comparisons, "$hemi:$het";
        push @comparisons, "$hom:$het,$wt";
        push @comparisons, "$hom,$het:$wt";
        push @comparisons, "$hemi:$het,$wt";
        push @comparisons, "$hemi,$het:$wt";
        push @comparisons, "$hom,$hemi:$het,$wt";
    }
    if ( !@comparisons ) {
        @comparisons =
            $wt  && $het  ? ("$het:$wt")
          : $wt  && $hom  ? ("$hom:$wt")
          : $wt  && $hemi ? ("$hemi:$wt")
          : $het && $hom  ? ("$hom:$het")
          : $het && $hemi ? ("$hemi:$het")
          : $hom && $hemi ? ("$hom:$hemi")
          : $sib && $mut  ? ("$mut:$sib")
          :                 ();
    }
}

foreach my $comparison (@comparisons) {
    my ( $all_exp, $all_con ) = split /:/xms, $comparison;
    confess "Experimental condition missing from $comparison for $output_dir"
      if !$all_exp;
    confess "Control condition missing from $comparison for $output_dir"
      if !$all_con;
    my ( $exp, $exp_name ) = split /=/xms, $all_exp;
    my ( $con, $con_name ) = split /=/xms, $all_con;
    if ( !$exp_name ) {
        $exp_name = $exp;
        $exp_name =~ s/,/_/xmsg;
    }
    if ( !$con_name ) {
        $con_name = $con;
        $con_name =~ s/,/_/xmsg;
    }
    my @exp = split /,/xms, $exp;
    my @con = split /,/xms, $con;
    my %rename;
    foreach my $condition (@exp) {
        confess "Unknown condition ($condition) in $comparison for $output_dir"
          if !$is_condition{$condition};
        $rename{$condition} = $exp_name;
    }
    foreach my $condition (@con) {
        confess "Unknown condition ($condition) in $comparison for $output_dir"
          if !$is_condition{$condition};
        $rename{$condition} = $con_name;
    }
    my $dir = File::Spec->catdir( $output_dir, $exp_name . '_vs_' . $con_name );
    if ( $interaction && @group_names ) {
        $dir .= '-interaction';
    }
    next if -e ( $dir . '.done' );
    make_path($dir);

    # Write new samples file
    my $new_samples_file = File::Spec->catfile( $dir, 'samples.txt' );
    open $samples_fh, '>', $new_samples_file;    ## no critic (RequireBriefOpen)
    ## no critic (RequireCheckedSyscalls)
    print {$samples_fh} "\tcondition";
    foreach my $group_name (@group_names) {
        print {$samples_fh} "\t" . $group_name;
    }
    print {$samples_fh} "\n";
    foreach my $sample (@all_samples) {
        my $condition = $condition_for{$sample};
        next if $remove_other_conditions && !exists $rename{$condition};

        printf {$samples_fh} "%s\t%s", $sample,
          exists $rename{ $condition_for{$sample} }
          ? $rename{ $condition_for{$sample} }
          : $condition_for{$sample};
        foreach my $group_name (@group_names) {
            print {$samples_fh} "\t" . $group_for{$group_name}{$sample};
        }
        print {$samples_fh} "\n";
    }
    ## use critic
    close $samples_fh;

    # Write new counts file
    my $new_counts_file = File::Spec->catfile( $dir, 'counts.txt' );
    open $counts_fh, '>', $new_counts_file;    ## no critic (RequireBriefOpen)
    my @header_samples;
    foreach my $sample (@all_samples) {
        my $condition = $condition_for{$sample};
        next if $remove_other_conditions && !exists $rename{$condition};
        push @header_samples, $sample;
    }
    printf {$counts_fh} "\t%s\n", ( join "\t", @header_samples );
    foreach my $i ( 0 .. ( scalar @genes ) - 1 ) {
        my @counts;
        foreach my $sample (@all_samples) {
            my $condition = $condition_for{$sample};
            next if $remove_other_conditions && !exists $rename{$condition};
            push @counts, $counts_for{$sample}->[$i];
        }
        printf {$counts_fh} "%s\t%s\n", $genes[$i], ( join "\t", @counts );
    }
    close $counts_fh;

    # Write R script
    my $design = 'condition';
    if (@group_names) {
        $design = sprintf '%s + %s', ( join ' + ', @group_names ), $design;
    }
    if ( $interaction && @group_names ) {
        foreach my $group_name (@group_names) {
            $design .= sprintf ' + %s:condition', $group_name;
        }
    }
    ## no critic (RequireBriefOpen)
    open my $r_fh, '>', File::Spec->catfile( $dir, 'deseq2.R' );
    ## use critic
    print {$r_fh} <<"EOF";    ## no critic (RequireCheckedSyscalls)
suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gplots))
countData <- read.table( "$new_counts_file", header=TRUE, row.names=1 )
samples <- read.table( "$new_samples_file", header=TRUE, row.names=1 )
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ $design)
dds <- DESeq(dds)
write.table(sizeFactors(dds), file="$dir/size-factors.txt", col.names=FALSE, quote=FALSE, sep="\\t")
write.table(counts(dds, normalized=TRUE), file="$dir/normalised-counts.txt", col.names=FALSE, quote=FALSE, sep="\\t")
res <- results(dds, contrast=c("condition", "$exp_name", "$con_name"), lfcThreshold=$lfc_threshold)
out <- data.frame(pvalue=res\$pvalue, padj=res\$padj, log2fc=res\$log2FoldChange, row.names=rownames(res))
write.table(out, file="$dir/output.txt", col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\\t")
pdf("$dir/qc.pdf")
plotMA(dds)
rld <- rlogTransformation(dds, blind=TRUE)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(as.matrix(dist(t(assay(rld)))), trace="none", col = rev(hmcol), margin=c(13, 13))
print(plotPCA(rld, intgroup=c("condition")))
plotDispEsts(dds)
dev.off()
write.table(data.frame(), file="$dir.done", col.names=FALSE)
quit()
EOF
    close $r_fh;

    # Run R script under LSF
    printf "Running %s\n", "$dir/deseq2.R";
    my $cmd = <<"EOF";
bsub -o $dir/deseq2.o -e $dir/deseq2.e -R'select[mem>$memory] rusage[mem=$memory]' -M$memory "$r_cmd $dir/deseq2.R"
EOF
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";
}

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
        'output_dir=s'            => \$output_dir,
        'counts_file=s'           => \$counts_file,
        'samples_file=s'          => \$samples_file,
        'comparisons=s@{,}'       => \@comparisons,
        'remove_other_conditions' => \$remove_other_conditions,
        'interaction'             => \$interaction,
        'lfc_threshold=i'         => \$lfc_threshold,
        'r_cmd=s'                 => \$r_cmd,
        'memory=i'                => \$memory,
        'strip_condition_prefix!' => \$strip_condition_prefix,
        'debug'                   => \$debug,
        'help'                    => \$help,
        'man'                     => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$counts_file ) {
        $counts_file = File::Spec->catfile( $output_dir, 'counts.txt' );
    }
    if ( !$samples_file ) {
        $samples_file = File::Spec->catfile( $output_dir, 'samples.txt' );
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

run_deseq2_rnaseq.pl

Run DESeq2 on RNA-seq counts

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes an RNA-Seq counts file and samples file and runs DESeq2 using
specified or all conditions.

=head1 EXAMPLES

    perl \
        run_deseq2_rnaseq.pl \
        --output_dir deseq2

    perl \
        run_deseq2_rnaseq.pl \
        --output_dir deseq2 \
        --counts_file counts.txt --samples_file samples.txt \
        --comparisons hom:het hom:wt het:wt hom=mut:het,wt=sib \
        --remove_other_conditions \
        --interaction \
        --lfc_threshold 1 \
        --memory 8000

=head1 USAGE

    run_deseq2_rnaseq.pl
        [--output_dir dir]
        [--counts_file file]
        [--samples_file file]
        [--comparisons comparison...]
        [--remove_other_conditions]
        [--interaction]
        [--lfc_threshold int]
        [--r_cmd command]
        [--memory int]
        [--strip_condition_prefix/--nostrip_condition_prefix]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--output_dir DIR>

Directory in which to create output directories (and which should contain
counts.txt and samples.txt if not specified explicitly).

=item B<--counts_file FILE>

RNA-Seq counts file (e.g. counts.txt).

=item B<--samples_file FILE>

DESeq2 samples file (e.g. samples.txt). The order of samples in the samples file
determines the order of the columns in the output.

=item B<--comparisons COMPARISONS>

Condition comparisons. Each comparison is a pair of experimental and control
conditions (in that order) separated by a colon (e.g. hom:wt). If multiple
conditions are to be combined then separate them with a comma (e.g. hom:wt,het).
To rename a condition, append an equals sign (e.g. hom=mut:het,wt=sib).

=item B<--remove_other_conditions>

Remove other conditions from the counts file prior to running DESeq2.

=item B<--interaction>

Include an interaction term in multi-factor designs.

=item B<--lfc_threshold>

Log2 fold change threshold for significance.

=item B<--r_cmd COMMAND>

Command for running R (defaults to Rscript).

=item B<--memory>

Memory to allocate for each LSF job.

=item B<--strip_condition_prefix>

Remove common prefix from conditions. This option is on by default.
Supplying the option --nostrip_condition_prefix turns this off.

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
