#!/usr/bin/env perl

# PODNAME: run_coseq_rnaseq.pl
# ABSTRACT: Run coseq on RNA-Seq counts

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2019-02-11

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
my $output_dir = 'coseq';
my $counts_file;
my $samples_file;
my $min_clusters = 2;
my $max_clusters = 16;    ## no critic (ProhibitMagicNumbers)
my $interval;
my $r_cmd                  = 'Rscript';
my $memory                 = 4000;        ## no critic (ProhibitMagicNumbers)
my $strip_condition_prefix = 1;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get samples
my %condition_for;
my %is_condition;
my @all_samples;
my %is_sample;
open my $samples_fh, '<', $samples_file;    ## no critic (RequireBriefOpen)
my $header = <$samples_fh>;

while ( my $line = <$samples_fh> ) {
    chomp $line;
    my ( $sample, $condition ) = split /\t/xms, $line;
    $condition_for{$sample}   = $condition;
    $is_condition{$condition} = 1;
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

# Get counts
my @genes;
my %counts_for;
open my $counts_fh, '<', $counts_file;    ## no critic (RequireBriefOpen)
$header = <$counts_fh>;
chomp $header;
$header =~ s/ \s+ count //xmsg;
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

my @ranges = ( [ $min_clusters, $max_clusters ] );
if ($interval) {
    my $start_clusters = $min_clusters;
    while ( $start_clusters < $max_clusters - $interval + 1 ) {
        push @ranges, [ $start_clusters, $start_clusters + $interval - 1 ];
        $start_clusters += $interval;
        if ( $start_clusters >= $max_clusters - $interval + 1 ) {
            push @ranges, [ $max_clusters - $interval + 1, $max_clusters ];
        }
    }
}
my @model_transforms =
  ( [qw(kmeans logclr)], [qw(kmeans clr)], [qw(Normal arcsin)], );

foreach my $range (@ranges) {
    my ( $min, $max ) = @{$range};
    foreach my $model_transform (@model_transforms) {
        my ( $model, $transform ) = @{$model_transform};
        my $signature = join q{-}, $model, $transform, $min, $max;
        my $dir       = File::Spec->catdir( $output_dir, $signature );
        next if -e ( $dir . '.done' );
        make_path($dir);

        # Write new samples file
        my $new_samples_file = File::Spec->catfile( $dir, 'samples.txt' );
        open $samples_fh, '>', $new_samples_file;
        ## no critic (RequireCheckedSyscalls)
        print {$samples_fh} "\tcondition\n";
        ## use critic
        foreach my $sample (@all_samples) {
            printf {$samples_fh} "%s\t%s\n", $sample, $condition_for{$sample};
        }
        close $samples_fh;

        # Write new counts file
        my $new_counts_file = File::Spec->catfile( $dir, 'counts.txt' );
        open $counts_fh, '>', $new_counts_file;  ## no critic (RequireBriefOpen)
        my @header_samples;
        foreach my $sample (@all_samples) {
            my $condition = $condition_for{$sample};
            push @header_samples, $sample;
        }
        printf {$counts_fh} "\t%s\n", ( join "\t", @header_samples );
        foreach my $i ( 0 .. ( scalar @genes ) - 1 ) {
            my @counts;
            foreach my $sample (@all_samples) {
                push @counts, $counts_for{$sample}->[$i];
            }
            printf {$counts_fh} "%s\t%s\n", $genes[$i], ( join "\t", @counts );
        }
        close $counts_fh;

        # Write R script
        ## no critic (RequireBriefOpen)
        open my $r_fh, '>', File::Spec->catfile( $dir, 'coseq.R' );
        ## use critic
        print {$r_fh} <<"EOF";    ## no critic (RequireCheckedSyscalls)
suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(coseq))
countData <- read.table( "$new_counts_file", header=TRUE, row.names=1, check.names=FALSE )
samples <- read.table( "$new_samples_file", header=TRUE, row.names=1 )
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ 1)
dds <- estimateSizeFactors(dds)
write.table(sizeFactors(dds), file="$dir/size-factors.txt", col.names=FALSE, quote=FALSE, sep="\\t")
write.table(counts(dds, normalized=TRUE), file="$dir/normalised-counts.txt", col.names=FALSE, quote=FALSE, sep="\\t")
run <- coseq(counts(dds, normalized=TRUE), K=$min:$max, model="$model", transformation="$transform", normFactors="none")
summary(run)
write.table(clusters(run), file="$dir/clusters.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\\t" )
pdf("$dir/clusters.pdf")
plot(run, conds=samples\$condition, collapse_reps="average")
dev.off()
write.table(data.frame(), file="$dir.done", col.names=FALSE)
quit()
EOF
        close $r_fh;

        # Run R script under LSF
        printf "Running %s\n", "$dir/coseq.R";
        my $cmd = <<"EOF";
bsub -q long -o $dir/coseq.o -e $dir/coseq.e -R'select[mem>$memory] rusage[mem=$memory]' -M$memory "$r_cmd $dir/coseq.R"
EOF
        WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";
    }
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
        'min_clusters=i'          => \$min_clusters,
        'max_clusters=i'          => \$max_clusters,
        'interval=i'              => \$interval,
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

    ## no critic (ProhibitMagicNumbers)
    if ( $max_clusters < $min_clusters + 10 ) {
        pod2usage(
            "--max_clusters must be at least 10 larger than --min_clusters\n");
    }

    if ( $interval && $interval <= 10 ) {
        pod2usage("--interval must be greater than 10\n");
    }
    ## use critic

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

run_coseq_rnaseq.pl

Run coseq on RNA-seq counts

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes an RNA-Seq counts file and samples file and runs coseq.

=head1 EXAMPLES

    perl \
        run_coseq_rnaseq.pl \
        --output_dir coseq

    perl \
        run_coseq_rnaseq.pl \
        --output_dir coseq \
        --counts_file counts.txt --samples_file samples.txt \
        --min_clusters 2 --max_clusters 32 \
        --interval 11 \
        --memory 8000

=head1 USAGE

    run_coseq_rnaseq.pl
        [--output_dir dir]
        [--counts_file file]
        [--samples_file file]
        [--min_clusters int]
        [--max_clusters int]
        [--interval int]
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

RNA-Seq counts file (e.g. counts.txt) or all file (e.g. all.tsv).

=item B<--samples_file FILE>

DESeq2 samples file (e.g. samples.txt).

=item B<--min_clusters INT>

Minimum number of clusters (defaults to 2).

=item B<--max_clusters INT>

Maximum number of clusters (defaults to 16).

=item B<--interval INT>

Run additional smaller intervals of specified size (minimum 11) between minimum
and maximum number of clusters.

=item B<--r_cmd COMMAND>

Command for running R (defaults to Rscript).

=item B<--memory INT>

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

This software is Copyright (c) 2019 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
