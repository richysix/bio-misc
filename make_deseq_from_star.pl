#!/usr/bin/env perl

# PODNAME: make_deseq_from_star.pl
# ABSTRACT: Make DESeq2 counts and samples files from STAR output

## Author     : ims@iansealy.com
## Maintainer : ims@iansealy.com
## Created    : 2019-09-13

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;
use File::Spec;
use Path::Tiny;

# Constants
Readonly our $COUNTS_FILE  => 'counts.txt';
Readonly our $SAMPLES_FILE => 'samples.txt';
Readonly our %COUNT_COL    => (
    unstranded   => 1,
    firststrand  => 2,
    secondstrand => 3,
);

# Default options
my @count_files;
my $strandedness = 'secondstrand';
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# If only one count file then assume it's a list of count files
if ( scalar @count_files == 1 ) {
    my $count_files = path( $count_files[0] )->slurp;
    chomp $count_files;
    @count_files = split /[\r\n]+/xms, $count_files;
}

# Get all sample names
my @samples = remove_common_prefix_suffix(@count_files);

# Get filehandles for all count files
my %fh_for;
foreach my $i ( 0 .. ( scalar @count_files ) - 1 ) {
    open my $fh, '<', $count_files[$i];    ## no critic (RequireBriefOpen)
    $fh_for{ $samples[$i] } = $fh;
}

# Output counts file
open my $counts_fh, '>', $COUNTS_FILE;     ## no critic (RequireBriefOpen)
printf {$counts_fh} "\t%s\n", join "\t", @samples;
GENE: while (1) {
    my $current_gene;
    my @counts;
    foreach my $sample (@samples) {
        my $line = readline $fh_for{$sample};
        last GENE if !defined $line;
        chomp $line;
        my @fields = split /\t/xms, $line;
        my $gene   = $fields[0];
        my $count  = $fields[ $COUNT_COL{$strandedness} ];
        if ( defined $current_gene && $current_gene ne $gene ) {
            confess sprintf 'Gene order inconsistent between files (%s vs %s)',
              $current_gene, $gene;
        }
        $current_gene = $gene;
        push @counts, $count;
    }
    next if $current_gene =~ m/\A N_/xms;
    printf {$counts_fh} "%s\t%s\n", $current_gene, join "\t", @counts;
}
close $counts_fh;

# Output samples file
open my $samples_fh, '>', $SAMPLES_FILE;
print {$samples_fh} "\tcondition\n";    ## no critic (RequireCheckedSyscalls)
foreach my $sample (@samples) {
    my $condition = $sample;
    $condition =~ s/(_rep)?_?\d+ \z//xms;
    printf {$samples_fh} "%s\t%s\n", $sample, $condition;
}
close $samples_fh;

# Close all filehandles
foreach my $sample (@samples) {
    close $fh_for{$sample};
}

sub remove_common_prefix_suffix {
    my @strings = @_;

    foreach ( 0 .. 1 ) {
        @strings = map { scalar reverse } @strings;
        my $common_idx = 0;
      CHAR: while (1) {
            my $char = substr $strings[0], $common_idx, 1;
            foreach my $string (@strings) {
                last CHAR if substr( $string, $common_idx, 1 ) ne $char;
            }
            $common_idx++;
        }
        if ($common_idx) {
            @strings = map { substr $_, $common_idx } @strings;
        }
    }

    return @strings;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'count_files=s@{1,}' => \@count_files,
        'strandedness=s'     => \$strandedness,
        'debug'              => \$debug,
        'help'               => \$help,
        'man'                => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !@count_files ) {
        pod2usage("--count_files must be specified\n");
    }

    if (   defined $strandedness
        && $strandedness ne 'unstranded'
        && $strandedness ne 'firststrand'
        && $strandedness ne 'secondstrand' )
    {
        pod2usage(
            "--strandedness must be unstranded, firststrand or secondstrand\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

make_deseq_from_star.pl

Make DESeq2 counts and samples files from STAR output

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a list of STAR ReadsPerGene files and creates merged
counts.txt and samples.txt files for DESeq2 in the current directory.

=head1 EXAMPLES

    perl \
        make_deseq_from_star.pl \
            --count_files \
                21somites_1/ReadsPerGene.out.tab \
                21somites_2/ReadsPerGene.out.tab

    perl \
        make_deseq_from_star.pl --count_files somites.fofn \
            --strandedness unstranded

=head1 USAGE

    make_deseq_from_star.pl
        [--count_files file...]
        [--strandedness unstranded|firststrand|secondstrand]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--count_files FILE...>

STAR ReadsPerGene files. If just one file is specified then it's assumed to be
a file of filenames.

=item B<--strandedness unstranded|firststrand|secondstrand>

Strandedness column to use from STAR ReadsPerGene files. Defaults to
secondstrand.

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

Ian Sealy <ims@iansealy.com>

=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2019 by Ian Sealy.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
