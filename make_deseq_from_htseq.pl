#!/usr/bin/env perl

# PODNAME: make_deseq_from_htseq.pl
# ABSTRACT: Make DESeq2 counts and samples files from htseq-count files

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-05-23

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;
use DBI;
use File::Spec;
use Path::Tiny;

# Constants
Readonly our $COUNTS_FILE  => 'counts.txt';
Readonly our $SAMPLES_FILE => 'samples.txt';
Readonly our $HOST         => 'seqw-db';
Readonly our $PORT         => '3379';
Readonly our $USER         => 'warehouse_ro';
Readonly our $PASS         => q{};
Readonly our $NAME         => 'sequencescape_warehouse';

# Default options
my @count_files;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connect to database
my $dsn = "dbi:mysql:database=$NAME;host=$HOST;port=$PORT";
my $dbh = DBI->connect( $dsn, $USER, $PASS );

# Cache sample names
my %name_for;
my $ary_ref = $dbh->selectall_arrayref(
    <<"SQL"
    SELECT name, supplier_name
    FROM   current_samples
    WHERE  (name LIKE '%STDY%' OR supplier_name LIKE '%STDY%')
    AND    supplier_name IS NOT NULL
SQL
);
foreach ( @{$ary_ref} ) {
    my ( $name, $supplier_name ) = @{$_};
    $name_for{$name} = $supplier_name;
}

# If only one count file then assume it's a list of count files
if ( scalar @count_files == 1 ) {
    my $count_files = path( $count_files[0] )->slurp;
    chomp $count_files;
    @count_files = split /[\r\n]+/xms, $count_files;
}

# Get filehandles for all count files
my %fh_for;
my @samples;
foreach my $count_file (@count_files) {
    my ( undef, undef, $filename ) = File::Spec->splitpath($count_file);
    my ($sample) = $filename =~ m/\A (.*) [.] \w+ \z/xms;
    confess sprintf q{Can't parse sample name from %s}, $filename if !$sample;
    if ( $sample =~ m/\A \d+STDY\d+ \z/xms ) {
        $sample = $name_for{$sample};
    }
    push @samples, $sample;
    open my $fh, '<', $count_file;    ## no critic (RequireBriefOpen)
    $fh_for{$sample} = $fh;
}

# Output counts file
open my $counts_fh, '>', $COUNTS_FILE;    ## no critic (RequireBriefOpen)
printf {$counts_fh} "\t%s\n", join "\t", @samples;
while (1) {
    my $current_gene;
    my @counts;
    foreach my $sample (@samples) {
        my $line = readline $fh_for{$sample};
        chomp $line;
        my ( $gene, $count ) = split /\t/xms, $line;
        if ( defined $current_gene && $current_gene ne $gene ) {
            confess sprintf 'Gene order inconsistent between files (%s vs %s)',
              $current_gene, $gene;
        }
        $current_gene = $gene;
        push @counts, $count;
    }
    last if $current_gene =~ m/\A _/xms || $current_gene =~ m/\A no_feature/xms;
    printf {$counts_fh} "%s\t%s\n", $current_gene, join "\t", @counts;
}
close $counts_fh;

# Output samples file
open my $samples_fh, '>', $SAMPLES_FILE;
print {$samples_fh} "\tcondition\n";    ## no critic (RequireCheckedSyscalls)
foreach my $sample (@samples) {
    my $condition = $sample;
    $condition =~ s/_?\d+ \z//xms;
    printf {$samples_fh} "%s\t%s\n", $sample, $condition;
}
close $samples_fh;

# Close all filehandles
foreach my $sample (@samples) {
    close $fh_for{$sample};
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'count_files=s@{1,}' => \@count_files,
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

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

make_deseq_from_htseq.pl

Make DESeq2 counts and samples files from htseq-count files

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a list of htseq-count output files and creates merged
counts.txt and samples.txt files for DESeq2 in the current directory.
Additionally, Sequencescape sample names (e.g. 3407STDY6401339) are replaced
with more readable sample names (e.g. zmp_ph250_E4_p2).

=head1 EXAMPLES

    perl \
        make_deseq_from_htseq.pl \
            --count_files 21somites_1.count 21somites_2.count

    perl \
        make_deseq_from_htseq.pl --count_files somites.fofn

=head1 USAGE

    make_deseq_from_htseq.pl
        [--count_files file...]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--count_files FILE...>

htseq-count output files. If just one file is specified then it's assumed to be
a file of filenames.

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

This software is Copyright (c) 2016 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
