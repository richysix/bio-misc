#!/usr/bin/env perl

# PODNAME: gather_rnaseq_count_qc_data.pl
# ABSTRACT: Gather count data for RNA-Seq experiments for QC

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-05-05

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;
use DBI;
use Sort::Naturally;

# Constants
Readonly our $HOST => 'seqw-db';
Readonly our $PORT => '3379';
Readonly our $USER => 'warehouse_ro';
Readonly our $PASS => q{};
Readonly our $NAME => 'sequencescape_warehouse';

Readonly our $SPIKE_PREFIX => 'ERCC';

Readonly our $OUTPUT_NAMES => join "\t",
  qw( sample genomic_region_counts spike_counts );

# Default options
my $input_file;
my @expts;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

my ( $genomic_count, $spike_count ) = get_counts($input_file);

# Connect to database
my $dsn = "dbi:mysql:database=$NAME;host=$HOST;port=$PORT";
my $dbh = DBI->connect( $dsn, $USER, $PASS );

printf "%s\n", $OUTPUT_NAMES;

foreach my $expt (@expts) {
    my $expt_q = $dbh->quote($expt);
    $expt_q =~ s/'\z/%'/xms;
    my $ary_ref = $dbh->selectall_arrayref(
        <<"SQL"
        SELECT name, supplier_name, public_name
        FROM   current_samples
        WHERE  (name LIKE $expt_q
                OR supplier_name LIKE $expt_q
                OR public_name LIKE $expt_q)
        AND    description NOT LIKE '%polyT%'
SQL
    );
    my @samples;
    foreach ( @{$ary_ref} ) {
        my ( $name, $supplier_name, $public_name ) = @{$_};
        if ( $name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $name;
        }
        elsif ( $supplier_name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $supplier_name;
        }
        elsif ( $public_name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $public_name;
        }
    }

    foreach my $sample ( sort { ncmp( $a, $b ) } @samples ) {
        my @values = ($sample);
        push @values, $genomic_count->{$sample};
        push @values, $spike_count->{$sample};
        @values = map { !defined $_ ? q{} : $_ } @values;
        printf "%s\n", join "\t", @values;
    }
}

if ( !@expts ) {
    foreach my $sample ( sort { ncmp( $a, $b ) } keys %{$genomic_count} ) {
        my @values = ($sample);
        push @values, $genomic_count->{$sample};
        push @values, $spike_count->{$sample};
        @values = map { !defined $_ ? q{} : $_ } @values;
        printf "%s\n", join "\t", @values;
    }
}

# Get counts from input file
sub get_counts {
    my ($file) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }
    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)

    # Get input headings
    my $header = <$fh>;
    my @headings = split /\t/xms, $header;

    # Get columns
    my @sample_cols;
    my %col_to_sample;
    my $col = -1;               ## no critic (ProhibitMagicNumbers)
    foreach my $heading (@headings) {
        $col++;
        if ( $heading =~ m/\A (\S+) \s count \z/xms ) {
            $col_to_sample{$col} = $1;
            push @sample_cols, $col;
        }
    }

    # Get counts
    my %genomic_count = map { $_ => 0 } values %col_to_sample;
    my %spike_count   = map { $_ => 0 } values %col_to_sample;
    while ( my $line = <$fh> ) {
        chomp $line;
        my @values = split /\t/xms, $line;
        foreach my $col (@sample_cols) {
            ## no critic (ProhibitMagicNumbers)
            if ( $values[4] =~ m/\A $SPIKE_PREFIX/xms ) {
                ## use critic
                $spike_count{ $col_to_sample{$col} } += $values[$col];
            }
            else {
                $genomic_count{ $col_to_sample{$col} } += $values[$col];
            }
        }
    }

    close $fh;

    return \%genomic_count, \%spike_count;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_file=s' => \$input_file,
        'expts=s@{1,}' => \@expts,
        'debug'        => \$debug,
        'help'         => \$help,
        'man'          => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$input_file ) {
        pod2usage("--input_file must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

gather_rnaseq_count_qc_data.pl

Gather count data for RNA-Seq experiments for QC

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file of RNA-Seq output produced by
https://gist.github.com/iansealy/2dca28d07c0764e014df or
https://gist.github.com/iansealy/b9cbc56bd1affe10d37a and an optional list of
experiment prefixes (e.g. zmp_ph80) and gathers count data for QC.

=head1 EXAMPLES

    perl \
        gather_rnaseq_count_qc_data.pl \
        --input_file all.tsv --expts zmp_ph80 zmp_ph81 > count_qc.tsv

    perl \
        gather_rnaseq_count_qc_data.pl \
        --input_file all.tsv > count_qc.tsv

=head1 USAGE

    gather_rnaseq_count_qc_data.pl
        [--input_file file]
        [--expts prefixes]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_file FILE>

RNA-Seq output file (e.g. all.tsv).

=item B<--expts PREFIXES>

Experiment prefixes (e.g. zmp_ph80).

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
