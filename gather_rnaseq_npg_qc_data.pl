#!/usr/bin/env perl

# PODNAME: gather_rnaseq_npg_qc_data.pl
# ABSTRACT: Gather NPG data for RNA-Seq experiments for QC

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

Readonly our @COMMON_NAMES =>
  qw( tag_index tag_sequence insert_size_quartile1 insert_size_median insert_size_quartile3 gc_percent_forward_read gc_percent_reverse_read sequence_mismatch_percent_forward_read sequence_mismatch_percent_reverse_read adapters_percent_forward_read adapters_percent_reverse_read tag_decode_percent tag_decode_count bam_num_reads bam_percent_mapped bam_percent_duplicate );
Readonly our $OUTPUT_NAMES => join q{, },
  ( qw( sample run lane ), @COMMON_NAMES );
Readonly our $FIELD_NAMES => join q{, },
  ( qw( id_run position ), @COMMON_NAMES );

# Default options
my @expts;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connect to database
my $dsn = "dbi:mysql:database=$NAME;host=$HOST;port=$PORT";
my $dbh = DBI->connect( $dsn, $USER, $PASS );

printf "%s\n", $OUTPUT_NAMES;

foreach my $expt (@expts) {
    my $expt_q = $dbh->quote($expt);
    $expt_q =~ s/'\z/%'/xms;
    my $ary_ref = $dbh->selectall_arrayref(
        <<"SQL"
        SELECT name, supplier_name
        FROM   current_samples
        WHERE  (name LIKE $expt_q OR supplier_name LIKE $expt_q)
        AND    description NOT LIKE '%polyT%'
SQL
    );
    my @samples;
    foreach ( @{$ary_ref} ) {
        my ( $name, $supplier_name ) = @{$_};
        if ( $name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $name;
        }
        elsif ( $supplier_name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $supplier_name;
        }
    }

    foreach my $sample ( sort { ncmp( $a, $b ) } @samples ) {
        my $sample_q = $dbh->quote($sample);
        $ary_ref = $dbh->selectall_arrayref(
            <<"SQL"
            SELECT   $FIELD_NAMES
            FROM     npg_plex_information npi, current_samples cs
            WHERE    npi.sample_id = cs.internal_id
            AND      (supplier_name = $sample_q OR name = $sample_q)
            ORDER BY id_run, position
SQL
        );
        foreach ( @{$ary_ref} ) {
            my @values = map { !defined $_ ? q{} : $_ } @{$_};
            printf "%s\n", join "\t", $sample, @values;
        }
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
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
    if ( !@expts ) {
        pod2usage("--expts must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

gather_rnaseq_npg_qc_data.pl

Gather NPG data for RNA-Seq experiments for QC

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a list of experiment prefixes (e.g. zmp_ph80) and gathers NPG
data for QC.

For samples run on multiple lanes, merge_duplicates_by_type.pl will need to also
be run.

=head1 EXAMPLES

    perl \
        gather_rnaseq_npg_qc_data.pl \
        --expts zmp_ph80 zmp_ph81 > npg_qc.tsv

    perl \
        merge_duplicates_by_type.pl \
        < npg_qc.tsv \
        --header \
        --key 1 \
        --mean 6 7 8 9 10 11 12 13 14 15 18 19 \
        --sum 16 17 \
        > npg_qc_uniq.tsv

=head1 USAGE

    gather_rnaseq_npg_qc_data.pl
        [--expts prefixes]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

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
