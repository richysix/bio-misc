#!/usr/bin/env perl

# PODNAME: get_ensemblgenomes_gene_annotation.pl
# ABSTRACT: Get annotation for each Ensembl Genomes gene

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-10-05

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Bio::EnsEMBL::LookUp;

# Default options
my $species = 'escherichia_coli_str_k_12_substr_dh10b';
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get Ensembl adaptor
my $lookup = Bio::EnsEMBL::LookUp->new();
my $dba    = $lookup->get_by_name_exact($species);
my $ga     = $dba->get_GeneAdaptor();

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Get all genes
my $genes = $ga->fetch_all();
foreach my $gene ( sort { $a->stable_id cmp $b->stable_id } @{$genes} ) {
    printf "%s\n", join "\t", $gene->stable_id, $gene->seq_region_name,
      $gene->seq_region_start, $gene->seq_region_end, $gene->seq_region_strand,
      $gene->biotype, ( $gene->external_name || $gene->stable_id ),
      ( $gene->description || q{} );
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'species=s' => \$species,
        'debug'     => \$debug,
        'help'      => \$help,
        'man'       => \$man,
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

get_ensemblgenomes_gene_annotation.pl

Get annotation for each Ensembl Genomes gene

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script dumps a list of Ensembl Genomes gene stable IDs along with some
annotation.

=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-90/ensembl/modules \
        -Ibranch-ensemblgenomes-37/ensemblgenomes-api/modules \
        get_ensemblgenomes_gene_annotation.pl

=head1 USAGE

    get_ensemblgenomes_gene_annotation.pl
        [--species species]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--species SPECIES>

Species (defaults to escherichia_coli_str_k_12_substr_dh10b).

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=head1 DEPENDENCIES

Ensembl Perl API - http://www.ensembl.org/info/docs/api/

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
