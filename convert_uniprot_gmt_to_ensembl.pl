#!/usr/bin/env perl

# PODNAME: convert_uniprot_gmt_to_ensembl.pl
# ABSTRACT: Convert UniProt IDs in GSEA GMT file to Ensembl IDs

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-11-13

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my $mapping_file;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get mapping
my %mapping;
open my $mapping_fh, '<', $mapping_file;
my $header = <$mapping_fh>;
while ( my $line = <$mapping_fh> ) {
    chomp $line;
    my ( $ensembl, undef, undef, $uniprot ) = split /\t/xms, $line;
    $mapping{$uniprot}{$ensembl} = 1;
}
close $mapping_fh;

# Iterate over STDIN
while ( my $line = <> ) {
    chomp $line;
    my ( $name, $description, @uniprot_ids ) = split /\t/xms, $line;
    my %ensembl_id;
    foreach my $uniprot_id (@uniprot_ids) {
        if ( !exists $mapping{$uniprot_id} ) {
            carp sprintf "No mapping for %s (%s / %s)\n", $uniprot_id, $name,
              $description;
        }
        else {
            my @ensembl = sort keys %{ $mapping{$uniprot_id} };
            if ( scalar @ensembl > 1 ) {
                carp sprintf "Multiple mappings (%s) for %s (%s / %s)\n",
                  ( join q{,}, @ensembl ), $uniprot_id, $name, $description;
            }
            foreach my $ensembl (@ensembl) {
                $ensembl_id{$ensembl} = 1;
            }
        }
    }
    if ( !%ensembl_id ) {
        carp sprintf "No Ensembl IDs for %s / %s\n", $name, $description;
    }
    else {
        printf "%s\t%s\t%s\n", $name, $description,
          ( join "\t", sort keys %ensembl_id );
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'mapping_file=s' => \$mapping_file,
        'debug'          => \$debug,
        'help'           => \$help,
        'man'            => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$mapping_file ) {
        pod2usage("--mapping_file must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

convert_uniprot_gmt_to_ensembl.pl

Convert UniProt IDs in GSEA GMT file to Ensembl IDs

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a GSEA GMT file on STDIN and convers UniProt IDs to Ensembl
IDs.

=head1 EXAMPLES

    perl \
        convert_uniprot_gmt_to_ensembl.pl \
        --mapping_file Mus_musculus.GRCm38.88.uniprot.tsv \
        < input.gmt > output.gmt

=head1 USAGE

    convert_uniprot_gmt_to_ensembl.pl
        [--mapping_file file]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--mapping_file FILE>

File downloaded from Ensembl that maps Ensembl IDs to UniProt IDs. For example:
ftp://ftp.ensembl.org/pub/release-88/tsv/mus_musculus/Mus_musculus.GRCm38.88.uniprot.tsv.gz

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
