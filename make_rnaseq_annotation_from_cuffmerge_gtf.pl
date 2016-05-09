#!/usr/bin/env perl

# PODNAME: make_rnaseq_annotation_from_cuffmerge_gtf.pl
# ABSTRACT: Make annotation for RNA-Seq from a Cuffmerge GTF file

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-05-09

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over STDIN and merge exons into genes
my %annotation_for;
while ( my $line = <> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    ## no critic (ProhibitMagicNumbers)
    my ( $chr, $start, $end, $strand ) = @fields[ 0, 3, 4, 6 ];
    ## use critic
    if ( $strand eq q{+} ) {
        $strand = 1;
    }
    elsif ( $strand eq q{-} ) {
        $strand = -1;    ## no critic (ProhibitMagicNumbers)
    }
    ## no critic (ProhibitMagicNumbers)
    my ($id) = $fields[8] =~ m/gene_id \s "([^"]+)"/xms;
    ## use critic
    if ( !exists $annotation_for{$id} ) {
        $annotation_for{$id} = [ $chr, $start, $end, $strand ];
    }
    else {
        if ( $start < $annotation_for{$id}->[1] ) {
            $annotation_for{$id}->[1] = $start;
        }
        if ( $end > $annotation_for{$id}->[2] ) {
            $annotation_for{$id}->[2] = $end;
        }
    }
}

foreach my $id ( sort keys %annotation_for ) {
    printf "%s\n", join "\t", $id, @{ $annotation_for{$id} }, 'unknown', $id,
      q{};
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'debug' => \$debug,
        'help'  => \$help,
        'man'   => \$man,
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

make_rnaseq_annotation_from_cuffmerge_gtf.pl

Make annotation for RNA-Seq from a Cuffmerge GTF file

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a Cuffmerge GTF file and produces annotation suitable for use
with RNA-Seq output produced by
https://gist.github.com/iansealy/a51daa2b30d19f062c1882b2b804c97a

=head1 EXAMPLES

    perl \
        make_rnaseq_annotation_from_cuffmerge_gtf.pl \
        < merged.gtf > annotation.txt

=head1 USAGE

    make_rnaseq_annotation_from_cuffmerge_gtf.pl
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

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
