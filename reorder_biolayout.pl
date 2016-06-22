#!/usr/bin/env perl

# PODNAME: reorder_biolayout.pl
# ABSTRACT: Reorder columns of BioLayout Express3D expression file

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-06-22

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;

# Constants
Readonly our $NUM_CORE_FIELDS => 12;
Readonly our $LINES_TO_CHECK  => 100;

# Default options
my @sort_by;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Get rows to sort by
my %unseen_sort_by = map { $_ => 1 } @sort_by;
my %sort_by_fields;
my @stored_lines;
my $line_count;
while ( my $line = <> ) {
    $line_count++;
    push @stored_lines, $line;
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $id = $fields[0];
    if ( exists $unseen_sort_by{$id} ) {
        splice @fields, 0, $NUM_CORE_FIELDS;
        $sort_by_fields{$id} = \@fields;
        delete $unseen_sort_by{$id};
    }
    last if !scalar keys %unseen_sort_by;
    if ( $line_count >= $LINES_TO_CHECK ) {
        confess sprintf 'Columns to sort by (%s) not seen in first %s lines',
          ( join q{, }, keys %unseen_sort_by ), $LINES_TO_CHECK;
    }
}

# Reorder columns
my $num_columns = scalar @{ $sort_by_fields{ $sort_by[0] } };
my @columns = sort { cmp_fields( $a, $b, \@sort_by, \%sort_by_fields ) }
  ( 0 .. $num_columns - 1 );

# Output reordered columns
foreach my $line (@stored_lines) {
    print_reordered_line( $line, \@columns );
}
while ( my $line = <> ) {
    print_reordered_line( $line, \@columns );
}

sub cmp_fields {
    my ( $a, $b, $sort_by, $sort_by_fields ) = @_;

    foreach my $row ( @{$sort_by} ) {
        if ( $sort_by_fields->{$row}->[$a] lt $sort_by_fields->{$row}->[$b] ) {
            return -1;    ## no critic (ProhibitMagicNumbers)
        }
        elsif ( $sort_by_fields->{$row}->[$a] gt $sort_by_fields->{$row}->[$b] )
        {
            return 1;
        }
    }

    return 0;
}

sub print_reordered_line {
    my ( $line, $columns ) = @_;

    chomp $line;
    my @fields = split /\t/xms, $line;
    my @core = splice @fields, 0, $NUM_CORE_FIELDS;
    @fields = @fields[@columns];
    printf "%s\n", join "\t", @core, @fields;

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'sort_by=s@{1,}' => \@sort_by,
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

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

reorder_biolayout.pl

Reorder columns of BioLayout Express3D expression file

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a BioLayout Express3D expression file and reorders the columns
by specified rows.

=head1 EXAMPLES

    perl \
        reorder_biolayout.pl \
        --sort_by ZFS_ID Sample < all.expression > all.reorder.expression

=head1 USAGE

    reorder_biolayout.pl
        [--sort_by rows]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--sort_by ROWS>

Row names to sort by.

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
