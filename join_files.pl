#!/usr/bin/env perl

# PODNAME: join_files.pl
# ABSTRACT: Join two files by common field

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2015-11-19

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

# Default options
my $file1;
my $file2;
my $file1_key  = 1;
my $file2_key  = 1;
my $has_header = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

## no critic (RequireBriefOpen)
open my $fh1, q{<}, $file1;
open my $fh2, q{<}, $file2;
## use critic

# Print header line if present
if ($has_header) {
    my $line1 = <$fh1>;
    my $line2 = <$fh2>;
    chomp $line1;
    chomp $line2;

    my @fields1 = split /\t/xms, $line1;
    my @fields2 = split /\t/xms, $line2;
    if ($debug) {
        ## no critic (RequireCarping)
        warn sprintf "Key field for %s: %s\n", $file1,
          $fields1[ $file1_key - 1 ];
        warn sprintf "Key field for %s: %s\n", $file2,
          $fields2[ $file2_key - 1 ];
        ## use critic
    }

    splice @fields2, $file2_key - 1, 1;
    printf "%s\t%s\n", $line1, join "\t", @fields2;
}

# Read all of second file into memory
my %lines_for;
while ( my $line = <$fh2> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $key = splice @fields, $file2_key - 1, 1;
    push @{ $lines_for{$key} }, join "\t", @fields;
}

# Iterate over first file and output joined lines
while ( my $line = <$fh1> ) {
    chomp $line;
    my @fields = split /\t/xms, $line;
    my $key = $fields[ $file1_key - 1 ];
    foreach my $rest_of_line ( @{ $lines_for{$key} } ) {
        printf "%s\t%s\n", $line, $rest_of_line;
    }
}

close $fh1;
close $fh2;

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'file1=s'     => \$file1,
        'file2=s'     => \$file2,
        'file1_key=i' => \$file1_key,
        'file2_key=i' => \$file2_key,
        'header'      => \$has_header,
        'debug'       => \$debug,
        'help'        => \$help,
        'man'         => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$file1 ) {
        pod2usage("--file1 must be specified\n");
    }
    if ( !$file2 ) {
        pod2usage("--file2 must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

join_files.pl

Join two tab-delimited files by common field

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes two tab-delimited files and joins them, inefficiently, by a
common field. The order of the first file is preserved.

=head1 EXAMPLES

    perl \
        join_files.pl \
        --file1 in1.tsv --file2 in2.tsv \
        > output.tsv

    perl \
        join_files.pl \
        --file1 in1.tsv --file2 in2.tsv \
        --file1_key 2 --file2_key 3 --header \
        > output.tsv

=head1 USAGE

    join_files.pl
        [--file1 file]
        [--file2 file]
        [--file1_key field]
        [--file2_key field]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--file1 FILE>

The first file.

=item B<--file2 FILE>

The second file.

=item B<--file1_key FIELD>

The field in the first file to join by.

=item B<--file2_key FIELD>

The field in the second file to join by.

=item B<--header>

Input has a header line.

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

This software is Copyright (c) 2015 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
