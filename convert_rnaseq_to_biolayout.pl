#!/usr/bin/env perl

# PODNAME: convert_rnaseq_to_biolayout.pl
# ABSTRACT: Convert RNA-Seq output into BioLayout Express3D format

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-04-25

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Readonly;

# Constants
Readonly our $NAME_FIELD  => 9;
Readonly our @CORE_FIELDS => ( 0 .. 10 );

# Default options
my $input_file;
my $samples_file;
my @metadata_files;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

my @sample_cols = output_header( $input_file, $samples_file, @metadata_files );
output_regions( $input_file, @sample_cols );

# Output header
sub output_header {
    ## no critic (ProhibitReusedNames)
    my ( $input_file, $samples_file, @metadata_files ) = @_;
    ## use critic

    # Get input headings
    open my $input_fh, '<', $input_file;
    my $header = <$input_fh>;
    chomp $header;
    my @headings = split /\t/xms, $header;
    close $input_fh;

    my @sample_cols;    ## no critic (ProhibitReusedNames)
    my $i = -1;         ## no critic (ProhibitMagicNumbers)
    foreach my $heading (@headings) {
        $i++;
        if ( $heading =~ m/\s normalised \s count \z/xms ) {
            push @sample_cols, $i;
        }
    }

    my @output_headings = ('ID');
    push @output_headings, @headings[@CORE_FIELDS];
    push @output_headings, @headings[@sample_cols];

    my @tmp_headings;
    foreach my $heading (@output_headings) {
        $heading =~ s/\s normalised \s count \z//xms;
        push @tmp_headings, $heading;
    }
    @output_headings = @tmp_headings;
    @output_headings = map { /\s/xms ? qq{"$_"} : $_ } @output_headings;

    printf "%s\n", join "\t", @output_headings;

    if ($samples_file) {
        output_samples_header($samples_file);
    }

    return @sample_cols;
}

# Output samples header
sub output_samples_header {
    my ($samples_file) = @_;    ## no critic (ProhibitReusedNames)

    open my $samples_fh, '<', $samples_file;
    my $header = <$samples_fh>;
    $header =~ s/\A \s+//xms;
    my @columns = split /\s+/xms, $header;
    close $samples_fh;

    # Get all sample headings
    my @all_sample_headings;
    foreach my $col ( 1 .. scalar @columns ) {
        my @sample_headings =
          ( $columns[ $col - 1 ], (q{}) x scalar @CORE_FIELDS );
        open my $samples_fh, '<', $samples_file;
        $header = <$samples_fh>;
        while ( my $line = <$samples_fh> ) {
            chomp $line;
            my @fields = split /\s+/xms, $line;
            push @sample_headings, $fields[$col];
        }
        close $samples_fh;
        @sample_headings = map { /\s/xms ? qq{"$_"} : $_ } @sample_headings;
        push @all_sample_headings, \@sample_headings;
    }

    # Concatenate all factors if more than one
    if ( scalar @all_sample_headings > 1 ) {
        my @concat_headings = ( join q{_}, @columns );
        push @concat_headings, (q{}) x scalar @CORE_FIELDS;
        foreach my $sample_idx (
            scalar @CORE_FIELDS + 1 .. scalar @{ $all_sample_headings[0] } - 1 )
        {
            my @values;
            foreach my $factor_idx ( 0 .. scalar @columns - 1 ) {
                push @values, $all_sample_headings[$factor_idx][$sample_idx];
            }
            push @concat_headings, join q{_}, @values;
        }
        unshift @all_sample_headings, \@concat_headings;
    }

    foreach my $headings (@all_sample_headings) {
        printf "%s\n", join "\t", @{$headings};
    }

    return;
}

# Output regions of interest
sub output_regions {
    my ( $file, @sample_cols ) = @_;    ## no critic (ProhibitReusedNames)

    # First pass to get duplicated gene names
    my %is_duplicate;
    open my $fh, '<', $file;
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = split /\t/xms, $line;
        $is_duplicate{ $fields[$NAME_FIELD] }++;
    }
    close $fh;
    %is_duplicate =
      map { $_ => 1 } grep { $is_duplicate{$_} > 1 } keys %is_duplicate;

    open $fh, '<', $file;    ## no critic (RequireBriefOpen)
    $header = <$fh>;
    while ( my $line = <$fh> ) {
        chomp $line;
        my @fields = split /\t/xms, $line;

        # Get ID by joining gene name and Ensembl ID if duplicated
        my $id =
          exists $is_duplicate{ $fields[$NAME_FIELD] }
          ? join q{:}, @fields[ $NAME_FIELD, 0 ]
          : $fields[$NAME_FIELD];

        my @output_fields = ($id);
        push @output_fields, @fields[@CORE_FIELDS];
        push @output_fields, @fields[@sample_cols];

        @output_fields = map { defined $_ ? $_       : q{-} } @output_fields;
        @output_fields = map { /\s/xms    ? qq{"$_"} : $_ } @output_fields;
        @output_fields =
          map {
            /\A ([-]?[\d.]+)(e[-]?\d+) \z/xms
              ? ( sprintf '%.8f', $1 ) . uc $2
              : $_
          } @output_fields;
        @output_fields =
          map { /\A [-]?\d+[.]\d+ \z/xms ? sprintf '%.8f', $_ : $_ }
          @output_fields;

        printf "%s\n", join "\t", @output_fields;
    }
    close $fh;

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_file=s'         => \$input_file,
        'samples_file=s'       => \$samples_file,
        'metadata_files=s@{,}' => \@metadata_files,
        'debug'                => \$debug,
        'help'                 => \$help,
        'man'                  => \$man,
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

convert_rnaseq_to_biolayout.pl

Convert RNA-Seq output into BioLayout Express3D format

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a tab-delimited file of RNA-Seq output produced by
https://gist.github.com/iansealy/2dca28d07c0764e014df or
https://gist.github.com/iansealy/b9cbc56bd1affe10d37a and converts to the format
required by BioLayout Express3D.

=head1 EXAMPLES

    perl \
        convert_rnaseq_to_biolayout.pl \
        --input_file all.tsv > all.expression

    perl \
        convert_rnaseq_to_biolayout.pl \
        --input_file all.tsv --samples_file samples.txt > all.expression

    perl \
        convert_rnaseq_to_biolayout.pl \
        --input_file all.tsv --samples_file samples.txt \
        --metadata_files \
            QC_all_stages_lab.tsv \
            QC_all_stages_npg.tsv \
            QC_all_stages_post_DETCT.tsv \
        > all.expression

=head1 USAGE

    convert_to_biolayout.pl
        [--input_file file]
        [--samples_file file]
        [--metadata_files files]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_file FILE>

RNA-Seq output file (e.g. all.tsv).

=item B<--samples_file FILE>

DESeq2 samples file (e.g. samples.txt).

=item B<--metadata_files FILES>

QC metadata files.

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
