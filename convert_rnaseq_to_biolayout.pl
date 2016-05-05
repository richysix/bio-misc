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
my $readable_samples;
my $config_file;
my @metadata_files;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

my @sample_cols = output_header(
    $input_file,  $samples_file, $readable_samples,
    $config_file, @metadata_files
);
output_regions( $input_file, @sample_cols );

# Output header
sub output_header {
    ## no critic (ProhibitReusedNames)
    my (
        $input_file,  $samples_file, $readable_samples,
        $config_file, @metadata_files
    ) = @_;
    ## use critic

    # Get input headings
    open my $input_fh, '<', $input_file;
    my $header = <$input_fh>;
    chomp $header;
    my @headings = split /\t/xms, $header;
    close $input_fh;

    my @sample_cols;    ## no critic (ProhibitReusedNames)
    my %sample_to_col;
    my $col = 0;
    my $i   = -1;       ## no critic (ProhibitMagicNumbers)
    foreach my $heading (@headings) {
        $i++;
        if ( $heading =~ m/\s normalised \s count \z/xms ) {
            push @sample_cols, $i;
            $heading =~ s/\s normalised \s count \z//xms;
            $sample_to_col{$heading} = ++$col;
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

    my @all_output_headings = ( \@output_headings );

    if ($samples_file) {
        push @all_output_headings, output_samples_header($samples_file);
    }
    if (@metadata_files) {
        push @all_output_headings,
          output_metadata_header( $config_file, \%sample_to_col,
            @metadata_files );
    }

    if ($readable_samples) {

        # Swap first two sample headings
        ( $all_output_headings[0], $all_output_headings[1] ) =
          ( $all_output_headings[1], $all_output_headings[0] );
        @tmp_headings =
          @{ $all_output_headings[0] }[ 0 .. scalar @CORE_FIELDS ];
        @{ $all_output_headings[0] }[ 0 .. scalar @CORE_FIELDS ] =
          @{ $all_output_headings[1] }[ 0 .. scalar @CORE_FIELDS ];
        @{ $all_output_headings[1] }[ 0 .. scalar @CORE_FIELDS ] =
          @tmp_headings;
    }

    foreach my $output_headings (@all_output_headings) {
        printf "%s\n", join "\t", @{$output_headings};
    }

    return @sample_cols;
}

# Output samples header
sub output_samples_header {
    my ($samples_file) = @_;    ## no critic (ProhibitReusedNames)

    open my $samples_fh, '<', $samples_file;    ## no critic (RequireBriefOpen)
    my $header = <$samples_fh>;
    chomp $header;
    if ( $header =~ m/\A \s/xms ) {
        $header =~ s/\A \s+//xms;
    }
    else {
        $header =~ s/\A \S+ \s+//xms;
    }
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

    # Convert to unique readable sample names
    my @headings = @{ $all_sample_headings[0] };
    my @readable = ('Sample');
    shift @headings;
    push @readable, (q{}) x scalar @CORE_FIELDS;
    my %heading_ordinal;
    foreach my $heading (@headings) {
        next if !$heading;
        push @readable, $heading . q{_} . ++$heading_ordinal{$heading};
    }
    unshift @all_sample_headings, \@readable;

    return @all_sample_headings;
}

# Output metadata header
sub output_metadata_header {
    ## no critic (ProhibitReusedNames)
    my ( $config_file, $sample_to_col, @metadata_files ) = @_;
    ## use critic

    my %round;
    my %skip;
    if ($config_file) {
        open my $fh, '<', $config_file;    ## no critic (RequireBriefOpen)
        while ( my $line = <$fh> ) {
            chomp $line;
            my ( $type, $value ) = split /\s+/xms, $line;
            next if !$value;
            if ( $value eq q{X} ) {
                $skip{$type} = 1;
            }
            elsif ( $value =~ m/\A [\d.]+ \z/xms ) {
                $round{$type} = $value;
            }
        }
        close $fh;
    }

    my @all_metadata_headings;
    foreach my $file (@metadata_files) {
        open my $fh, '<', $file;
        my $header = <$fh>;
        chomp $header;
        my @columns = split /\t/xms, $header;
        shift @columns;    # Ignore sample name column
        close $fh;

        # Get all metadata headings
        foreach my $col ( 1 .. scalar @columns ) {
            next if exists $skip{ $columns[ $col - 1 ] };
            my @headings =
              ( $columns[ $col - 1 ], (q{}) x scalar @CORE_FIELDS );
            open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
            $header = <$fh>;
            while ( my $line = <$fh> ) {
                chomp $line;
                my @fields = split /\t/xms, $line;
                my $sample = $fields[0];
                confess sprintf "Sample %s unknown\n", $sample
                  if !exists $sample_to_col->{$sample};
                if ( exists $round{ $columns[ $col - 1 ] } ) {
                    my $round_to = $round{ $columns[ $col - 1 ] };
                    $fields[$col] =
                      int( ( $fields[$col] + $round_to / 2 ) / $round_to ) *
                      $round_to;
                }
                $headings[ $sample_to_col->{$sample} + scalar @CORE_FIELDS ] =
                  $fields[$col];
            }
            close $fh;
            @headings = map { /\s/xms ? qq{"$_"} : $_ } @headings;
            push @all_metadata_headings, \@headings;
        }
    }

    return @all_metadata_headings;
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
        'readable_samples'     => \$readable_samples,
        'config_file=s'        => \$config_file,
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
    if ( $readable_samples && !$samples_file ) {
        pod2usage("--readable_samples specified without --samples_file\n");
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
        --input_file all.tsv --samples_file samples.txt --readable_samples \
        --config_file fields.txt \
        --metadata_files \
            QC_all_stages_lab.tsv \
            QC_all_stages_npg.tsv \
            QC_all_stages_post_DETCT.tsv \
        > all.expression

=head1 USAGE

    convert_to_biolayout.pl
        [--input_file file]
        [--samples_file file]
        [--readable_samples]
        [--config_file file]
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

=item B<--readable_samples>

Force sample names to be readable.

=item B<--config_file FILE>

Metadata config file. One data type per row with name in first column and either
an integer (specifying how to round values) or X (indicating types to ignore) in
second column.

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
