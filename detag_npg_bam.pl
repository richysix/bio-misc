#!/usr/bin/env perl

# PODNAME: detag_npg_bam.pl
# ABSTRACT: Detag NPG BAM files

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-09-11

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use Bio::DB::HTS;
use Text::Fuzzy;

# Default options
my $tag_file;
my @bam_files;
my $edit_distance = 2;
my $output_prefix = q{};
## no critic (ProhibitMagicNumbers)
my $random_start = 1;
my $random_end   = 10;
my $inread_start = 11;
my $inread_end   = 18;
## use critic
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Change to zero-based numbering
$random_start--;
$inread_start--;
my $random_length = $random_end - $random_start;
my $inread_length = $inread_end - $inread_start;

# Read tags and create filehandles for all output files
my %fh_for;
open $fh_for{q{}}{1},  '>', $output_prefix . 'notag.1.fastq';
open $fh_for{q{}}{2},  '>', $output_prefix . 'notag.2.fastq';
open $fh_for{q{ }}{1}, '>', $output_prefix . 'multitag.1.fastq';
open $fh_for{q{ }}{2}, '>', $output_prefix . 'multitag.2.fastq';
open my $tag_fh, '<', $tag_file;
while ( my $line = <$tag_fh> ) {
    chomp $line;
    my ( $sample, $tag ) = split /\s+/xms, $line;
    confess sprintf 'Tag (%s) seen twice in tag file', $tag
      if exists $fh_for{$tag};
    open $fh_for{$tag}{1}, '>', $output_prefix . $sample . '.1.fastq';
    open $fh_for{$tag}{2}, '>', $output_prefix . $sample . '.2.fastq';
}
close $tag_fh;
my @tags = sort keys %fh_for;
my @tag_lengths = sort map { length } grep { m/\w/xms } @tags;
confess 'Tag lengths vary' if $tag_lengths[0] != $tag_lengths[-1];
my $tag_length = $tag_lengths[-1];

# Read BAM files and write FASTQ files
my %unknown;
my $multi = 0;
foreach my $bam_file (@bam_files) {
    my $hts = Bio::DB::HTS->new( -bam => $bam_file );
    confess 'BAM file must be sorted by name'
      if $hts->header->text =~ m/SO:coordinate/xms;
    my $iterator = $hts->features(
        ## no critic (RequireInterpolationOfMetachars)
        -filter => q{return if $a->get_tag_values('SUPPLEMENTARY')},
        ## use critic
        -iterator => 1,
    );
    while ( my $read1 = $iterator->next_seq ) {
        my $read2 = $iterator->next_seq;
        if ( $read2->get_tag_values('FIRST_MATE') ) {
            ( $read1, $read2 ) = ( $read2, $read1 );
        }
        my $nearest_tag;
        my $barcode   = $read1->aux_get('BC');
        my $random    = $read1->aux_get('br');
        my $read1_seq = $read1->query->dna;
        confess sprintf
          'Got short barcode (%s) and random bases (%s) for read (%s)',
          $barcode, $random, $read1->qname
          if $random && length $barcode < $tag_length;
        if ( length $barcode < $tag_length ) {
            $random = substr $read1_seq, $random_start, $random_length,
              q{-} x $random_length;
            $barcode .= substr $read1_seq, $inread_start, $inread_length,
              q{-} x $inread_length;
            $read1_seq =~ s/\-//xmsg;
        }
        my $description  = q{};
        my $tf           = Text::Fuzzy->new( $barcode, max => $edit_distance );
        my @nearest_tags = $tf->nearestv( \@tags );
        if ( !@nearest_tags ) {
            $nearest_tag = q{};               # No tag
            $description = q{ } . $barcode;
            $unknown{$barcode}++;
        }
        elsif ( scalar @nearest_tags > 1 ) {
            $nearest_tag = q{ };                                 # Multiple tags
            $description = q{ } . join q{|}, sort @nearest_tags;
            $multi++;
        }
        else {
            $nearest_tag = $nearest_tags[0];
        }
        ## no critic (ProhibitMagicNumbers)
        printf { $fh_for{$nearest_tag}{1} } "@%s#%s/%d%s\n%s\n+\n%s\n",
          $read1->qname, $random, 1, $description, $read1_seq, join q{},
          map { chr $_ + 33 } @{ $read1->qscore };
        printf { $fh_for{$nearest_tag}{2} } "@%s#%s/%d\n%s\n+\n%s\n",
          $read2->qname, $random, 2, $read2->query->dna, join q{},
          map { chr $_ + 33 } @{ $read2->qscore };
        ## use critic
    }
}

# Close all filehandles
foreach my $tag ( keys %fh_for ) {
    close $fh_for{$tag}{1};
    close $fh_for{$tag}{2};
}

# Output QC stats
printf "Number of reads matching multiple tags: %d\n", $multi;
print "Unknown tags:\n";    ## no critic (RequireCheckedSyscalls)
my $max;
foreach my $tag ( reverse sort { $unknown{$a} <=> $unknown{$b} } keys %unknown )
{
    if ( !defined $max ) {
        $max = $unknown{$tag};
    }
    last if $unknown{$tag} < $max / 100;    ## no critic (ProhibitMagicNumbers)
    my $tf           = Text::Fuzzy->new($tag);
    my @nearest_tags = $tf->nearestv( \@tags );
    my $distance     = $tf->last_distance();
    printf "%s\t%d\t%s\t%d\n", $tag, $unknown{$tag},
      ( join q{|}, @nearest_tags ), $distance;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'tag_file=s'      => \$tag_file,
        'bam_files=s@{,}' => \@bam_files,
        'edit_distance=i' => \$edit_distance,
        'output_prefix=s' => \$output_prefix,
        'random_start=i'  => \$random_start,
        'random_end=i'    => \$random_end,
        'inread_start=i'  => \$inread_start,
        'inread_end=i'    => \$inread_end,
        'debug'           => \$debug,
        'help'            => \$help,
        'man'             => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$tag_file ) {
        pod2usage("--tag_file must be specified\n");
    }
    if ( !@bam_files ) {
        pod2usage("--bam_files must be specified\n");
    }

    return;
}

__END__
=pod

=encoding UTF-8

=head1 NAME

detag_npg_bam.pl

Detag NPG BAM files

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a name-sorted BAM file and splits it by barcode (BC) using a
supplied list of possible barcodes.

=head1 EXAMPLES

    perl \
        detag_npg_bam.pl \
        --tag_file tags.txt \
        --bam_files sorted.bam

    mkdir foo
    perl \
        detag_npg_bam.pl \
        --tag_file tags.txt \
        --bam_files sorted.bam \
        --edit_distance 2 \
        --output_prefix foo/bar \
        --random_start 9 \
        --random_end 18 \
        --inread_start 1 \
        --inread_end 8

=head1 USAGE

    detag_npg_bam.pl
        [--tag_file file]
        [--bam_files file...]
        [--edit-distance int]
        [--output_prefix path]
        [--random_start int]
        [--random_end int]
        [--inread_start int]
        [--inread_end int]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--tag_file FILE>

File containing sample and barcode pairs (separated by whitespace).

=item B<--bam_files FILE...>

Name-sorted BAM files containing reads with BC tags.

=item B<--edit_distance INT>

Maximum edit distance between tag and barcode.

=item B<--output_prefix PATH>

Prefix added to all output filenames.

=item B<--random_start INT>

Base position where random bases start in read 1.

=item B<--random_end INT>

Base position where random bases end in read 1.

=item B<--inread_start INT>

Base position where in-read tag starts in read 1.

=item B<--inread_end INT>

Base position where in-read tag ends in read 1.

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
