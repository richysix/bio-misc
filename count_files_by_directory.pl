#!/usr/bin/env perl

# PODNAME: count_files_by_directory.pl
# ABSTRACT: Count files by directory

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2017-09-06

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

# Iterate over STDIN and count files per directory
my %count_for;
while ( my $line = <> ) {
    chomp $line;
    my @dirs = split /\//xms, $line;
    pop @dirs;
    my $full_dir = q{};
    foreach my $dir (@dirs) {
        $full_dir .= $dir;
        $count_for{$full_dir}++;
        $full_dir .= q{/};
    }
}

# Output counts
$count_for{q{/}} = $count_for{q{}};
delete $count_for{q{}};
foreach my $dir (
    ## no critic(ProhibitReverseSortBlock)
    sort { $count_for{$b} <=> $count_for{$a} || $a cmp $b }
    ## use critic
    keys %count_for
  )
{
    printf "%d\t%s\n", $count_for{$dir}, $dir;
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

count_files_by_directory.pl

Count files by directory

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a file on STDIN containing filenames and counts the number of
files per directory. The final counts are the number in that directory and
below.

=head1 EXAMPLES

    perl \
        keep_best_duplicate.pl \
        < files.txt > counts.txt

=head1 USAGE

    count_files_by_directory.pl
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

This software is Copyright (c) 2017 by Genome Research Ltd.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
