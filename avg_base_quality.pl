#!/usr/bin/env perl

# PODNAME: avg_base_quality.pl
# ABSTRACT: Calculate average base quality

## Author     : ian.sealy@sanger.ac.uk
## Maintainer : ian.sealy@sanger.ac.uk
## Created    : 2016-12-16

use warnings;
use strict;
use autodie;
use Getopt::Long;
use Pod::Usage;
use Carp;
use version; our $VERSION = qv('v0.1.0');

use List::Util qw( sum );

# Default options
my $overall = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Iterate over STDIN and add sum
my @avg_quals;
while ( my $line = <> ) {
    chomp $line;
    $line =~ s/\s//xmsg;
    my @quals;
    foreach my $base ( split //xms, $line ) {
        push @quals, ord($base) - 33;    ## no critic (ProhibitMagicNumbers);
    }
    my $row_avg = q{};
    if ( scalar @quals ) {
        $row_avg = sum(@quals) / scalar @quals;
    }
    push @avg_quals, $row_avg;
}
if ( !$overall ) {
    foreach my $row_avg (@avg_quals) {
        if ( length $row_avg ) {
            printf "%.1f\n", $row_avg;
        }
        else {
            printf "\n";
        }
    }
}
else {
    @avg_quals = grep { length } @avg_quals;
    my $overall_avg = q{};
    if ( scalar @avg_quals ) {
        $overall_avg = sum(@avg_quals) / scalar @avg_quals;
        printf "%.1f\n", $overall_avg;
    }
    else {
        printf "\n";
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'overall' => \$overall,
        'debug'   => \$debug,
        'help'    => \$help,
        'man'     => \$man,
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

avg_base_quality.pl

Calculate average base quality

=head1 VERSION

version 0.1.0

=head1 DESCRIPTION

This script takes a file of base qualities on STDIN and calculates the average
per row or, optionally, overall.

=head1 EXAMPLES

    perl \
        avg_base_quality.pl \
        < input.txt > output.txt

=head1 USAGE

    avg_base_quality.pl
        [--overall]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--overall>

Just calculate overall average base quality.

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
