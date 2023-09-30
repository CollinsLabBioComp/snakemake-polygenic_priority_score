#!/usr/bin/env perl
# $Id: merge.plx 5970 2016-02-06 18:39:02Z pchines $

use strict;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

merge.plx - merge files together based on common key column

=head1 SYNOPSIS

merge.plx [options] file1name:col file2name:col

=head1 DESCRIPTION

This program will merge together delimited text files; both files must use
the same delimiter, in this version.  In order to merge the
files there must be one column in each file where the values correspond to
one another.  The output will contain a line with the concatenated values
from both files where the files all have the same key value in common.  Thus,
merge performs an intersection between the files, and outputs a single line
for each intersection, containing all of the data from all of the files.

Typically each key value will appear only once in each file.  In this version
of the merge program, however, multiple values may be handled in one of three
ways.  See options B<-a>, B<-f> and B<-l> below.  If none of these options
are specified and multiple values for the same key are encountered, the
program will terminate with an error.

=cut

use vars qw(@FILES %KEY_COLUMN $DELIM $CASE_INSENSITIVE $IGNORE_BLANKS
            $USE_FIRST $USE_LAST $USE_ALL $OUTER_JOIN
            $OrderFrom $KEEP_ALL_COLS $KeyIncluded);

process_commandline();
my $rh_cumulative;
my $ra_order;
foreach my $file (@FILES) {
    my ($rh,$ra) = read_file($file, $KEY_COLUMN{$file});
    $rh_cumulative = intersect_and_concatenate($rh_cumulative, $rh);
    if ($OrderFrom eq $file) {
        $ra_order = $ra;
    }
}
if (!$ra_order) {
    my @order = sort keys %$rh_cumulative;
    $ra_order = \@order;
}
elsif ($OUTER_JOIN) {
    my %got;
    @got{@$ra_order} = ();
    for my $key (sort keys %$rh_cumulative) {
        if (!exists $got{$key}) {
            push @$ra_order, $key;
            $got{$key} = 1;
        }
    }
}
foreach my $key (@$ra_order) {
    if (ref $rh_cumulative->{$key}) {
        print "$_\n" foreach (@{ $rh_cumulative->{$key} });
    }
    elsif (defined $rh_cumulative->{$key}) {
        print $rh_cumulative->{$key}, "\n";
    }
}

## End MAIN

sub intersect_and_concatenate {
    my ($rh_cum, $rh_add) = @_;
    my $nfields = get_field_count($rh_add);
    my $cfields = get_field_count($rh_cum);
    my %new;
    if ($rh_cum) {
        my %allkeys;
       @allkeys{keys %$rh_cum} = ();
        @allkeys{keys %$rh_add} = () if $OUTER_JOIN > 1;
        foreach my $key (keys %allkeys) {
            if (exists $rh_cum->{$key} && exists $rh_add->{$key}) {
                $new{$key} = concatenate($rh_cum->{$key}, $rh_add->{$key});
            }
            elsif (exists $rh_cum->{$key} && $OUTER_JOIN) {
                $new{$key} = concatenate($rh_cum->{$key},
                    $DELIM x ($nfields-1) );
            }
            elsif ($OUTER_JOIN > 1) {
                $new{$key} = concatenate($DELIM x ($cfields-1),
                        $rh_add->{$key});
            }
        }
    }
    else {
        %new = %$rh_add;
    }
    return \%new;
}

sub concatenate {
    my ($val1, $val2) = @_;
    my $newval = [];
    if (ref($val1) && ref($val2)) {
        foreach my $v1 (@$val1) {
            foreach my $v2 (@$val2) {
                push @$newval, join($DELIM, $v1, $v2);
            }
        }
    }
    elsif (ref $val1) {
        foreach my $v1 (@$val1) {
            push @$newval, join($DELIM, $v1, $val2);
        }
    }
    elsif (ref $val2) {
        foreach my $v2 (@$val2) {
            push @$newval, join($DELIM, $val1, $v2);
        }
    }
    else {
        $newval = join($DELIM, $val1, $val2);
    }
    return $newval;
}

sub read_file {
    my ($file, $ra_col) = @_;
    my %data;
    open(FILE, "<$file") || die "Can't open $file, $!\n";
    local $/ = guess_newline(*FILE);
    my $nfields = count_fields(*FILE);
    if ($nfields == 1) {
        warn "Warning: only one column found in '$file';\n"
           . "         perhaps delimiter is wrong?\n";
    }
    if ($nfields < max(@$ra_col)) {
        die "Error: there are only $nfields columns in '$file';\n"
          . "       can't use column @$ra_col to merge\n";
    }
    if ($KeyIncluded) {
        $nfields -= @$ra_col;
    }
    my @keys;
    while (<FILE>) {
        chomp;
        my @data = split /$DELIM/;
        my $key = join($DELIM, map { $data[$_-1]||'' } @$ra_col);
        if ($CASE_INSENSITIVE) {
            $key = uc($key);
        }
        if ($IGNORE_BLANKS) {
            $key = s/^\s+//;
            $key = s/\s+$//;
            $key = s/\s*$DELIM\s*/$DELIM/g;
        }
        if ($KeyIncluded) {
            for my $col (sort { $b<=>$a } @$ra_col) {
                splice(@data,$col-1,1);
            }
            $_ = join($DELIM, @data);
        }
        my $line = fixed_fields($_, $nfields);
        if ($data{$key}) {
            if ($USE_FIRST) {
                # Leave current data alone
            }
            elsif ($USE_LAST) {
                $data{$key} = $line;
            }
            elsif ($USE_ALL) {
                if (ref $data{$key}) {
                    push @{ $data{$key} }, $line;
                }
                else {
                    $data{$key} = [ $data{$key}, $line ];
                }
            }
            else {
                die "Encountered key value '$key' multiple times at line "
                    .__LINE__." of file '$file'.\n"
                    ."Use option -a, -f or -l to determine how to handle.\n";
            }
        }
        else {
            $data{$key} = $line;
            push @keys, $key;
        }
    }
    close FILE;
    $KeyIncluded = 1 if !$KEEP_ALL_COLS;
    return (\%data, \@keys);
}

sub guess_newline {
    my ($fh) = @_;
    my $newline;
    my $text;
    if (read($fh, $text, 4096)) {
        if ($text =~ m/(\015?\012|\015)/) {
            $newline = $1;
        }
    }
    seek $fh, 0, 0;
    return $newline || "\n";
}

sub count_fields {
    my ($fh) = @_;
    my $max_n = 0;
    while (<$fh>) {
        my $n = eval "tr/$DELIM/$DELIM/ + 1";
        $max_n = $n if $n > $max_n;
    }
    seek $fh, 0, 0;
    return $max_n || 1;
}

sub get_field_count {
    my ($rh) = @_;
    my ($text) = values %$rh;
    if (ref $text) {
        $text = $text->[0];
    }
    return 0 if !defined $text;
    my $nfields = eval "(\$text =~ tr/$DELIM/$DELIM/) + 1";
    return $nfields;
}

sub fixed_fields {
    my ($line, $nfields) = @_;
    my $f = eval "(\$line =~ tr/$DELIM/$DELIM/) + 1";
    $line .= $DELIM x ($nfields - $f);
    return $line;
}

sub max {
    my $max = shift;
    for (@_) {
        if ($_ > $max) {
            $max = $_;
        }
    }
    return $max;
}

sub process_commandline {
    my ($help, $version, $order_from);
    $DELIM = "\t";
    $USE_ALL = $USE_FIRST = $USE_LAST = 0;
    $OUTER_JOIN = 0;
    $OrderFrom = '';
    GetOptions(
        'help'      => \$help,
        'version'   => \$version,
        'delimiter=s' => \$DELIM,
        'insensitive' => \$CASE_INSENSITIVE,
        'blanks'    => \$IGNORE_BLANKS,
        'all'       => \$USE_ALL,
        'first'     => \$USE_FIRST,
        'last'      => \$USE_LAST,
        'outer+'    => \$OUTER_JOIN,
        'keep'      => \$KEEP_ALL_COLS,
        'sort|n:i'  => \$order_from,
        ) || pod2usage(1);
    if ($help)    { pod2usage(verbose => 2); }
    if ($version) { die $0, ", ", q$Revision: 5970 $, "\n"; }
    pod2usage(1) if $USE_ALL + $USE_FIRST + $USE_LAST > 1;

    @FILES = @ARGV;
    if (@FILES < 2) {
        pod2usage("Need at least two files to merge\n");
    }
    for my $file (@FILES) {
        pod2usage("Must supply column numbers of key columns")
            if $file !~ /:\d+(,\d+)*$/;
        pod2usage("No other ':' allowed in file names")
            if $file =~ /:.*:\d+(,\d+)*$/;
    }
    %KEY_COLUMN = split(/:/, join(":", @FILES));
    my $n;
    while (my ($k,$v) = each %KEY_COLUMN) {
        my @c = split /,/, $v;
        $KEY_COLUMN{$k} = \@c;
        $n ||= @c;
        if (@c != $n) {
            die "All files must have same number of key fields\n";
        }
    }
    @FILES = map { s/:\d+(,\d+)*$//;$_ } @FILES;
    if (defined $order_from) {
        if ($order_from) {
            --$order_from;
        }
        $OrderFrom = $FILES[$order_from]
            || die "Can't take order from file #" . ($order_from+1) . "; only "
                . scalar(@FILES) . " files provided.\n";
    }
}

=head1 OPTIONS

=over 4

=item B<--all>

Where there are multiple lines that share the same key value, keep all of
them in the output.  This results in a Cartesian product, i.e. if there are
two rows in the first file and three in the second file, there will be six
rows in the output.  Mutually exclusive with options B<-f> and B<-l>.

=item B<--blanks>

Ignore leading and trailing blanks in key fields.

=item B<--delimiter> string

Set the column delimiter (defaults to tab).

=item B<--first>

Where there are multiple lines that share the same key value, keep only the
first of these encountered in the input.  Mutually exclusive with options
B<-a> and B<-l>.

=item B<--help>

Display full manual page.

=item B<--insensitive>

Match key columns without regard to case.

=item B<--keep>

Keep all key columns.  The (new as of version 1.14) default behavior is to
keep the first instance of this column only.

=item B<--last>

Where there are multiple lines that share the same key value, keep only the
last of these encountered in the input.  Mutually exclusive with options
B<-a> and B<-f>.

=item B<--outer>

Perform an outer join, i.e. include all records from earlier files, even if
a matching record does not appear in later files.

Use this option twice, eg.  C<merge.plx -o -o> in order to include all
records from all files, even if no match is found in the others.

=item B<--sort> [FILE_NUMBER]

Sort results by the order of the specified file.  If no file number is
provided, "1" is assumed, to sort by the order of the first file.

=back

=head1 TODO

Think about renaming to join.plx and supporting command line options used by
textutils join command (which works fine, except that it requires files to be
sorted), and only allows two files at a time.

=head1 AUTHOR

Peter Chines - pchines@nhgri.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
FUSION Research Group" as the citation. 

=cut
