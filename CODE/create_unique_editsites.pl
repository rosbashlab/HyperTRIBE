#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
# remove duplicate editing sites from overlapping genes on same strand
# this will reduce your list of target genes, so use this code with Caution
my ($bedgraph_file)= $ARGV[0];
my $hash = {};

open(my $BEDFILE, "<", $bedgraph_file) 
    or die "unable to open file $bedgraph_file";
while ( my $line = <$BEDFILE> ) {
    chomp $line;
#    print "$line\n";
    next if ($line=~/^[\s\n]/);    
    # skip header trackname line
    next if ($line=~/^trackname/);   
    my @arr = split(/\t/, $line);
    my $loci = $arr[23];
#    my $loci = $arr[26];
#    $loci =~s/\s//g;
    if (exists $hash->{$loci}) {
	next;
    } else {
	$hash->{$loci}++;
	print "$line\n";
    }
 #   $line =~s/ //g;
 #   print "\n|$loci|\n$line\n";
#    die;
#    push @{$arr2d}, \@arr;
 
}
close $BEDFILE;

exit;
