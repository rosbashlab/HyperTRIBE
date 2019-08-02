#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;


my $filename= pop @ARGV;
my $hash = {};

#the col name is zero based

while (@ARGV) {    
    my $col = shift @ARGV;
    my $threshold = shift @ARGV;
    $hash->{$col} = $threshold;
}


#now open the file and apply the filters

open(my $INFILE, "<", $filename) 
    or die "unable to open file $filename";

#might need to add an option for header later...
#my $first_line = <$INFILE>;

#print "$first_line";
while ( my $line = <$INFILE> ) {
#	next if ($line=~/^track/);
    chomp $line;
    my @arr = split(/\t/, $line);
    
    my $SKIP = 0;    
    foreach my $col (keys %$hash) {
	#	
	my $threshold_val = $hash->{$col};
	if ( $arr[$col] < $threshold_val ) {
	    $SKIP=1;	   
	}
    }
    
    #skip the line if one of the filter does not pass the test 
    next if ($SKIP); 
    print "$line\n";
    
}

close $INFILE;
