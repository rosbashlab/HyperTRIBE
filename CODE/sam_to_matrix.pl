#!/usr/bin/env perl
use strict;
use warnings;

#author Joe Rodriguez 2011

$| = 1;
my $file = $ARGV[0];
my $exp = $ARGV[1];
my $timepoint = $ARGV[2];

open (OUT,">$file\.matrix\.wig") or die $!;
open(FH,$file) or die $!;
my $matrix = {};
my ($left,$right) = (0,0);
my $prevchr = "NA";
my $reads;
while(<FH>){
    chomp;
    my $line = $_;
    if($line =~ /^\@/){next;}
    my($ID,$chr,$start,$matchinfo,$seq) = (split(/\t/,$line))[0,2,3,5,9];
    
    if(($prevchr ne "NA") && ($chr ne $prevchr)){
	foreach my $pos (sort {$a <=> $b}keys %{$matrix}){
	    my $A = $matrix->{$pos}->{"A"};
	    my $T = $matrix->{$pos}->{"T"};
	    my $C = $matrix->{$pos}->{"C"};
	    my $G = $matrix->{$pos}->{"G"};
	    my $N = $matrix->{$pos}->{"N"};
	    my $total = $A + $T + $C + $G + $N;
	    print OUT "$exp\t$timepoint\t$prevchr\t$pos\t$A\t$T\t$C\t$G\t$N\t$total\n";
	    delete($matrix->{$pos});
	}
	$matrix = {};
    }
    
    my @datapts = split('',$seq);
    if($matchinfo =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/){
	my($left,$span1,$middle1,$span2,$middle2,$span3,$right) = ($1,$2,$3,$4,$5,$6,$7);
	for (my $j = 0; $j <= $#datapts; $j++){
	    my $base = $datapts[$j];	
	    if($j < $left){
		my $bp = $start + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	    elsif($j < ($left + $middle1)){
		my $bp = $start + $span1 + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	    elsif($j < ($left + $middle1 + $middle2)){
		my $bp = $start + $span1 + $span2 + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	    elsif($j < ($left + $middle1 + $middle2 + $right)){
		my $bp = $start + $span1 + $span2 + $span3 + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	}
    }
    elsif($matchinfo =~ /^(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M$/){
	my($left,$span1,$middle,$span2,$right) = ($1,$2,$3,$4,$5);
	for (my $j = 0; $j <= $#datapts; $j++){
	    my $base = $datapts[$j];	
	    if($j < $left){
		my $bp = $start + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	    elsif($j < ($left + $middle)){
		my $bp = $start + $span1 + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	    elsif($j < ($left + $middle + $right)){
		my $bp = $start + $span1 + $span2 + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	}
    }
    elsif($matchinfo =~ /^(\d+)M(\d+)N(\d+)M$/){
	my($left,$span,$right) = ($1,$2,$3);
	for (my $j = 0; $j <= $#datapts; $j++){
	    my $base = $datapts[$j];	
	    if($j < $left){
		my $bp = $start + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	    elsif($j < ($left + $right)){
		my $bp = $start + $span + $j;
		$matrix->{$bp}->{$base}++;
		foreach my $b ("A","T","C","G","N"){
		    $matrix->{$bp}->{$b} += 0;
		}
	    }
	}
    }
    elsif($matchinfo =~ /^(\d+)M$/){
	for (my $j = 0; $j <= $#datapts; $j++){
	    my $base = $datapts[$j];	
	    my $bp = $start + $j;
	    $matrix->{$bp}->{$base}++;
	    foreach my $b ("A","T","C","G","N"){
		$matrix->{$bp}->{$b} += 0;
	    }
	}
    }
    else{
#		print "GREPOUT $line\n";
    }
    foreach my $pos (sort {$a <=> $b}keys %{$matrix}){
	if($pos < $start){
	    my $A = $matrix->{$pos}->{"A"};
	    my $T = $matrix->{$pos}->{"T"};
	    my $C = $matrix->{$pos}->{"C"};
	    my $G = $matrix->{$pos}->{"G"};
	    my $N = $matrix->{$pos}->{"N"};
	    my $total = $A + $T + $C + $G + $N;
	    print OUT "$exp\t$timepoint\t$prevchr\t$pos\t$A\t$T\t$C\t$G\t$N\t$total\n";
	    delete($matrix->{$pos});
	}
	else{
	    last;
			}
    }
    
    
    $prevchr = $chr;
    
#print "Processed\t$ID\n";
}
close FH;

foreach my $pos (sort {$a <=> $b}keys %{$matrix}){
    my $A = $matrix->{$pos}->{"A"};
    my $T = $matrix->{$pos}->{"T"};
    my $C = $matrix->{$pos}->{"C"};
    my $G = $matrix->{$pos}->{"G"};
    my $N = $matrix->{$pos}->{"N"};
    my $total = $A + $T + $C + $G + $N;
    print OUT "$exp\t$timepoint\t$prevchr\t$pos\t$A\t$T\t$C\t$G\t$N\t$total\n";
}
close OUT;
