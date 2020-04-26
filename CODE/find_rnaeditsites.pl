#!/usr/bin/env perl
use strict;
use warnings;
use DBI; 
use Getopt::Std;

$| = 1;

# Author: First version written by Joe Rodriguez in 2011. In subsequent years, this script has been overhauled by Reazur Rahman
# current version written by Reazur Rahman on 01/11/2018

# Identifies RNA editing sites from RNA sequencing data relative to a genomic dataset.  
# Bed file data containing read counts of each nucleotide per genomic coordinate must be uploaded to mysql database.
# Editing sites will be written to an user defined output file (-o flag) 
# The gDNA/wtRNA is used as reference to call edit sites

my $USAGE = "find_rnaeditsites.pl -a annotationfile -t RNAtablename -e RNAexp -c tp -o outfile -g gDNA_tablename -j gDNA_exp -k gDNA_tp\n OR find_rnaeditsites.pl -a annotationfile -t RNAtablename -e RNAexp -c tp -o outfile -g wtRNA_tablename -j wtRNA_exp -k wtRNA_tp\n"; 
my %option;
getopts( 'a:t:e:g:c:o:j:k:h', \%option );

my ($annotationfile, $tablename, $exp, $gDNAtablename, $OUTFILE, $tp, $gexp, $gDNAtp);

#MYSQL CONFIG VARIABLES
my $host = "localhost";
my $database = "dmseq";
my $user = "root"; #mysql username
my $password = ""; #mysql password, if any
#my $password = "password"; #mysql password, if any

#connect to the mysql database
my $dsn = "DBI:mysql:$database:$host"; 
my $dbh = DBI->connect( $dsn, $user, $password, {RaiseError=>1, PrintError=>0} ) or die $DBI::errstr;

my ($sth);
#-----------------
#MYSQL TABLE VARIABLES
my $mincovthresh = '9'; 	#Minimum genome coverage must be greater than this number.  

# if help is need to remember how to run the 
if ( $option{h} ) { print STDERR "$USAGE\n\n"; exit 1; }

if ( $option{a} && $option{t} && $option{g} && $option{e} && $option{o} && $option{j}) {
    $annotationfile = $option{a};
    $tablename = $option{t};
    $exp = $option{e};
    $gDNAtablename = $option{g};
    $OUTFILE = $option{o};
    $gexp = $option{j};    
    if (exists $option{k}) {
	$gDNAtp =  $option{k};
    } else {
	die "value for gDNAtp is not being set";
    }
    if (exists $option{c}) {
        $tp =  $option{c};
    } else {
        die "value for tp is not being set";
    }
} else {
    print "annotationfile = $option{a};\n
    tablename = $option{t};\n
    exp = $option{e};\n
    gDNAtablename = $option{g};\n
    OUTFILE = $option{o};\n
    tp = $option{c};\n
    gexp = $option{j};\n 
    gDNAtp = $option{k};\n";

    die "all the variables are not being passed properly at n04rnaeditID.pl";
}

#MYSQL EDITING VARIABLES
my %gstrings = ("A"=>"g.Acount","T"=>"g.Tcount","C"=>"g.Ccount","G"=>"g.Gcount");
my %rstrings = ("A"=>"a.Acount","T"=>"a.Tcount","C"=>"a.Ccount","G"=>"a.Gcount");

my $noneditbase = "A";
my $editbase = "G";
my $noneditbaseREV = "T";
my $editbaseREV = "C";

my %complement = ("A"=>"T","T"=>"A","C"=>"G","G"=>"C");

my $filesuffix = $noneditbase . $editbase;

# read the annotatation file is in refflat format, with "all field from selected tables"
# parse the annotation file and return a hash of genes
my ($genes) = &parse_annotation_file($annotationfile);

open(OUT,">$OUTFILE") or die $!;

#OUTPUT ANY INSTANCE OF AN EDITING EVENT
print OUT "Chr\tEdit_coord\tName\tType\tA_count\tT_count\tC_count\tG_count\tTotal_count\tA_count_gDNA/wtRNA\tT_count_gDNA/wtRNA\tC_count_gDNA/wtRNA\tG_count_gDNA/wtRNA\tTotal_count_gDNA/wtRNA\tEditbase_count\tTotal_count\tEditbase_count_gDNA/wtRNA\tTotal_count_gDNA/wtRNA\n";

foreach my $chr (keys %{$genes}){
    foreach my $gene (keys %{$genes->{$chr}}){
	my $exoncount = keys %{$genes->{$chr}->{$gene}};
	my $currexoncount = 0;
	my @genearray;
	my $strand;
	my ($left,$right);
	foreach my $start (sort {$a <=> $b} keys %{$genes->{$chr}->{$gene}}){
	    $currexoncount = $currexoncount + 1;
	    my $end = $genes->{$chr}->{$gene}->{$start}->{END};
	    my $type = $genes->{$chr}->{$gene}->{$start}->{TYPE};
	    my $ID = $genes->{$chr}->{$gene}->{$start}->{ID};
	    $strand = $genes->{$chr}->{$gene}->{$start}->{STRAND};
	    if(!defined($left)){
		$left = $start;
		$right = $end;
	    }
	    else{
		$left = min($left,$start);
		$right = max($right,$end);
	    }	
	    push(@genearray,[$start,$end]);
	}
	my @gene_data_TP;
	my $coding_sum = {};
	my $gene_size = 0;
	my $data = {};
	my $Rconsdata = {}; #rnadata
	my $Gconsdata = {}; #gDNAdata
	my $exonbases = {};
	for(my $i =0; $i <= ($#genearray); $i++){ #EACH EXON
	    my ($s0,$e0) = ($genearray[$i]->[0],$genearray[$i]->[1]);
	    foreach my $bp ($s0..$e0){
		$exonbases->{$bp} = 1;
	    }
	}
	
	my $myquery1 = "SELECT a.*, g.* from $tablename as a, $gDNAtablename as g where a.experiment IN ('" . $exp . "') and a.time = '" . $tp . "' and a.chr = '" . $chr . "' and a.coord between '" . $left . "' and '" . $right . "'";
	$myquery1 .=     " and g.experiment='" . $gexp . "' and g.time='" . $gDNAtp . "' and g.chr=a.chr and g.coord=a.coord and g.totalcount>'" . $mincovthresh . "'";
#	if($strand eq "+"){$myquery1 .= " and $gstrings{$editbase} < '1' and ($gstrings{$noneditbase}/g.Totalcount) >= '0.8' and $rstrings{$editbase} > '0'";}
#	elsif($strand eq "-"){$myquery1 .= " and $gstrings{$editbaseREV} < '1' and ($gstrings{$noneditbaseREV}/g.Totalcount) >= '0.8' and $rstrings{$editbaseREV} > '0'";}

	if($strand eq "+"){$myquery1 .= " and ($gstrings{$editbase}/g.Totalcount) < '0.005' and ($gstrings{$noneditbase}/g.Totalcount) >= '0.8' and $rstrings{$editbase} > '0'";}
	elsif($strand eq "-"){$myquery1 .= " and ($gstrings{$editbaseREV}/g.Totalcount) < '0.005' and ($gstrings{$noneditbaseREV}/g.Totalcount) >= '0.8' and $rstrings{$editbaseREV} > '0'";}
	
# execute the query
	eval {
	    $sth = $dbh->prepare($myquery1); 
	    $sth->execute;
	};
	
	if ($@) {
	    print STDERR "Error in executing mysql query: $@";
	}

	&print_query_result($sth, $exonbases, $strand, $data, $Gconsdata, $chr, $gene);
    } #foreach my $gene
} #chr

#close connection to database
$sth->finish; 
$dbh->disconnect;

close OUT;
print "Wrote $OUTFILE\n";		
exit;


sub max {
    my $n1 = shift;
    my $n2 = shift;
    if($n1 > $n2){
	return $n1;
}
    else{
	return $n2;
    }
}

sub min {
    my $n1 = shift;
    my $n2 = shift;
    if($n1 < $n2){
	return $n1;
    }
    else{
	return $n2;
    }
}


#---------------------------------
# Print results from the SQL query of edit sites
sub print_query_result {
    my ($execute1, $exonbases, $strand, $data, $Gconsdata, $chr, $gene) = @_; 
 
    while (my @results = $execute1->fetchrow()) {
	my $coordtype = "EXON";
	my $coord = "$results[3]";
	# if the coordinate is not present in an exon, then it is classified as intron. since it is within gene bounds, cannot be intergenic regions
	if(!exists($exonbases->{$coord})){
	    $coordtype = "INTRON";
	}	
	
	my $exp = $results[0];
	my $tp = $results[1];
	my $A = $results[4];
	my $T = $results[5];
	my $C = $results[6];					
	my $G = $results[7];
	my $total = $results[9];
	
	my $gA = $results[14];
	my $gT = $results[15];
	my $gC = $results[16];					
	my $gG = $results[17];
	my $gtotal = $results[19];
	#RNA counts for that position
	$data->{$coord}->{'A'} = $A;
	$data->{$coord}->{'T'} = $T;
	$data->{$coord}->{'C'} = $C;
	$data->{$coord}->{'G'} = $G;
	$data->{$coord}->{'tot'} = $total;
	#gDNA counts for that position
	$Gconsdata->{$coord}->{'A'} = $gA;
	$Gconsdata->{$coord}->{'T'} = $gT;
	$Gconsdata->{$coord}->{'C'} = $gC;
	$Gconsdata->{$coord}->{'G'} = $gG;
	$Gconsdata->{$coord}->{'tot'} = $gtotal;		
	
	my ($e,$etot);
	my $edit = $editbase;
	
	if($strand eq "-"){
	    $edit = $complement{$editbase};
	}
	
	my $string;
	foreach my $base ('A','T','C','G','tot'){						
	    $data->{$coord}->{$base} += 0;
	    $Gconsdata->{$coord}->{$base} += 0;
	}
	
	foreach my $base ('A','T','C','G','tot'){						
	    $string .= "$data->{$coord}->{$base}\t";
	}
	
	$e = "$data->{$coord}->{$edit}";
	$etot = "$data->{$coord}->{'tot'}";
	
	foreach my $base ('A','T','C','G','tot'){
	    $string .= "$Gconsdata->{$coord}->{$base}\t";
	}
	my $eG = $Gconsdata->{$coord}->{$edit};
	my $eGtot = $Gconsdata->{$coord}->{'tot'};
	print OUT "$chr\t$coord\t$gene\t$coordtype\t$string$e\t$etot\t$eG\t$eGtot\n";	
    } # end while ...   
    return;
}


#--------------------
#parse annotation file in refflat format and return a hash with parsed data
sub parse_annotation_file {
    my ($annotationfile) = @_;
    #read the annotatation file should be in refflat format, with "all field from selected tables"
    open(FH,$annotationfile) or die $!;
    my $genes = {};
    while(<FH>){
	chomp;
	my $line = $_;
	next if ($line=~/^[\#\n]/); # skip empty lines or lines with comments
	next if (($line=~/strand/i) || ($line=~/^name/)); # skip header line
	my($gene,$chr,$strand,$sString,$eString) = (split(/\t/,$line))[0,2,3,9,10]; #fly
	next if ($chr=~/chrUextra/); #the annotations in chrUextra are not considered
	next if ($gene=~/^His.*:/); #remove the Histone genes
	#    if ($gene=~/^His.*:/) {
	#	print STDERR "$gene\n";
	#	next;
	#    }
	next if ($gene=~/^snmRNA*:/); #remove the snmRNA
	
	my @starts = split(/\,/,$sString);
	my @ends = split(/\,/,$eString);
	for (my $i = 0; $i <= $#starts;$i++){
	    if(exists($genes->{$chr}->{$gene}->{$starts[$i]})){
		if($genes->{$chr}->{$gene}->{$starts[$i]}->{END} < $ends[$i]){
		    $genes->{$chr}->{$gene}->{$starts[$i]}->{END} = $ends[$i];
		}
	    } else{
		$genes->{$chr}->{$gene}->{$starts[$i]} = {END=>$ends[$i],TYPE=>"CODING",ID=>$gene,STRAND=>$strand};
	    }
	}
    }
    close FH;
    
    return $genes;
    
}
