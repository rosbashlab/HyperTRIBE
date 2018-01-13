#!/usr/bin/env perl
use strict;
use warnings;

#Author Reazur Rahman, January 2018
#usage perl summarize_results.pl present_both.bedgraph > summary_output.xls

#------------------------
# Output  file format
# col 1: gene name
# col 2: number of editing sites for gene
# col 3: Avg editing percentage
# col 4: edit_percentage_read concatenation for each replicate, separated by comma (13%_29r,13%_29r). ";" is used to separate between different editing sites
# col 5: gene feature concatenation, separated by comma (EXON,INTRON). ";" is used to separate between different editing sites
# col 6: Edit identifier concatenation, separated by comma. ";" is used to separate between different editing sites
#
#--------------------
#IMPORTANT PLEASE READ: gene name for all replicates at for a given edit site should be identical, otherwise these are sites from overlapping genes. These cross overlapping sites are skipped because they are counted twice. Hence, the number of edit sites reported by this script may not equal the number of lines in the bedgraph file.
# ------ input bedgraph column header --------
#1. Chr name
#2. Start coordinate
#3. End coordinate
#4. Editing percentage
#5. Concatenation of editing percentage and reads (Total nucleotide count for RNA)
#6. Chr name
#7. Edit Coordinate
#8. gene Name
#9. Type
#10. A count
#11. T count
#12. C count
#13. G count
#14. Total nucleotide count
#15. A count from gDNA/wtRNA
#16. T count from gDNA/wtRNA
#17. C count from gDNA/wtRNA
#18. G count from gDNA/wtRNA
#19. Total count from gDNA/wtRNA
#20. Editbase Count
#21. Total nucleotide count
#22. Editbase count from gDNA/wtRNA
#23. Total nucleotide count from gDNA/wtRNA.
#24. Identifier (chr name and coordinate)
#------------------------

#open the input bedgraph file
my ($bedgraph_file)= $ARGV[0];
open(my $BEDFILE, "<", $bedgraph_file) 
    or die "unable to open file $bedgraph_file";

my $hash = {};
my $offset = 24; 
while ( my $line = <$BEDFILE> ) {
    chomp $line;
#    print "$line\n";
    next if ($line=~/^[\s\n]/);    
    # skip header trackname line
    next if ($line=~/^trackname/);   

    my @arr = split(/\t/, $line);
    my $arr_length =  @arr;
    my $num_replicate = int( $arr_length/$offset);
    my $num_loop = $num_replicate -1;
#    print "num_rep= $num_replicate; num_loop=$num_loop\n";
    my $gene_hash = {};
    my $edit_percentage_read_str = "";
    my $site_identifier = "";
    my $sum_editing_percentage = 0;
    my $gene_features = "";
    foreach my $tp (0..$num_loop) {
	my $gene_name = $arr[($offset)*$tp + 7];
	my $feature = $arr[($offset)*$tp + 8];
	my $editpercentage = $arr[($offset)*$tp + 3];
	my $edit_read_str = $arr[($offset)*$tp + 4];
	my $edit_identifier = $arr[($offset)*$tp + 23];
#create str and hash for printing
	$gene_hash->{$gene_name}++;
	$site_identifier = $edit_identifier;
	$edit_percentage_read_str .= "$edit_read_str,";
	$sum_editing_percentage += $editpercentage;
	$gene_features .= "$feature,";
#if the number of gene name is one, then add to number of edit sites

#	print "gene_name: $gene_name; editpercentage: $editpercentage; edit_read_str: $edit_read_str; edit_identifier: $edit_identifier\n";	
    }
    
#if the number of gene is more than, it means it is picking the overlapping genes, so it can be skipped
    next if (keys %$gene_hash> 1);    
#if the number of gene name is one, then add to number of edit sites
    my $gene_name;
    foreach my $key (keys %$gene_hash) {
	$gene_name=$key;
    }
    my $avg_editing = $sum_editing_percentage/$num_replicate;
    $avg_editing= sprintf( "%.1f", $avg_editing);
    $edit_percentage_read_str=~s/,$//;
    $gene_features=~s/,$//;
    my $str_arr = [$gene_name, $avg_editing, $site_identifier, $edit_percentage_read_str, $gene_features];
    
    push @{$hash->{$gene_name}}, $str_arr; 
    #$hash_count->{$gene_name}=++;
    
#    print "$gene_name, $avg_editing, $site_identifier, $edit_percentage_read_str\n";
#    die;
}
close $BEDFILE;
my $print_arr = [];
foreach my $gene (sort keys %$hash ) {
    my $gene_name = $gene;
    my $line_arr = $hash->{$gene};
    my $num_edit_sites = @$line_arr;
    my ($avg_editing,  $site_identifier, $edit_percentage_read_str, $gene_features, $sum_editing);
    foreach my $line (@$line_arr) {
	$sum_editing += "$line->[1]";
#	$avg_editing .= "$line->[1];";
	$site_identifier .= "$line->[2];";
	$edit_percentage_read_str .= "$line->[3];";
	$gene_features .= "$line->[4];";
    }
   $avg_editing= $sum_editing/$num_edit_sites;
   $avg_editing= sprintf( "%.1f", $avg_editing);
#    $avg_editing=~s/;$//;
    $site_identifier=~s/;$//;
    $edit_percentage_read_str=~s/;$//;
    $gene_features=~s/;$//;
    #print "$gene_name\t$num_edit_sites\t$avg_editing\t$gene_features\t$edit_percentage_read_str\t$site_identifier\n";
    push @$print_arr, [$gene_name, $num_edit_sites, $avg_editing, $gene_features, $edit_percentage_read_str, $site_identifier];
    
}
#print header for the file
print "Gene_name\tNum_edit_sites\tAvg_edting\tFeatures\tEdit_percent_read_str\tIdentifier_str\n";
foreach my $line (sort { $b->[1] <=> $a->[1] } @$print_arr) {
#    print "@$line\n";
    my $print_str = join("\t", @$line) ;
    print "$print_str\n";
}

exit;
