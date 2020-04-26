#!/bin/sh

HyperTRIBE_DIR="/home/analysis/editing/HyperTRIBE/CODE"

# The combination of tablename, expt name and tp value is used to extract the base composition at a given location of the genome for each experimental condition. the find_rnaeditsites.pl script uses these variables to extract and compare the base composition between the rna library and gDNA library to call the edit sites
#---------------------------
# edit the following varibales as need
annotationfile="/home/analysis/genome/dm6/refFlat.txt"

gDNAtablename="protocol_testDB"
gDNAexp="gDNA"
gDNAtp="25"
RNAtablename="protocol_testDB"
RNAexp="s2mRNA"

#edit the timepoint array as needed
timepoint=(2 3 4 5)
#timepoint=( 2 )
# one or more samples can be run by altering this array. 
#---------------

for tp in ${timepoint[@]}
do
  outfile1=$RNAexp"_"$gDNAtp"_"$tp"_A2G.txt"
  perl $HyperTRIBE_DIR/find_rnaeditsites.pl -a $annotationfile -t $RNAtablename -e $RNAexp -c $tp -o $outfile1 -g $gDNAtablename -j $gDNAexp -k $gDNAtp  

  # Previously used Threshold_editsites_20reads.py script is replaced with filter_by_threshold_without_header.pl 

  # convert to bedgraph format
  prefix=${outfile1%.txt*}
  python $HyperTRIBE_DIR/convert_editsites_to_bedgraph.py $outfile1
  mv $outfile1".bedgraph" $prefix".bedgraph"

  # apply edit % and read threshold, the index is zero based in the perl script
  # 4th col is edit thresold: $threshold
  # 21st column has a read threshold of 10
  # create a 5% threshold edit file
  edit_threshold=5
  read_threshold=20
  prefix=${outfile1%.txt*}
  out_bedgraph=$prefix"_"$edit_threshold"%.bedgraph"
  perl $HyperTRIBE_DIR/filter_by_threshold_without_header.pl 3 $edit_threshold 20 $read_threshold $prefix".bedgraph" > $out_bedgraph


# the following code has been replace with the code above. Empirical data from our lab has shown that keeping the threshold at 5% provides a reliable target list. But, users can change it as needed.

#  python $HyperTRIBE_DIR/Threshold_editsites_20reads.py $outfile1
#  outfile2=$outfile1".threshold"
#
#  python $HyperTRIBE_DIR/convert_editsites_to_bedgraph.py $outfile2 
#  prefix=${outfile2%.txt.threshold*}
#  mv $outfile2".bedgraph" $prefix".bedgraph"
#  echo $outfile2 $prefix".bedgraph"

# ----------------------------
# remove duplicate editing sites from overlapping genes on the same strand
# this will reduce your list of target genes, so use this code with Caution
#  perl $HyperTRIBE_DIR/create_unique_editsites.pl $prefix".bedgraph" > tmp.txt
#  mv tmp.txt $prefix".bedgraph"



done
