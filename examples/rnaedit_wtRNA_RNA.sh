#!/bin/sh

HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"
annotationfile="/path_from_root/HyperTRIBE/annotations/refFlat.txt"
# The combination of tablename, expt name and tp value is used to extract the base composition at a given location of the genome for each experimental condition. the find_rnaeditsites.pl script uses these variables to extract and compare the base composition between the rna library and gDNA library to call the edit sites
#---------------------------
# edit the following varibales as need
#annotationfile="/path_from_root/HyperTRIBE/annotations/refFlat.txt"
wtRNAtablename="exampleDB"
wtRNAexp="wtRNA"
wtRNAtp="1"
RNAtablename="exampleDB"
RNAexp="HyperTRIBE"

#edit the timepoint array as needed
timepoint=(2)
# one or more samples can be run by altering this array. 
#---------------

for tp in ${timepoint[@]}
do
  outfile1=$RNAexp"_"$wtRNAtp"_"$tp"_A2G.txt"
  perl $HyperTRIBE_DIR/find_rnaeditsites.pl -a $annotationfile -t $RNAtablename -e $RNAexp -c $tp -o $outfile1 -g $wtRNAtablename -j $wtRNAexp -k $wtRNAtp  
# filter editsites based on cut off 20 reads, it is hard coded in the python script
# to change the threshold, chnage the value in the python script TotalCountThreshold = 10
#
  python $HyperTRIBE_DIR/Threshold_editsites_20reads.py $outfile1
  outfile2=$outfile1".threshold"

  python $HyperTRIBE_DIR/convert_editsites_to_bedgraph.py $outfile2 
  prefix=${outfile2%.txt.threshold*}
  mv $outfile2".bedgraph" $prefix".bedgraph"
#  echo $outfile2 $prefix".bedgraph"

done
