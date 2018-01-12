#!/bin/sh
HyperTRIBE_DIR="/home/analysis/editing/HyperTRIBE/CODE"

#the sam file from the previous step
samfile=$1
#tablename for mysql table
tablename=$2
#expt name (identifier for an experiment), choose something short
expt=$3
#unique replicate number or timepoint for an experiment
tp=$4

echo "Load SAM file to MySQL tables with the following parameters:"
echo "SAMFILE: $samfile"
echo "TABLENAME: $tablename"
echo "EXPT NAME: $expt"
echo "TP: $tp"

prefix=${samfile%.sam*}
#create the matrix file
perl $HyperTRIBE_DIR/sam_to_matrix.pl $samfile $expt $tp 
matrix_file=$prefix".matrix"
mv $samfile".matrix.wig" $matrix_file

#load the matrixfile to mysql database.
#the mysql database associated variable for this perl script needs to be updated before running this step 
perl $HyperTRIBE_DIR/load_matrix_data.pl -t $tablename -d $matrix_file


