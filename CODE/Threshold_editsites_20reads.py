#!/usr/bin/env python

# Set the input file name
# (The program must be run from within the directory 
# that contains this data file)
Usage = """

Used for Thresholding of (Input:  RNAedit04_output)
Aoife 11_July_2013

Command line:  python Threshold_editsites_20reads.py [input_file1] [input_file2] [input_file3] [etc] &

Removes 'edited sites' below threshold
1) remove entries <20 reads  (column 8 - Total)
2) remove entries < 0.10 edited (edit/total) (column 14/column 15) (editbasecount/Total)


Appends '.threshold' to the output files

To change thresholds change these variables in the script:
                TotalCountThreshold = 20
                PercentEditThreshold = 0.10 (10%)

Structure of RNAedit04_output(input file):
0	1	2	3	4	5	6	7	8	9	10	11	12	13
chr	coord	CG	type	A	T	C	G	Total	AgDNA	TgDNA	CgDNA	GgDNA	TotalgDNA	

14		15	16			17
editbasecount	Total	editbaseGenomecount	GenomeTot


0               1       2       3       4       5 
chr3RHet	2506585	CG40198	INTRON	0	0


(Modification of :convert_editOut_to_bed_2.py)
"""

import sys

# These are the Thresholds, entries below these values will be removed.  Edit here if necessary
TotalCountThreshold = 20
PercentEditThreshold = 0.1

if len(sys.argv)<2:  # the first line of the arguement is always the name of the program
        print Usage
else:
        FileList = sys.argv[1:]
        
        for InfileName in FileList:
                Infile = open (InfileName, 'r')
                
                OutFileName = InfileName + '.threshold'
                OutFile = open(OutFileName, 'w')
                
                # Initialize the counter used to keep track of line numbers
                LineNumber = 0
                # Loop through each line in the file, skipping the header line - assumes one header line.
                
                for Line in Infile:
                        if LineNumber == 0:
                                #print Line
                                OutFile.write(Line + "")
                        # clean file of unneeded char/data 
                        # Remove the line-ending characters



                        if LineNumber > 0:
                                Line=Line.strip('\n')
                                ElementList = Line.split('\t')

                                Total = ElementList[8]
                                percentEdit = float(ElementList[14])/float(ElementList[15])
                                #print '.....................'
                                #print 'total'
                                #print Total
                                #print 'percentEdit'
                                #print percentEdit

                                
                                if int(Total) >= int(TotalCountThreshold) and percentEdit >= float(PercentEditThreshold):
                                        #print ' ..........passes threshold '
                                        outputLine = "\t".join(ElementList)
                                        #print outputLine

                                #create new outputfile
                                        OutFile.write(outputLine + "\n")
                
                                #else:
					# print '..........doesnt pass threshold'
                        # Index the counter used to keep track of line numbers
                        LineNumber = LineNumber + 1
                        
                # After the loop is completed, close the file
                Infile.close()
                OutFile.close()
