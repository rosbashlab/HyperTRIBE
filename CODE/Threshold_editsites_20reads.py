#!/usr/bin/env python

Usage = """

Used for Thresholding of find_rnaeditsites.pl output with read and edit % threshold 
Original script by Aoife Mcmahon: 11_July_2013
Updated by Reazur Rahman: 2019

Command line:  python Threshold_editsites_20reads.py [input_file1] [input_file2] [input_file3] [etc] &

Removes 'edited sites' below threshold
1) remove entries <20 reads  (column 8 - Total)
2) remove entries < 0.10 edited (edit/total) (column 14/column 15) (editbasecount/Total)


Appends '.threshold' to the output files

To change thresholds change these variables in the script:
                TotalCountThreshold = 20
                PercentEditThreshold = 0.10 (10%)

"""

import sys

# These are the Thresholds, entries below these values will be removed.  Edit here if necessary
TotalCountThreshold = 20
PercentEditThreshold = 0.1

if len(sys.argv)<2:  # the first line of the arguement is always the name of the program
        print(Usage)
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
                        # print the header line, go the next element of the loop
                        if (Line.startswith("Chr\t")):
                                OutFile.write(Line)
                                continue

                        Line=Line.strip('\n')
                        ElementList = Line.split('\t')
                        
                        Total = ElementList[8]
                        percentEdit = float(ElementList[14])/float(ElementList[15])
                                        
                        if int(Total) >= int(TotalCountThreshold) and percentEdit >= float(PercentEditThreshold):
                                #print ' ..........passes threshold '
                                outputLine = "\t".join(ElementList)
                                OutFile.write(outputLine + "\n")
                
                
                # After the loop is completed, close the file
                Infile.close()
                OutFile.close()
