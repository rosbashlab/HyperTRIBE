#!/usr/bin/env python

# Set the input file name
# (The program must be run from within the directory 
# that contains this data file)
Usage = """
convert_editsites_to_bedgraph.py
converts outoput from find_rnaeditsites.pl to bedgraph format. Creates a bed track, to display labelled edit sites (% and number of total reads

"""

import sys

if len(sys.argv)<2:  # the first line of the 'arguement is always the name of the program
        print Usage
else:
        FileList = sys.argv[1:]
        for InfileName in FileList:
                #InfileName = InfileName
                print InfileName

                Infile = open (InfileName, 'r')
        
                OutfileName = InfileName + '.bedgraph'
                Outfile = open(OutfileName, 'w')
        
        #       Initialize the counter used to keep track of line numbers
        #       no header so counter starts at one
                LineNumber = 0
        # Loop through each line in the file


	#	Outfile.write('trackname = "' + InfileName+ '"' + ' description = "test description" visibility=2 color=0,0,0 useScore=0' + "\n")
                for Line in Infile:
			if LineNumber >= 1:
				#print LineNumber
                # clean file of unneeded char/data 
                # Remove the line-ending characters
 	                      	Line=Line.strip('\n')
				#print Line
        	               	ElementList = Line.split('\t')

   				#print ElementList
                
                #calc % editing and add to rows
        	                G = float(ElementList [14])
              		        Total = float(ElementList [15])
                        	percentG = (float(G/Total))*100
                #print 'percent G: ' + str(percentG)
                #        	ElementList.insert(18,percentG)
		# create name in format X%_xr
				name = str(int(percentG)) + '%_' + str(ElementList[15]) + 'r' 
				#print name 
		# create uniq ID
				ID = ElementList[0] + '_' + ElementList[1]
				#print ID
				ElementList.insert(18, ID)

				#print LineNumber, ":", ElementList			
			
				#print str(ElementList).strip('')
				strList = "\t".join(str(e) for e in ElementList)
				#print strList

				
                #create bed format
                #        	forbed = str(ElementList[0]) + "\t" + str(ElementList[1]) + "\t" + str(ElementList[1]) + "\t" + str(ElementList[18]) + "\t" + str(ElementList[20]) + "\t" + strList
                                forbed = str(ElementList[0]) + "\t" + str(ElementList[1]) + "\t" + str(ElementList[1]) + "\t" + str(int(percentG)) + "\t" + str(name) + "\t" + strList
                        	Outfile.write(forbed + "\n")
	
                # Index the counter used to keep track of line numbers
                        LineNumber = LineNumber + 1
                
        # After the loop is completed, close the file
                Infile.close()
                Outfile.close()

		print "Done, bed file created"
