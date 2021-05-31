#!/usr/bin/env python

Usage = """
convert_editsites_to_bedgraph.py edit_sites_A2G.txt
converts outoput from find_rnaeditsites.pl to bedgraph format. Creates a bed track, to display labelled edit sites (% and number of total reads

"""

import sys
# the first line of the 'arguement is always the name of the program
if len(sys.argv)<2:  
        print (Usage)
else:
        FileList = sys.argv[1:]
        for InfileName in FileList:
                print (InfileName)

                Infile = open (InfileName, 'r')
                OutfileName = InfileName + '.bedgraph'
                Outfile = open(OutfileName, 'w')
        
                #	Outfile.write('trackname = "' + InfileName+ '"' + ' description = "test description" visibility=2 color=0,0,0 useScore=0' + "\n")
                # Loop through each line in the file
                for Line in Infile:
                        #skip the header
                        if (Line.startswith("Chr\t")): 
                                continue
                        else:
                                Line=Line.strip('\n')
                                ElementList = Line.split('\t')
                                
                                #calc % editing and add to rows
                                G = float(ElementList [14])
                                Total = float(ElementList [15])
                                percentG = (float(G/Total))*100
                                #round the editing percentage to 1 decimal place
                                percentG=round(percentG,1)
                                # print 'percent G: ' + str(percentG)
                                # create name in format X%_xr
                                name = str(percentG) + '%_' + str(ElementList[15]) + 'r' 
                                # print name 
                                # create uniq ID
                                ID = ElementList[0] + '_' + ElementList[1]
                                # append to end of list
                                ElementList.append(ID)
                                
                                strList = "\t".join(ElementList)
                                # print strList
                                
                                # create bed format
                                for_bedgraph = str(ElementList[0]) + "\t" + str(ElementList[1]) + "\t" + str(ElementList[1]) + "\t" + str(percentG) + "\t" + str(name) + "\t" + strList
                                Outfile.write(for_bedgraph + "\n")
                                
                # After the loop is completed, close the file
                Infile.close()
                Outfile.close()
                            
                print("Done, bed file created")
