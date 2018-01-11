Example Run
===========
This page provides documentation for an example run of with sample dataset based on alignment to *Drosophila* (dm6). Load the dataset into mysql table and identify RNA editing sites. This assumes that load_matrix_data.pl and find_rnaeditsites.pl has been setup properly. 
::

    #path of the examples files
    cd /path_from_root/HyperTRIBE/examples/
    #unzip the compressed files
    gunzip *.gz
    # prepare the annotation files for drosophila dm6
    cd /path_from_root/HyperTRIBE/annotations/
    gunzip *.gz

Open *load_table.sh* to update the HyperTRIBE_DIR variable
::

    #path of HyperTRIBE code
    HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"

Open *rnaedit_wtRNA_RNA.sh* to update the HyperTRIBE_DIR variable
::

    #path of HyperTRIBE code
    HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"
    #set the location of the annotation file
    annotationfile="/path_from_root/HyperTRIBE/annotations/refFlat.txt

Now, load the data to MySQL and run the shell scripts. This script converts the SAM file to matrix and then loads the data into MySQL tables.
::

    #upload data to mysql tables, one for RNA and one gDNA
    #1. First argument alignment SAM file
    #2. mysql tablename
    #3. Expt name
    #4. Timepoint/Replicate Integer
    ./load_table.sh S2_wtRNA_chr2L.sort.sam exampleDB wtRNA 1
    ./load_table.sh HyperTRIBE_rep1_chr2L.sort.sam exampleDB HyperTRIBE 2
    # 

Now identify the editing sites with rnaedit_wtRNA_RNA.sh. If you changed the arguments of the loadPtable.sh script, please update the relevant variables in rnaedit_wtRNA_RNA.sh before running that script.
::

    ./rnaedit_wtRNA_RNA.sh

This should produce an output file called *HyperTRIBE_1_2_A2G.bedgraph*, where 3712 RNA editing sites are identified.  You can also log into mysql to see how the files, *S2_wtRNA_chr2L.sort.matrix* & *HyperTRIBE_rep1_chr2L.sort.matrix*, have populated the mysql tables.
