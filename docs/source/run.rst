Run
===

HyperTRIBE identifies RNA editing sites (A to I change, where I is read as G) by comparing RNA sequence from transcriptome with genomic DNA sequence (gDNA) or wild type RNA (wtRNA) from the same background. There are major steps in the pipeline:

1. Trim and align sequence libraries to Genome

2. Load alignments data to MySQL database

3. Find RNA edit sites using gDNA-RNA or wtRNA-RNA approaches

Optional: Download Test Datasets
--------------------------------
Demonstrate the workflow of HyperTRIBE computational analysis for both gDNA-RNA and wtRNA-RNA approaches using this dataset. Download a dataset of five sequencing libraries from NCBI GEO: GSE102814. Use SRA accession number of each library, then use SRA Toolkit to download the files in sra format, and then subsequently convert it to fastq format. The SRA accessions for the sequence libraries along with identifiers are: 1) S2 Genomic DNA: SRR3177714; 2) S2 WT mRNA: SRR6426146; 3) Hrp48 HyperTRIBE Replicate 1: SRR5944748; 4) Hrp48 HyperTRIBE Replicate 2: SRR5944749; 5) HyperADARcd alone: SRR5944750.
::

    # Download the SRA file:
    prefetch SRR3177714
    # Convert SRA to fastq format:
    fastq-dump SRR3177714.sra
    # This produces a fastq file called SRR3177714.fastq. Rename the file:
    mv SRR3177714.fastq S2_gDNA.fastq

Repeat these steps for the other files (mentioned above) to produce these fastq files: S2_wtRNA.fastq, HyperTRIBE_rep1.fastq, HyperTRIBE_rep2.fastq, and HyperADARcd_rep1.fastq 


1. Trim and align sequence libraries to Genome
----------------------------------------------
Reads are trimmed to remove low quality bases from genomic DNA and RNA libraries (fastq files). 6 nucleotide from 5' end of the read is removed to reduce potential errors from random hexamer mispriming. Align to reference genome or transcriptome, and they are sorted by using samtools. Run the “trim_and_align_gDNA.sh” shell script for gDNA library.

A. Trim and align gDNA Libraries
::

    nohup /path_from_root/HyperTRIBE/CODE/trim_and_align_gDNA.sh S2_gDNA.fastq &

The output: "S2_gDNA.sort.sam” records the alignment to the genome in SAM format. Similarly, run the “trim_and_align.sh” shell script for each RNA library

B. RNA Libraries
::

    nohup /path_from_root/HyperTRIBE/CODE/trim_and_align.sh S2_wtRNA.fastq &

The output: “S2_wtRNA.sort.sam” records the alignment to the transcriptome in SAM format. These alignments are processed by Picard to remove PCR duplicates. Repeat this step for the other RNA libraries.

2. Load alignment data to MySQL database
----------------------------------------
Convert alignment in SAM format to a matrix format, which records nucleotide frequency at each position in the transcriptome/genome. This data is then uploaded to a MySQL table based on the arguments provided.

Run the “load_table.sh” script with four arguments: 1. SAM file name 2. MySQL table name 3. Experiment name (unique identifier for the sequence library; include alphabets and digits) 4. An integer reflecting replicate/time-point. The combination of MySQL table name, experiment name, and replicate/time-point integer must be unique for each library. 
::

    nohup /path_from_root/HyperTRIBE/CODE/load_table.sh sam_filename mysql_tablename expt_name replicate/timepoint &
    #1. sam_filename
    #2. mysql_tablename
    #3. expt_name (unique identifier for the experiment, include alphabets and digits). Create short names. The expt_name cannot be longer than 20 characters.
    #4. replicate or timepoint: This is has to be an integer

Example for gDNA library
::

    nohup /path_from_root /HyperTRIBE/CODE/load_table.sh S2_gDNA.sort.sam s2_gDNA s2_gDNA 25 &

Example for RNA library
::
    
    nohup /path_from_root/HyperTRIBE/CODE/load_table.sh S2_wtRNA.sort.sam testRNA rnalibs 2 &

Repeat for the other RNA libraries with appropriate arguments, for example: Hrp48_HyperTRIBE_rep1.sort.sam testRNA rnalibs 3; Hrp48_HyperTRIBE_rep2.sort.sam testRNA rnalibs 4; HyperADARcd_rep1.sort.sam testRNA rnalibs 5

**Carefully note the last three arguments for all the load_table.sh runs for each library because these values are needed in the next step of the analysis. The expt_name cannot be longer than 20 characters**

Optional: Confirm MySQL Tables contains the data
------------------------------------------------
Edit the perl script “diagnose_mysql_table.pl” so that it can interact with the MySQL database. Provide mysql username (in this case “root”), password and name of database which was created earlier (line 21-25)
::

    perl diagnose_mysql_table.pl –t testRNA wtRNA.matrix HyperTRIBE_rep1.matrix HyperTRIBE_rep1.matrix H yperADARcd_rep1.matrix

This script calculates the total number of lines in all the matrix files and compares it to total number of entries in the given MySQL table. These numbers should be equal, otherwise you might have data contamination. 

3.A. Find RNA edit sites using gDNA-RNA approaches
--------------------------------------------------
Find RNA editing sites by comparing the gDNA and RNA sequence in mysql tables. For each position in the genome, the nucleotide composition is looked up using the tablename, expt name and replicate/timepoint integer. 

Copy rnaedit_gDNA_RNA.sh to your working directory
::

    cd /directory_of_choice/
    cp /location_from_root/TRIBE/CODE/rnaedit_gDNA_RNA.sh .

Edit the following variables in *rnaedit_gDNA_RNA.sh* by using a text editor like “nano". Based on how we processed the data above, it should look like this (line 3, 8-13)
::

    #open and edit the following variables in the script as needed
    HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"
    annotationfile="/path_from_root/HyperTRIBE/annotation/dm6_refFlat.txt "
    gDNAtablename="s2_gDNA"
    gDNAexp="s2_gDNA"
    gDNAtp="25"
    RNAtablename="testRNA"
    RNAexp="rnalibs"
    timepoint=(2 3 4 5)
    #the timepoint array allows you run multiple libraries one after another, if desired

Now, run the updated shell script from current directory
::

    ./rnaedit_gDNA_RNA.sh

First "find_rnaeditsites.pl" is used conduct a comparison of gDNA against RNA for each nucleotide in the transcriptome to call a set of editing sites (minimum coverage of nucleotide in reference table is hard coded to be 9 nucleotide). Then, “filter_by_threshold_without_header.pl” ensures that the editing sites have at least 5% editing and at least a coverage of 20 reads. The output for this shell script is a list of editing sites in bedgraph format, which is generated for each pairwise comparison. In this case there will be four bedgraph files with editing sites for: 1) S2_wtRNA: rnalibs_25_2_A2G.bedgraph; 2) HyperTRIBE_rep1: rnalibs_25_3_A2G.bedgraph; 3) HyperTRIBE_rep2: rnalibs_25_4_A2G.bedgraph; and 4) HyperADARcd_rep1: rnalibs_25_5_A2G.bedgraph



3.B. Find RNA edit sites using wtRNA-RNA approaches
---------------------------------------------------
Find RNA editing sites by using the wtRNA-RNA approach as an alternative to previous step. 

Copy rnaedit_wtRNA_RNA.sh to your working directory
::

    cd /directory_of_choice/
    cp /path_from_root/HyperTRIBE/CODE/rnaedit_wtRNA_RNA.sh .

Edit the following variables in *rnaedit_wtRNA_RNA.sh* by using a text editor like “nano". Based on how we processed the data above, it should look like this (line 3, 8-13)
::

    #open and edit the following variables in the script as needed
    HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"
    annotationfile="/path_from_root/HyperTRIBE/annotation/dm6_refFlat.txt "
    wtRNAtablename=" testRNA "
    wtRNAexp="rnalibs"
    wtRNAtp="2"
    RNAtablename="testRNA"
    RNAexp="rnalibs"
    timepoint=(3 4 5)
    #the timepoint array allows you run multiple libraries one after another, if desired

Now, run the updated shell script from current directory
::

    ./rnaedit_wtRNA_RNA.sh

First "find_rnaeditsites.pl" is used conduct a comparison of wtRNA against RNA for each nucleotide in the transcriptome to call a set of editing sites. Then, “filter_by_threshold_without_header.pl” ensures that the editing sites have at least 5% editing and at least a coverage of 20 reads. The output for this shell script is a list of editing sites in bedgraph format, which is generated for each pairwise comparison. In this case there will be three bedgraph files with editing sites for: 1) HyperTRIBE_rep1: rnalibs_2_3_A2G.bedgraph; 2) HyperTRIBE_rep2: rnalibs_2_4_A2G.bedgraph; and 3) HyperADARcd_rep1: rnalibs_2_5_A2G.bedgraph.


4. Post-processing of editing outputsOutputs
--------------------------------------------
For gDNA-RNA approach, identify high confidence set of HyperTRIBE editing sites. Use bedtools intersect to find the overlap between two HyperTRIBE replicates
::

    bedtools intersect -wa -wb -f 0.9 -r -a rnalibs_25_3_A2G.bedgraph -b rnalibs_25_4_A2G.bedgraph > present_both.bedgraph
    #Remove background (S2 wtRNA) editing sites:
    bedtools intersect -wa -v -f 0.9 -r -a present_both.bedgraph -b rnalibs_25_2_A2G.bedgraph > temp.bed
    #Remove HyperADARcd editing sites:
    bedtools intersect -wa -v -f 0.9 -r -a temp.bed -b rnalibs_25_5_A2G.bedgraph > HyperTRIBE_1_2_gDNA.bedgraph


For wtRNA-RNA approach, identify high confidence set of HyperTRIBE editing sites. Use bedtools to find the overlap between two HyperTRIBE replicates
::

    bedtools intersect -wa -wb -f 0.9 -r -a rnalibs_2_3_A2G.bedgraph -b rnalibs_2_4_A2G.bedgraph > present_both_wtRNA.bedgraph
    #Remove HyperADARcd editing sites:
    bedtools intersect -wa -v -f 0.9 -r -a present_both_wtRNA.bedgraph -b rnalibs_2_5_A2G.bedgraph > HyperTRIBE_1_2_wtRNA.bedgraph

The editing sites in bedgraph format can be visualized on IGV.

5. Column descriptions for bedgraph output files
------------------------------------------------
Description of column header in the bedgraph files are provided below: 
1. Chr name
2. Start coordinate
3. End coordinate
4. Editing percentage
5. Concatenation of editing percentage and reads (Total reads in HyperTRIBE RNA)
6. Chr name
7. Edit Coordinate
8. Name
9. Type
10. A count
11. T count
12. C count
13. G count
14. Total nucleotide count
15. A count from gDNA/wtRNA
16. T count from gDNA/wtRNA
17. C count from gDNA/wtRNA
18. G count from gDNA/wtRNA
19. Total count from gDNA/wtRNA
20. Editbase Count
21. Total nucleotide count
22. Editbase count from gDNA/wtRNA
23. Total nucleotide count from gDNA/wtRNA
24. Identifier (chr name and coordinate)

6. Create list of Edited Transcripts
------------------------------------
Summarize editing results by creating a list of transcripts that are marked by editing. Create gene list and summary for gDNA-RNA approach:
::

    perl /path_from_root/HyperTRIBE/CODE/summarize_results.pl HyperTRIBE_1_2_gDNA.bedgraph > HyperTRIBE_results_gDNA.xls

Create gene list and summary for wtRNA-RNA approach:
::

    perl /path_from_root/HyperTRIBE/CODE/summarize_results.pl HyperTRIBE_1_2_wtRNA.bedgraph > HyperTRIBE_results_wtRNA.xls

Description of the column headers for output file.
1) gene name
2) number of editing sites for gene
3) Avg editing percentage
4) edit_percentage_read concatenation for each replicate, separated by "," (13%_29r,13%_29r). ";" is used to separate between different editing sites
5) gene feature concatenation, separated by comma (EXON,INTRON). ";" is used to separate between different editing sites
6) Edit identifier concatenation, separated by comma for each editing site and  ";" is used to separate between different editing sites
