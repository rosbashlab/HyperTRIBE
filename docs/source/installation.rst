Installation
============

Installation of HyperTRIBE's computational pipeline involves installing a set of softwares, downloading a set of annotation files, and updating the shells scripts to provide their locations.


Software Dependencies (tested version)
--------------------------------------
- Trimmomatics (v. 0.36)
- Bowtie2 (v. 2.1.0, 2.2.9)
- STAR (v. 2.5.2b)
- samtools (v. 1.3.1)
- bedtools suite (v. 2.16.2)
- Perl (5.8.8, 5.12.5, 5.22.1) 
- Perl module DBI.pm (1.631, 1.636) and DBD:mysql (4.042)
- MySQL database (MySQL, MariDB)
- Python (3.0 or above)
- SRA Toolkit 

Operating systems used: RHEL 5.11 and RHEL 7.2. HyperTRIBE is likely to work seamlessly with the latest version of these software/packages as well.

Source Code
-----------
Download the source code from github.
::

   cd /path_from_root/desired_location
   git clone https://github.com/rosbashlab/HyperTRIBE

Install MySQL (MariaDB) database and Perl modules
-------------------------------------------------
Install MariaDB, here is a helpful `resource <http://hypertribe.readthedocs.io/en/latest/mariadb.html>`_. Once the root password is set, log in to MySQL using the root password. 
::

    mysql -h localhost -u root -p

If MySQL is installed on a different server, then use the IP address of the server to replace “localhost”.

HyperTRIBE software requires the creation of MySQL database. The scripts *load_matrix_data.pl* and *find_rnaeditsites.pl* uses a MySQL database called **dmseq** by default. Feel free to create a different MySQL database for this purpose, but in that case please update the "$database" variable on line 24 in each perl scripts.
::

    CREATE DATABASE dmseq;   

After exiting MySQL, install MariaDB-devel
::

    yum install mariadb-devel

Install the perl modules using cpan or cpanm unless they are already installed.
::

    #install DBI.pm
    cpan DBI
    #install DBD::mysql
    cpan DBD::mysql    


Update Perl Scripts
-------------------
Perl Script variables have to updated so that they can communicate with the MySQL database tables. Edit "load_matrix_data.pl" (line 22-25) and "find_rnaeditsites.pl" (line 23-26) to provide: 1. Provide MySQL username (in this case “root”); 2. Password (password for databases); 3. Name of database, which was created earlier. If the MySQL database is hosted on a different machine, change the “$host” variable from “localhost” to its IP address.


Download Annotations and Genome Sequences
-----------------------------------------
HyperTRIBE needs genome sequence in fasta format and transcriptome annotation in gtf format and ucsc table browser native format. Users need to download these files based on the organism and genome build of interest.The two annotation files are used at different step of the pipeline.

**Download RefSeq annotation for Drosophila** (dm6) from `UCSC Genome Browser <https://genome.ucsc.edu/index.html>`_ Tools=TableBrowser; clade=Insect; genome=D.melanogaster; assembly=dm6; group:Genes and Gene Predictions; track=RefSeq Genes; table=refFlat; output format=all field from selected tables; output file: choose a filename (dm6_refFlat.txt).

**Optional: Recommended procedure for getting the refseq annotation file when the Refseq Genes track is not available as an option in Table Browser.** Here is an example of his alternate approach downloading the refseq annotation for human genome from `UCSC Genome Browser <https://genome.ucsc.edu/index.html>`_ Tools:TableBrowser; clade=Mammal; genome=human; assembly=hg38; group:Genes and Gene Predictions; track=NCBI RefSeq; table=RefSeq Curated; output format=all field from selected tables; output file: hg38_refseq.txt. Since this annotation format is different from the refFlat format available for fruit fly, the columns have to be rearranged using awk.
::

    awk '{print $13"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' hg38_refseq.txt > hg38_ncbi_refseq_curated.txt


**Download RefSeq annotation in GTF format** from Illumina’s `Igenome page <https://support.illumina.com/sequencing/sequencing_software/igenome.html>`_ . As an laternative, this annotation may be downloaded from UCSC browser using the instruction above, except a small change, choose output file format as *GTF – gene transfer format*.

As part of HyperTRIBE code distribution these annotation files are provided as an example of required annotation file.


These `annotation files  <https://github.com/rosbashlab/HyperTRIBE/tree/master/annotation>`_ for Drosophila (dm6) are automatically downloaded when  HyperTRIBE source code is cloned. Uncompress the annotation files, which creates a directory with all the annotation files.
::

    cd /path_from_root/HyperTRIBE/annotation
    #uncompress the annotaion files as needed
    gunzip genes.gtf.gz

Download the genome sequence from UCSC genome browser. Create bowtie2 indices for the genome and STAR indices for the transcriptome. The *genes.gtf* file is used during the creation of STAR indices. 
::

    #Create bowtie2 indices.
    cd /location_of_genome/
    bowtie2-build genome.fa genome
    
Create STAR indices
::

     STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_indices --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
     #star indices are created at /location_of_genome/star_indices

Update Shell Scripts
--------------------
Update Shell Scripts to indicate location of software location, annotation and genome indices. Use a text editor, like nao or emacs, to open shell scripts and update the location of HyperTRIBE code, annotation files, Bowtie2 and STAR indices.

Edit these variables in **trim_and_align.sh**
::

    star_indices="/path_from_root/star_indices"
    TRIMMOMATIC_JAR="/path_from_root/trimmomatic.jar"
    PICARD_JAR="/path_from_root/picard.jar"

If you want to use a different trimmer or aligner, feel free to change the code

Edit these variables in **trim_and_align_gDNA.sh**
::

    #location of bowtie2 indices
    bowtie_indexes="/path_from_root/genome"
    TRIMMOMATIC_JAR="/path_from_root/trimmomatic.jar"


Edit **load_table.sh**
::

    #location of HyperTRIBE code
    HyperTRIBE_DIR="/path_from_root/HyperTRIBE/CODE"

**Congratulations!!! Now, you are ready to run HyperTRIBE.**



